#include "MFP_collisions.H"
#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"

std::string Collisions::tag = "collisions";
bool Collisions::registered = GetActionFactory().Register(Collisions::tag, ActionBuilder<Collisions>);

Collisions::Collisions(){}
Collisions::~Collisions(){}

Collisions::Collisions(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    const sol::table state_names = def["states"];
    Vector<std::string> names;

    for (const auto& key_value_pair : state_names) {
        std::string state_name = key_value_pair.second.as<std::string>();
        State& istate = MFP::get_state(state_name);

        switch (istate.get_type()) {
        case State::StateType::Hydro: {
            HydroState* hydro = static_cast<HydroState*>(&istate);
            species.push_back(hydro);
            state_indexes.push_back(istate.global_idx);
            hydro->associated_actions.push_back(action_idx);
            names.push_back(state_name);
            break;
        }
        default:
            Abort("An invalid state has been defined for the Plasma5 source "+name);
        }
    }

    // grab any collision cross sections
    sol::table cross_sections = def["cross_sections"].get_or(sol::table());

    if (cross_sections.empty()) {

        // only need ccs with a neutral species
        bool has_neutral = true;
        for (const HydroState* s : species) {
            for (const Real& q : s->charge) {
                if (q != 0.0) {
                    has_neutral = false;
                    break;
                }
            }
            if (!has_neutral) break;
        }
        if (has_neutral)
            Abort("Error: collisions enabled with a neutral species without hydro_ccs being defined");
    }

    ccs.resize(species.size());
    for (int a=0; a<species.size(); ++a) {
        ccs[a].resize(species.size(),0);
    }

    // load the hydro cross sections
    std::pair<bool, int> index_a, index_b;
    Real sigma;
    for (const auto& a : cross_sections) {
        index_a = findInVector(names, a.first.as<std::string>());
        if (!index_a.first) {continue;}
        for (const auto& b : a.second.as<sol::table>()) {
            index_b = findInVector(names, b.first.as<std::string>());
            if (!index_b.first) {continue;}
            sigma = b.second.as<Real>()/(MFP::x_ref*MFP::x_ref); // non-dimensionalise
            ccs[index_a.second][index_b.second] = sigma;
            ccs[index_b.second][index_a.second] = sigma; // cross load as sigma_ab = sigma_ba
        }
    }


    return;
}

void Collisions::calc_time_derivative(MFP* mfp, Vector<std::pair<int, MultiFab> > &dU, const Real time, const Real dt)
{
    BL_PROFILE("Collisions::calc_time_derivative");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    size_t n_species = species.size();

    struct SpeciesData {
        Real rho, m, q, T, q2, n;
        Array<Real,3> vel;
    };

    Vector<Vector<Real>> U(n_species);
    Vector<SpeciesData> SD(n_species);
    Vector<MultiFab*> species_data;
    for (size_t i=0; i<n_species;++i) {
        const HydroState& hstate = *species[i];
        species_data.push_back(&(mfp->get_data(hstate.data_idx,time)));
        U[i].resize(hstate.n_cons());
    }

    Vector<Array4<Real>> species4(n_species);
    Vector<Array4<Real>> species_dU4(n_species);

    const Real Debye = MFP::Debye;
    const Real n0 = MFP::n0;

    // coefficients
    const Real c1 = 0.02116454531141366/(n0*pow(Debye,4)); // sqrt(2)/(12*pi^(3/2))*...
    const Real c2 = 0.4135669939329334; // (2/(9*pi))**(1/3.0)
    const Real c3 = 2.127692162140974*n0; // (4/3)*n0*sqrt(8/pi)
    const Real lnC = 10.0; // Coulomb logarithm


    for (int n=0; n<n_species; ++n) {
        dU[species[n]->data_idx].first = 1;
    }

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(species[0]->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif

        for (int n=0; n<n_species; ++n) {
            species4[n] = species_data[n]->array(mfi);
            species_dU4[n] = dU[species[n]->data_idx].second.array(mfi);
        }


        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) {
                        continue;
                    }
#endif

                    // calculate all of the data up-front
                    for (int n=0; n<n_species; ++n) {
                        SpeciesData& sd = SD[n];

                        const Array4<Real>& sp4 = species4[n];
                        Vector<Real>& Un = U[n];
                        for (size_t l=0; l<Un.size(); ++l) {
                            Un[l] = sp4(i,j,k,l);
                        }

                        sd.rho = Un[+HydroDef::ConsIdx::Density];
                        sd.m = species[n]->get_mass_from_cons(Un);
                        sd.q = species[n]->get_charge_from_cons(Un);
                        sd.T = species[n]->get_temperature_from_cons(Un);
                        sd.q2 = sd.q*sd.q;
                        sd.n = sd.rho/sd.m;

                        const Real i_rho = 1/sd.rho;
                        sd.vel[0] = Un[+HydroDef::ConsIdx::Xmom]*i_rho;
                        sd.vel[1] = Un[+HydroDef::ConsIdx::Ymom]*i_rho;
                        sd.vel[2] = Un[+HydroDef::ConsIdx::Zmom]*i_rho;

                    }

                    for (int a = 0; a < n_species; ++a) {

                        SpeciesData& sda = SD[a];

                        for (int b = a+1; b < n_species; ++b) {

                            SpeciesData& sdb = SD[b];

                            const Real m_ab = (sda.m*sdb.m)/(sda.m + sdb.m);


                            Real du2 = 0.0;
                            Array<Real,3> du;
                            for (int i=0; i<3; ++i) {
                                du[i] = sdb.vel[i] - sda.vel[i];
                                du2 += du[i]*du[i];
                            }

                            // collision frequency
                            Real nu;
                            if ((sda.q != 0) && (sdb.q != 0)) {
                                Real coeff_1 = c1*((sda.q2*sdb.q2*lnC)/(m_ab*sda.m));
                                nu = sdb.n*coeff_1*pow(c2*du2 + sda.T/sda.m + sdb.T/sdb.m, -1.5);
                            } else {
                                Real coeff_2 = (sdb.m*ccs[a][b])/(sda.m + sdb.m);
                                nu = sdb.n*coeff_2*c3*sqrt(sda.T/sda.m + sdb.T/sdb.m);
                            }

                            // collisions

                            Real Q = m_ab*sda.n*nu*du2 + 3*sda.rho*nu/(sda.m + sdb.m)*(sdb.T - sda.T);

                            Array<Real,3> R;
                            for (int n=0; n<3; ++n) {
                                R[n] = sda.rho*nu*du[n];
                                Q += R[n]*sda.vel[n];
                            }

                            for (int n=0; n<3; ++n) {
                                species_dU4[a](i,j,k,+HydroDef::ConsIdx::Xmom+n) += dt*R[n];
                                species_dU4[b](i,j,k,+HydroDef::ConsIdx::Xmom+n) -= dt*R[n];
                            }

                            species_dU4[a](i,j,k,+HydroDef::ConsIdx::Eden) += dt*Q;
                            species_dU4[b](i,j,k,+HydroDef::ConsIdx::Eden) -= dt*Q;
                        }
                    }
                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
