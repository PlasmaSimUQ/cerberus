#include "MFP_collisions.H"
#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"

Collisions::Collisions(){}
Collisions::~Collisions(){}

Collisions::Collisions(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    init_species_info(def);

    const sol::table state_names = def["states"];
    Vector<std::string> names;

    for (const auto& key_value_pair : state_names) {
        std::string state_name = key_value_pair.second.as<std::string>();
        State& istate = MFP::get_state(state_name);

        switch (istate.get_type()) {
        case State::StateType::Hydro: {
            HydroState* hydro = static_cast<HydroState*>(&istate);
            states.push_back(hydro);
            state_indexes.push_back(istate.global_idx);
            hydro->associated_actions.push_back(action_idx);
            names.push_back(state_name);
            break;
        }
        default:
            Abort("An invalid state has been defined for the Collisions source "+name);
        }
    }

    size_t num_states = states.size();
    size_t n_species = species_info.size();


    std::pair<bool, int > found;
    for (int i=0; i<n_species; ++i) {
        std::string& sname = species_info[i].name;

        // search through the states
        for (int s=0; s<num_states; ++s) {
            found = findInVector(states[s]->gas->comp_names, sname);
            if (found.first) {
                SpeciesInfo& info = species_info[i];
                info.state_idx = s;
                info.alpha_idx = found.second;
                info.m = states[s]->gas->mass[info.alpha_idx];
                info.q = states[s]->gas->charge[info.alpha_idx];
                break;
            }
        }

        if (!found.first)
            Abort("Couldn't find sub-component for '"+sname+"'");
    }

    return;
}

void Collisions::init_species_info(const sol::table &def)
{
    const sol::table state_names = def["states"];
    sol::table species_names = def.get_or("species",sol::table());

    // assume that we have single species state
    if (!species_names.valid()) {
        species_names = state_names;
    }

    int idx = 0;
    for (const auto& item : species_names) {
        SpeciesInfo info;
        info.idx = idx; ++idx;
        info.name = item.second.as<std::string>();
        species_info.push_back(info);
    }

}


void Collisions::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("Collisions::get_data");

    Vector<Array<int,2>> options(states.size());

    for (size_t i=0; i<states.size();++i) {
        options[i] = {states[i]->global_idx, 0};
    }

    Action::get_data(mfp, options, update, time);

}


void Collisions::calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt)
{
    BL_PROFILE("Collisions::calc_time_derivative");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    size_t n_states = states.size();
    size_t n_species = species_info.size();

    for (int n=0; n<n_states; ++n) {
        update[states[n]->data_idx].dU_status = UpdateData::Status::Changed;
    }


    Vector<Vector<Real>> U(n_states), Q(n_states), alpha(n_states);
    Vector<Vector<Array<Real,+HydroDef::ConsIdx::NUM>>> accounting(n_states);
    Vector<MultiFab*> state_data;
    for (size_t i=0; i<n_states;++i) {
        const HydroState& hstate = *states[i];
        state_data.push_back(&update[hstate.data_idx].U);
        U[i].resize(hstate.n_cons());
        Q[i].resize(hstate.n_prim());
        accounting[i].resize(hstate.gas->n_species());
        alpha[i].resize(hstate.gas->n_species());
    }

    Vector<Array4<Real>> state4(n_species);
    Vector<Array4<Real>> state_dU4(n_species);


    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(states[0]->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif

        for (int n=0; n<n_states; ++n) {
            state4[n] = state_data[n]->array(mfi);
            state_dU4[n] = update[states[n]->data_idx].dU.array(mfi);
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

                    // get stuff from all states (U, Q, alpha)
                    for (int si=0; si<n_states; ++si) {

                        // get the vector of conserved properties
                        Vector<Real>& Un = U[si];
                        const Array4<const Real>& U4 = state4[si];
                        for (size_t l=0; l<Un.size(); ++l) {
                            Un[l] = U4(i,j,k,l);
                        }

                        if (Un[+HydroDef::ConsIdx::Density] < states[si]->effective_zero) {
                            std::fill(Q[si].begin(), Q[si].end(), 0.0);
                        } else {
                            states[si]->gas->cons2prim(U[si], Q[si]);
                            states[si]->gas->get_alpha_fractions_from_prim(Q[si], alpha[si]);
                        }
                    }



                    // calculate all of the data up-front
                    for (int n=0; n<n_species; ++n) {
                        SpeciesInfo& info = species_info[n];

                        // zero out the delta terms
                        info.delta.fill(0.0);

                        Vector<Real>& Un = U[n];
                        Vector<Real>& Qn = Q[n];
                        Real& a = alpha[n][info.alpha_idx];

                        if (Un[+HydroDef::ConsIdx::Density] > states[info.state_idx]->effective_zero) {
                            info.rho = a*Un[+HydroDef::ConsIdx::Density];
                            info.m = states[n]->gas->mass[info.alpha_idx];
                            info.q = states[n]->gas->charge[info.alpha_idx];
                            info.T = states[n]->gas->get_temperature_from_prim(Qn);
                            info.q2 = info.q*info.q;
                            info.n = info.rho/info.m;

                            const Real i_rho = 1/info.rho;
                            info.vel[0] = a*Un[+HydroDef::ConsIdx::Xmom]*i_rho;
                            info.vel[1] = a*Un[+HydroDef::ConsIdx::Ymom]*i_rho;
                            info.vel[2] = a*Un[+HydroDef::ConsIdx::Zmom]*i_rho;
                        } else {
                            info.rho = 0.0;
                            info.m = 0.0;
                            info.q = 0.0;
                            info.T = 0.0;
                            info.q2 = 0.0;
                            info.n = 0.0;
                            info.vel[0] = 0.0;
                            info.vel[1] = 0.0;
                            info.vel[2] = 0.0;
                        }
                    }

                    // calculate the mass, momentum, and energy transfer for each species
                    calc_update();

                    //////////////////////////////////////////////////////////////////
                    // now accumulate the changes across species and update the states

                    // get the original amount of stuff across all states
                    for (int si=0; si<n_states; ++si) {
                        for (int ai=0; ai<alpha[si].size(); ++ai) {
                            for (int ci=0; ci<+HydroDef::ConsIdx::NUM; ++ci) {
                                accounting[si][ai][ci] = U[si][ci]*alpha[si][ai]; // original amount of each conserved quantity
                            }
                        }
                    }

                    // get the updated amounts of stuff for all sub-components that have been updated by the reactor
                    for (int isp=0; isp<n_species; ++isp) {
                        SpeciesInfo& info = species_info[isp];
                        for (int ci=0; ci<+HydroDef::ConsIdx::NUM; ++ci) {
                            accounting[info.state_idx][info.alpha_idx][ci] += dt*info.delta[ci]; // updated amount
                        }
                    }

                    for (int si=0; si<n_states; ++si) {

                        // apply update for all conserved quantities
                        for (int ci=0; ci<+HydroDef::ConsIdx::NUM; ++ci) {

                            // accumulate how much stuff is held by all the sub-components in this state
                            Real sum = 0.0;
                            for (const auto& acc : accounting[si]) {
                                sum += acc[ci];
                            }

                            // now apply update
                            state_dU4[si](i,j,k,ci) += sum - state4[si](i,j,k,ci);
                        }

                        // apply the update for tracers
                        // note that we can simply use the density of the sub-component which we have already calculated
                        for (int ti=0; ti<states[si]->n_tracers; ++ti) {
                            state_dU4[si](i,j,k,+HydroDef::ConsIdx::NUM + ti) += accounting[si][ti][+HydroDef::ConsIdx::Density]  - state4[si](i,j,k,+HydroDef::ConsIdx::NUM + ti);
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
