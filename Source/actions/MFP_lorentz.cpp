#include "MFP_lorentz.H"
#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"
#include "Eigen"
#include "Dense"

std::string Lorentz::tag = "Lorentz";
bool Lorentz::registered = GetActionFactory().Register(Lorentz::tag, ActionBuilder<Lorentz>);

Lorentz::Lorentz(){}
Lorentz::~Lorentz(){}

Lorentz::Lorentz(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    const sol::table state_names = def["states"];

    for (const auto& key_value_pair : state_names) {
        std::string state_name = key_value_pair.second.as<std::string>();
        State& istate = MFP::get_state(state_name);

        switch (istate.get_type()) {
        case State::StateType::Field: {
            if (field != nullptr) Abort("Only one field state can be set for the Plasma5 source "+name);
            field = static_cast<FieldState*>(&istate);
            state_indexes.push_back(istate.global_idx);
            field->associated_actions.push_back(action_idx);
            break;
        }
        case State::StateType::Hydro: {
            HydroState* hydro = static_cast<HydroState*>(&istate);
            species.push_back(hydro);
            state_indexes.push_back(istate.global_idx);
            hydro->associated_actions.push_back(action_idx);
            break;
        }
        default:
            Abort("An invalid state has been defined for the Plasma5 source "+name);
        }
    }

    return;
}

void Lorentz::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("Lorentz::get_data");

    Vector<Array<int,2>> options(species.size()+1);

    options[0] = {field->global_idx, 0};

    for (size_t i=0; i<species.size();++i) {
        options[i+1] = {species[i]->global_idx, 0};
    }

    Action::get_data(mfp, options, update, time);
}

void Lorentz::calc_time_derivative(MFP* mfp, Vector<UpdateData> &update, const Real time, const Real dt)
{
    BL_PROFILE("Lorentz::calc_time_derivative");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    size_t n_species = species.size();

    Vector<Vector<Real>> U(species.size());
    Vector<MultiFab*> species_data;
    for (size_t i=0; i<species.size();++i) {
        const HydroState& hstate = *species[i];
        U[i].resize(hstate.n_cons());
        update[hstate.data_idx].dU_status = UpdateData::Status::Changed;
    }

    Vector<Array4<Real>> species4(n_species);
    Vector<Array4<Real>> species_dU4(n_species);

    // define some 'registers'

    Real pD, pB;
    Real Bx, By, Bz;
    Real Ex, Ey, Ez;
    Real ep;

    Real q, m, r;
    Real rho, mx, my, mz;
    Real u, v, w;

    const Real Larmor = MFP::Larmor;
//    const Real Debye = MFP::Debye;
    const Real lightspeed = MFP::lightspeed;

//    const Real f1 = dt*Larmor/(Debye*Debye*lightspeed);

    // get charge and current density
    Real charge_density, current_x, current_y, current_z;

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(field->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif

        Array4<Real> const& field4 = update[field->data_idx].U.array(mfi);
//        Array4<Real> const& field_dU4 = update[field->data_idx].second.array(mfi);

        for (int n=0; n<n_species; ++n) {
            species4[n] = species_data[n]->array(mfi);
            species_dU4[n] = update[species[n]->data_idx].dU.array(mfi);
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


                    ep = field4(i,j,k,+FieldDef::ConsIdx::ep);

                    // magnetic field
                    Bx = field4(i,j,k,+FieldDef::ConsIdx::Bx);
                    By = field4(i,j,k,+FieldDef::ConsIdx::By);
                    Bz = field4(i,j,k,+FieldDef::ConsIdx::Bz);
                    pB = field4(i,j,k,+FieldDef::ConsIdx::psi);

                    // electric field
                    Ex = field4(i,j,k,+FieldDef::ConsIdx::Dx)/ep;
                    Ey = field4(i,j,k,+FieldDef::ConsIdx::Dy)/ep;
                    Ez = field4(i,j,k,+FieldDef::ConsIdx::Dz)/ep;
                    pD = field4(i,j,k,+FieldDef::ConsIdx::phi);

                    // get charge and current density
                    charge_density = 0.0;
                    current_x = 0.0;
                    current_y = 0.0;
                    current_z = 0.0;

                    for (size_t n = 0; n < n_species; ++n) {

                        const Array4<Real>& sp4 = species4[n];
                        Vector<Real>& UU = U[n];
                        for (size_t l=0; l<UU.size(); ++l) {
                            UU[l] = sp4(i,j,k,l);
                        }

                        rho =   sp4(i,j,k,+HydroDef::ConsIdx::Density);
                        mx =    sp4(i,j,k,+HydroDef::ConsIdx::Xmom);
                        my =    sp4(i,j,k,+HydroDef::ConsIdx::Ymom);
                        mz =    sp4(i,j,k,+HydroDef::ConsIdx::Zmom);


                        m = species[n]->get_mass_from_cons(UU);
                        q = species[n]->get_charge_from_cons(UU);

                        r = q/m;

                        charge_density += rho*r;
                        current_x += r*mx;
                        current_y += r*my;
                        current_z += r*mz;

                        u = mx/rho;
                        v = my/rho;
                        w = mz/rho;

                        species_dU4[n](i,j,k,+HydroDef::ConsIdx::Xmom) += dt*(rho*r/Larmor)*(lightspeed*Ex + v*Bz - w*By);
                        species_dU4[n](i,j,k,+HydroDef::ConsIdx::Ymom) += dt*(rho*r/Larmor)*(lightspeed*Ey + w*Bx - u*Bz);
                        species_dU4[n](i,j,k,+HydroDef::ConsIdx::Zmom) += dt*(rho*r/Larmor)*(lightspeed*Ez + u*By - v*Bx);
                        species_dU4[n](i,j,k,+HydroDef::ConsIdx::Eden) += dt*(rho*r*lightspeed)/Larmor*(u*Ex + v*Ey + w*Ez);
                    }
                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
