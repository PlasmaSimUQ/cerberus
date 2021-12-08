#include "MFP_plasma5.H"
#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"
#include "Eigen"
#include "Dense"

std::string Plasma5::tag = "plasma5";
bool Plasma5::registered = GetActionFactory().Register(Plasma5::tag, ActionBuilder<Plasma5>);

Plasma5::Plasma5(){}
Plasma5::~Plasma5(){}

Plasma5::Plasma5(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    std::string solver_name = def["solver"];

    if (solver_name == "explicit") {
        solver = TimeIntegrator::ForwardsEuler;
    } else if (solver_name == "implicit") {
        solver = TimeIntegrator::BackwardsEuler;
    } else {
        Abort("Plasma5 source needs to know what type of solver to use, options are ['explicit','implicit']");
    }

    const sol::table state_names = def["states"];

    for (const auto& key_value_pair : state_names) {
        std::string state_name = key_value_pair.second.as<std::string>();
        State& istate = MFP::get_state(state_name);

        switch (istate.get_type()) {
        case State::StateType::Field: {
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

void Plasma5::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("Plasma5::get_data");

    Vector<Array<int,2>> options(species.size()+1);

    options[0] = {field->global_idx, 0};

    for (size_t i=0; i<species.size();++i) {
        options[i+1] = {species[i]->global_idx, 0};
    }

    Action::get_data(mfp, options, update, time);

}

void Plasma5::calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt)
{
    BL_PROFILE("Plasma5::solve");

    switch(solver) {
    case TimeIntegrator::ForwardsEuler :
        explicit_solve(mfp, update, time, dt);
        break;
    case TimeIntegrator::BackwardsEuler :
        implicit_solve(mfp, update, time, dt);
        break;
    default:
        Abort("How did we get here?");
    }

}

void Plasma5::explicit_solve(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt)
{
    BL_PROFILE("Plasma5::explicit_solve");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    size_t n_species = species.size();

    update[field->data_idx].dU_status = UpdateData::Status::Changed;

    Vector<Vector<Real>> U(species.size());
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
    const Real Debye = MFP::Debye;
    const Real lightspeed = MFP::lightspeed;

    const Real c_h = field->div_speed;
    const Real div_damp = field->div_damping;
    const Real c_d = div_damp != 0.0 ? c_h/div_damp : 0.0;

    const Real f1 = dt*Larmor/(Debye*Debye*lightspeed);
    const Real f2 = dt*c_h*c_h*f1/lightspeed;
    const Real f3 = c_d != 0.0 ? dt*c_h*c_h/(c_d*c_d) : 0.0;

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
        Array4<Real> const& field_dU4 = update[field->data_idx].dU.array(mfi);

        for (int n=0; n<n_species; ++n) {
            const int data_idx = species[n]->data_idx;
            species4[n] = update[data_idx].U.array(mfi);
            species_dU4[n] = update[data_idx].dU.array(mfi);
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


                        m = species[n]->gas->get_mass_from_cons(UU);
                        q = species[n]->gas->get_charge_from_cons(UU);

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

                    // electric field and divergence constraint sources

                    field_dU4(i,j,k,+FieldDef::ConsIdx::Dx) += -f1*current_x;
                    field_dU4(i,j,k,+FieldDef::ConsIdx::Dy) += -f1*current_y;
                    field_dU4(i,j,k,+FieldDef::ConsIdx::Dz) += -f1*current_z;

                    // source for D divergence correction
                    field_dU4(i,j,k,+FieldDef::ConsIdx::phi) += -field4(i,j,k,+FieldDef::ConsIdx::phi)*f3 + f2*charge_density;

                    // source for B divergence correction
                    field_dU4(i,j,k,+FieldDef::ConsIdx::psi) += -field4(i,j,k,+FieldDef::ConsIdx::psi)*f3;
                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}

void Plasma5::implicit_solve(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt)
{
    BL_PROFILE("Plasma5::implicit_solve");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);


    size_t n_species = species.size();

    update[field->data_idx].dU_status = UpdateData::Status::Changed;

    Vector<Vector<Real>> U(species.size());
    for (size_t i=0; i<species.size();++i) {
        const HydroState& hstate = *species[i];
        U[i].resize(hstate.n_cons());
        update[hstate.data_idx].dU_status = UpdateData::Status::Changed;
    }

    Vector<Array4<Real>> species4(n_species);
    Vector<Array4<Real>> species_dU4(n_species);


    // define the linear system matrix/vector and the associated solver
    int N = 3*(n_species+1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(N,N);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd c(N);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> cphqr;

    Real trace_B, trace_C, mc;

    Real pD, pB;
    Real Bx, By, Bz;
    Real Ex, Ey, Ez;
    Real ep;

    Real charge_density;
    Real q, m, r, g;
    Real rho, mx, my, mz;
    Vector<Real> R(n_species);

    const Real Larmor = MFP::Larmor;
    const Real Debye = MFP::Debye;
    const Real lightspeed = MFP::lightspeed;
    const Real c_h = field->div_speed;
    const Real div_damp = field->div_damping;
    const Real c_d = div_damp != 0.0 ? c_h/div_damp : 0.0;

    const Real f1 = dt*Larmor/(Debye*Debye*lightspeed);
    const Real f2 = dt*c_h*c_h*f1/lightspeed;
    const Real f3 = c_d != 0.0 ? dt*c_h*c_h/(c_d*c_d) : 0.0;

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
        Array4<Real> const& field_dU4 = update[field->data_idx].dU.array(mfi);

        for (int n=0; n<n_species; ++n) {
            const int data_idx = species[n]->data_idx;
            species4[n] = update[data_idx].U.array(mfi);
            species_dU4[n] = update[data_idx].dU.array(mfi);
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


                    // load in electric field
                    b(n_species*3 + 0) = Ex;
                    b(n_species*3 + 1) = Ey;
                    b(n_species*3 + 2) = Ez;

                    charge_density = 0.0;

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

                        m = species[n]->gas->get_mass_from_cons(UU);
                        q = species[n]->gas->get_charge_from_cons(UU);

                        r = q/m;
                        R[n] = r;

                        charge_density += rho*r;

                        charge_density += rho*r;
                        trace_B = -dt*r*rho*lightspeed/Larmor;
                        trace_C = dt*r*Larmor/(lightspeed*Debye*Debye);
                        mc = dt*r/Larmor;

                        // diagonal elements
                        A(3*n + 0, 3*n + 1) = -Bz*mc;
                        A(3*n + 0, 3*n + 2) =  By*mc;

                        A(3*n + 1, 3*n + 0) =  Bz*mc;
                        A(3*n + 1, 3*n + 2) = -Bx*mc;

                        A(3*n + 2, 3*n + 0) = -By*mc;
                        A(3*n + 2, 3*n + 1) =  Bx*mc;

                        // end column elements
                        A(3*n + 0, 3*n_species + 0) = trace_B;
                        A(3*n + 1, 3*n_species + 1) = trace_B;
                        A(3*n + 2, 3*n_species + 2) = trace_B;

                        // end row elements
                        A(3*n_species + 0, 3*n + 0) = trace_C;
                        A(3*n_species + 1, 3*n + 1) = trace_C;
                        A(3*n_species + 2, 3*n + 2) = trace_C;

                        // b vector (momentum components)
                        b(n*3 + 0) = mx;
                        b(n*3 + 1) = my;
                        b(n*3 + 2) = mz;

                    }

                    // solve linear system for updated momentum and electric field
                    cphqr.compute(A);
                    c = cphqr.solve(b);


                    // the updates for the D field
                    // note that the linear system has solved for the updated field but we want the delta value
                    // hence we calculate delta = new - old
                    field_dU4(i,j,k,+FieldDef::ConsIdx::Dx) += c(n_species*3 + 0)*ep - field4(i,j,k,+FieldDef::ConsIdx::Dx);
                    field_dU4(i,j,k,+FieldDef::ConsIdx::Dy) += c(n_species*3 + 1)*ep - field4(i,j,k,+FieldDef::ConsIdx::Dy);
                    field_dU4(i,j,k,+FieldDef::ConsIdx::Dz) += c(n_species*3 + 2)*ep - field4(i,j,k,+FieldDef::ConsIdx::Dz);

                    // new electric field
                    Ex = c(n_species*3 + 0);
                    Ey = c(n_species*3 + 1);
                    Ez = c(n_species*3 + 2);

                    // get the updates for hydro species
                    for (int n=0; n<n_species; ++n) {

                        Array4<Real>& sp4 = species4[n];

                        mx = c(n*3 + 0);
                        my = c(n*3 + 1);
                        mz = c(n*3 + 2);

                        species_dU4[n](i,j,k,+HydroDef::ConsIdx::Xmom) += mx - sp4(i,j,k,+HydroDef::ConsIdx::Xmom);
                        species_dU4[n](i,j,k,+HydroDef::ConsIdx::Ymom) += my - sp4(i,j,k,+HydroDef::ConsIdx::Ymom);
                        species_dU4[n](i,j,k,+HydroDef::ConsIdx::Zmom) += mz - sp4(i,j,k,+HydroDef::ConsIdx::Zmom);

                        species_dU4[n](i,j,k,+HydroDef::ConsIdx::Eden) += dt*(R[n]*lightspeed/Larmor)*(Ex*mx + Ey*my + Ez*mz);

                    }

                    // source for D divergence correction
                   field_dU4(i,j,k,+FieldDef::ConsIdx::phi) += -field4(i,j,k,+FieldDef::ConsIdx::phi)*f3 + f2*charge_density;

                   // source for B divergence correction
                   field_dU4(i,j,k,+FieldDef::ConsIdx::psi) += -field4(i,j,k,+FieldDef::ConsIdx::psi)*f3;
                }
            }
        }
        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}

Real Plasma5::get_allowed_time_step(MFP* mfp)
{

    if (solver != TimeIntegrator::ForwardsEuler) return std::numeric_limits<Real>::max();

    constexpr Real PI = acos(-1);

    const Real* dx = mfp->Geom().CellSize();

    BL_PROFILE("Plasma5::explicit_solve");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    MultiFab* field_data = &(mfp->get_new_data(field->data_idx));


    size_t n_species = species.size();

    Vector<Vector<Real>> U(species.size());
    Vector<MultiFab*> species_data;
    for (size_t i=0; i<species.size();++i) {
        const HydroState& hstate = *species[i];
        species_data.push_back(&(mfp->get_new_data(hstate.data_idx)));
        U[i].resize(hstate.n_cons());
    }

    Vector<Array4<Real>> species4(n_species);



    // define some 'registers'

    Real Bx, By, Bz;
    Real q, m, r;
    Real rho;


    const Real Larmor = MFP::Larmor;
    const Real Debye = MFP::Debye;
    Real D2 = Debye*Debye;

    Real max_speed = MFP::lightspeed;
    max_speed = std::max(max_speed, field->div_speed);

    Real dt = dx[0]/max_speed;
    for (int d=1; d<AMREX_SPACEDIM; ++d) {
        dt = std::min(dt, dx[d]/max_speed);
    }

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

        Array4<Real> const& field4 = field_data->array(mfi);

        for (int n=0; n<n_species; ++n) {
            species4[n] = species_data[n]->array(mfi);
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

                    // magnetic field
                    Bx = field4(i,j,k,+FieldDef::ConsIdx::Bx);
                    By = field4(i,j,k,+FieldDef::ConsIdx::By);
                    Bz = field4(i,j,k,+FieldDef::ConsIdx::Bz);

                    const Real B = std::sqrt(Bx*Bx + By*By + Bz*Bz);

                    for (size_t n = 0; n < n_species; ++n) {

                        const Array4<Real>& sp4 = species4[n];
                        Vector<Real>& UU = U[n];
                        for (size_t l=0; l<UU.size(); ++l) {
                            UU[l] = sp4(i,j,k,l);
                        }

                        rho =   sp4(i,j,k,+HydroDef::ConsIdx::Density);

                        m = species[n]->gas->get_mass_from_cons(UU);
                        q = species[n]->gas->get_charge_from_cons(UU);

                        r = q/m;


                        const Real omega_p = 10*std::sqrt(rho*r*r/D2)/(2*PI);

                        if (omega_p > 0) {
                            dt = std::min(dt, 1/omega_p);
                        }

                        const Real omega_c = 10*(std::abs(r)*B)/(Larmor*2*PI);

                        if (omega_c > 0) {
                            dt = std::min(dt, 1/omega_c);
                        }

                    }
                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }

    return dt;
}
