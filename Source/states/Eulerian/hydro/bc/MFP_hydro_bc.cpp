#ifdef AMREX_USE_EB

    #include "MFP_hydro_bc.H"

    #include "MFP_field.H"
    #include "MFP_hydro.H"
    #include "MFP_transforms.H"

//-----------------------------------------------------------------------------

HydroBoundaryEB::HydroBoundaryEB() {}
HydroBoundaryEB::~HydroBoundaryEB() {}

//-----------------------------------------------------------------------------

void HydroState::set_eb_bc(const sol::table& bc_def)
{
    BL_PROFILE("HydroState::set_eb_bc");

    std::string bc_type = bc_def.get<std::string>("type");

    if (bc_type == HydroSlipWall::tag) {
        eb_bcs.push_back(
          std::unique_ptr<HydroBoundaryEB>(new HydroSlipWall(global_idx, flux_solver.get())));
    } else if (bc_type == HydroNoSlipWall::tag) {
        if (!viscous) {
            Abort("Requested EB bc of type '" + bc_type +
                  "' without defining 'viscosity' for state '" + name + "'");
        }
        eb_bcs.push_back(std::unique_ptr<HydroBoundaryEB>(
          new HydroNoSlipWall(global_idx, flux_solver.get(), viscous.get(), bc_def)));
    } else if (bc_type == DirichletWall::tag) {
        eb_bcs.push_back(std::unique_ptr<HydroBoundaryEB>(
          new DirichletWall(global_idx, flux_solver.get(), bc_def)));
    } else if (bc_type == MultiStateWall::tag) {
        const std::string& em_name = bc_def["em_state"].get_or<std::string>("");

        if (em_name.empty()) {
            Abort("Requested EB bc of type '" + bc_type +
                  "' without defining 'em_state' for state '" + name + "'");
        }

        const State& em_state = get_state(em_name);

        if (em_state.get_type() != StateType::Field) {
            Abort("Requested EB bc of type '" + bc_type + "' for state '" + name + "', but '" +
                  em_state.name + "' is not of type 'StateType::Field'");
        }

        eb_bcs.push_back(std::unique_ptr<HydroBoundaryEB>(
          new MultiStateWall(*this, em_state, flux_solver.get(), bc_def)));
    } else {
        Abort("Requested EB bc of type '" + bc_type + "' which is not compatible with state '" +
              name + "'");
    }
}

//-----------------------------------------------------------------------------

std::string DirichletWall::tag = "dirichlet";

DirichletWall::DirichletWall() {}
DirichletWall::~DirichletWall() {}

DirichletWall::DirichletWall(int idx, RiemannSolver* flux, const sol::table& bc_def)
{
    BL_PROFILE("DirichletWall::DirichletWall");

    state_idx = idx;
    flux_solver = flux;

    // grab the wall state from the lua definition
    // only get stuff that is defined, anything else will
    // be filled in from the fluid side of the wall
    int i = 0;
    for (const auto& name : HydroState::prim_names) {
        sol::object val = bc_def[name];
        if (val.valid()) { wall_value.push_back({i, val.as<Real>()}); }
        ++i;
    }

    normal_flux.resize(flux->get_n_flux());
    cell_state.resize(+HydroDef::PrimIdx::NUM);
    wall_state.resize(+HydroDef::PrimIdx::NUM);
}

void DirichletWall::solve(Array<Array<Real, 3>, 3>& wall_coord,
                          Array<Real, AMREX_SPACEDIM> wall_centre,
                          const Vector<Array4<const Real>>& all_prim,
                          const int i,
                          const int j,
                          const int k,
                          const Real* dx,
                          Array<Vector<Real>, AMREX_SPACEDIM>& F)
{
    BL_PROFILE("DirichletWall::solve");

    //
    // get the inviscid flux
    //

    const Array4<const Real>& p4 = all_prim[state_idx];

    // grab the values we need
    for (size_t n = 0; n < +HydroDef::PrimIdx::NUM; ++n) { cell_state[n] = p4(i, j, k, n); }

    transform_global2local(cell_state, wall_coord, HydroState::prim_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    std::copy(cell_state.begin(), cell_state.end(), wall_state.begin());

    for (const auto& pair : wall_value) { wall_state[pair.first] = pair.second; }

    flux_solver->solve(cell_state, wall_state, normal_flux, nullptr);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, HydroState::cons_vector_idx);

    // split it up into components
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int n = 0; n < normal_flux.size(); ++n) {
            F[d][n] = wall_coord[0][d] * normal_flux[n];
        }
    }
}

//-----------------------------------------------------------------------------

std::string HydroSlipWall::tag = "slip_wall";

HydroSlipWall::HydroSlipWall() {}
HydroSlipWall::~HydroSlipWall() {}

HydroSlipWall::HydroSlipWall(int idx, RiemannSolver* flux)
{
    BL_PROFILE("HydroSlipWall::HydroSlipWall");

    state_idx = idx;
    flux_solver = flux;

    normal_flux.resize(flux->get_n_flux());
    cell_state.resize(+HydroDef::PrimIdx::NUM);
    wall_state.resize(+HydroDef::PrimIdx::NUM);
}

void HydroSlipWall::solve(Array<Array<Real, 3>, 3>& wall_coord,
                          Array<Real, AMREX_SPACEDIM> wall_centre,
                          const Vector<Array4<const Real>>& all_prim,
                          const int i,
                          const int j,
                          const int k,
                          const Real* dx,
                          Array<Vector<Real>, AMREX_SPACEDIM>& F)
{
    BL_PROFILE("HydroSlipWall::solve");

    //
    // get the inviscid flux
    //

    const Array4<const Real>& p4 = all_prim[state_idx];

    // grab the values we need
    for (size_t n = 0; n < +HydroDef::PrimIdx::NUM; ++n) { cell_state[n] = p4(i, j, k, n); }

    transform_global2local(cell_state, wall_coord, HydroState::prim_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    std::copy(cell_state.begin(), cell_state.end(), wall_state.begin());

    wall_state[+HydroDef::PrimIdx::Xvel] *= -1;

    Real shk = 1.0;  // consider there to be a shock present if we have a switching flux
    flux_solver->solve(cell_state, wall_state, normal_flux, &shk);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, HydroState::cons_vector_idx);

    // split it up into components
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int n = 0; n < +HydroDef::ConsIdx::NUM; ++n) {
            F[d][n] = wall_coord[0][d] * normal_flux[n];
        }
    }

    return;
}

//-----------------------------------------------------------------------------

std::string HydroNoSlipWall::tag = "no_slip_wall";

HydroNoSlipWall::HydroNoSlipWall() {}
HydroNoSlipWall::~HydroNoSlipWall() {}

HydroNoSlipWall::HydroNoSlipWall(int idx,
                                 RiemannSolver* flux,
                                 HydroViscous* visc,
                                 const sol::table& bc_def)
{
    BL_PROFILE("HydroNoSlipWall::HydroNoSlipWall");

    state_idx = idx;
    flux_solver = flux;
    viscous = visc;

    // grab the wall velocity from the lua definition
    wall_velocity[0] = 0.0;
    wall_velocity[1] = bc_def["v1"].get_or(0.0);
    wall_velocity[2] = bc_def["v2"].get_or(0.0);

    // grab a wall temperature
    wall_temp = bc_def["T"].get_or(-1.0);

    normal_flux.resize(flux->get_n_flux());
    cell_state.resize(+HydroDef::PrimIdx::NUM);
    wall_state.resize(+HydroDef::PrimIdx::NUM);
}

void HydroNoSlipWall::solve(Array<Array<Real, 3>, 3>& wall_coord,
                            Array<Real, AMREX_SPACEDIM> wall_centre,
                            const Vector<Array4<const Real>>& all_prim,
                            const int i,
                            const int j,
                            const int k,
                            const Real* dx,
                            Array<Vector<Real>, AMREX_SPACEDIM>& F)
{
    BL_PROFILE("HydroNoSlipWall::solve");

    //
    // get the inviscid flux
    //
    const Array4<const Real>& p4 = all_prim[state_idx];

    // grab the values we need
    for (size_t n = 0; n < +HydroDef::PrimIdx::NUM; ++n) { cell_state[n] = p4(i, j, k, n); }

    transform_global2local(cell_state, wall_coord, HydroState::prim_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    std::copy(cell_state.begin(), cell_state.end(), wall_state.begin());

    wall_state[+HydroDef::PrimIdx::Xvel] *= -1;

    Real shk = 1.0;  // consider there to be a shock present if we have a switching flux
    flux_solver->solve(cell_state, wall_state, normal_flux, &shk);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, HydroState::cons_vector_idx);

    // split it up into components
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int n = 0; n < normal_flux.size(); ++n) {
            F[d][n] = wall_coord[0][d] * normal_flux[n];
        }
    }

    //
    // get the viscous flux
    //

    if (viscous) {
        // get the wall state
        Array<Real, 3> global_wall_velocity = wall_velocity;
        Array<int, 1> idx_list = {0};
        // convert to global coords
        transform_local2global(global_wall_velocity, wall_coord, idx_list);

        Array<Real, 3>& wall_normal = wall_coord[0];

        const Array<int, 3> slope_idx = {+HydroDef::PrimIdx::Xvel,
                                         +HydroDef::PrimIdx::Yvel,
                                         +HydroDef::PrimIdx::Zvel};

        Array<Real, 3> velocity_slope = EulerianState::calc_wall_normal_slopes(slope_idx,
                                                                               global_wall_velocity,
                                                                               wall_normal,
                                                                               wall_centre,
                                                                               p4,
                                                                               i,
                                                                               j,
                                                                               k);

        // convert cell_state (flux vector form) to primitive vector form

        Real T, mu, kappa;
        viscous->get_coeffs(cell_state, T, mu, kappa);

        Array<Real, 3> dxinv = {0, 0, 0};
        AMREX_D_TERM(dxinv[0] = 1 / dx[0];, dxinv[1] = 1 / dx[1];, dxinv[2] = 1 / dx[2];)

        //
        // transform them to d/dx, d/dy and d/dz given transverse derivatives are zero
        Real dudx = velocity_slope[0] * wall_normal[0] * dxinv[0];
        Real dudy = velocity_slope[0] * wall_normal[1] * dxinv[1];
        Real dudz = velocity_slope[0] * wall_normal[2] * dxinv[2];
        //
        Real dvdx = velocity_slope[1] * wall_normal[0] * dxinv[0];
        Real dvdy = velocity_slope[1] * wall_normal[1] * dxinv[1];
        Real dvdz = velocity_slope[1] * wall_normal[2] * dxinv[2];
        //
        Real dwdx = velocity_slope[2] * wall_normal[0] * dxinv[0];
        Real dwdy = velocity_slope[2] * wall_normal[1] * dxinv[1];
        Real dwdz = velocity_slope[2] * wall_normal[2] * dxinv[2];

        Real divu = dudx + dvdy + dwdz;

        Real tautmp = 2.0 / 3.0 * mu * divu;

        Real tauxx = mu * 2 * dudx + tautmp;
        Real tauyy = mu * 2 * dvdy + tautmp;

        Real tauxy = mu * (dudy + dvdx);
        Real tauxz = mu * (dudz + dwdx);
        Real tauyz = mu * (dwdy + dvdz);

        Real xvel = wall_velocity[0] * wall_coord[0][0] + wall_velocity[1] * wall_coord[1][0] +
                    wall_velocity[2] * wall_coord[2][0];
        Real yvel = wall_velocity[0] * wall_coord[0][1] + wall_velocity[1] * wall_coord[1][1] +
                    wall_velocity[2] * wall_coord[2][1];
        Real zvel = wall_velocity[0] * wall_coord[0][2] + wall_velocity[1] * wall_coord[1][2] +
                    wall_velocity[2] * wall_coord[2][2];

        F[0][+HydroDef::ConsIdx::Xmom] -= tauxx;
        F[0][+HydroDef::ConsIdx::Ymom] -= tauxy;
        F[0][+HydroDef::ConsIdx::Zmom] -= tauxz;
        F[0][+HydroDef::ConsIdx::Eden] -= xvel * tauxx + yvel * tauxy + zvel * tauxz;

        F[1][+HydroDef::ConsIdx::Xmom] -= tauxy;
        F[1][+HydroDef::ConsIdx::Ymom] -= tauyy;
        F[1][+HydroDef::ConsIdx::Zmom] -= tauyz;
        F[1][+HydroDef::ConsIdx::Eden] -= xvel * tauxy + yvel * tauyy + zvel * tauyz;

    #if AMREX_SPACEDIM == 3
        Real tauzz = mu * 2 * dwdz + tautmp;
        F[2][+HydroDef::ConsIdx::Xmom] -= tauxz;
        F[2][+HydroDef::ConsIdx::Ymom] -= tauyz;
        F[2][+HydroDef::ConsIdx::Zmom] -= tauzz;
        F[2][+HydroDef::ConsIdx::Eden] -= xvel * tauxz + yvel * tauyz + zvel * tauzz;
    #endif

        // if we haven't specified a wall temperature
        // assume adiabatic (dTdn = 0)
        if (wall_temp > 0.0) {
            const Array<int, 1> temp_idx = {+HydroDef::PrimIdx::Temp};
            const Array<Real, 1> wall_temp_value = {wall_temp};

            Array<Real, 1> temp_slope = EulerianState::calc_wall_normal_slopes(temp_idx,
                                                                               wall_temp_value,
                                                                               wall_normal,
                                                                               wall_centre,
                                                                               p4,
                                                                               i,
                                                                               j,
                                                                               k);

            // calculate the temperature gradient
            Real dTdx = temp_slope[0] * wall_normal[0] * dxinv[0];
            Real dTdy = temp_slope[0] * wall_normal[1] * dxinv[1];

            F[0][+HydroDef::ConsIdx::Eden] -= kappa * dTdx;
            F[1][+HydroDef::ConsIdx::Eden] -= kappa * dTdy;

    #if AMREX_SPACEDIM == 3
            Real dTdz = temp_slope * wall_normal[2] * dxinv[2];
            F[2][+HydroDef::ConsIdx::Eden] -= kappa * dTdz;
    #endif
        }
    }
}

//-----------------------------------------------------------------------------

std::string MultiStateWall::tag = "multi_state_wall";

MultiStateWall::MultiStateWall() {}
MultiStateWall::~MultiStateWall() {}

MultiStateWall::MultiStateWall(const State& this_state,
                               const State& em_state,
                               RiemannSolver* flux,
                               const sol::table& bc_def)
{
    BL_PROFILE("MultiStateWall::MultiStateWall");

    // these idx get changed to data idx later
    state_idx = this_state.global_idx;
    em_state_idx = em_state.global_idx;
    flux_solver = flux;

    // you should grab any state information that you need here!

    normal_flux.resize(flux->get_n_flux());
    cell_hydro_state.resize(+HydroDef::PrimIdx::NUM);
    cell_em_state.resize(+FieldDef::ConsIdx::NUM);
    wall_state.resize(+HydroDef::PrimIdx::NUM);
}

void MultiStateWall::solve(Array<Array<Real, 3>, 3>& wall_coord,
                           Array<Real, AMREX_SPACEDIM> wall_centre,
                           const Vector<Array4<const Real>>& all_prim,
                           const int i,
                           const int j,
                           const int k,
                           const Real* dx,
                           Array<Vector<Real>, AMREX_SPACEDIM>& F)
{
    BL_PROFILE("MultiStateWall::solve");

    //
    // get the inviscid flux
    //

    // grab the em values we need and transform them to local coords
    const Array4<const Real>& em_p4 = all_prim[em_state_idx];
    for (size_t n = 0; n < +FieldDef::ConsIdx::NUM; ++n) { cell_em_state[n] = em_p4(i, j, k, n); }
    transform_global2local(cell_em_state, wall_coord, FieldState::vector_idx);

    // grab the hydro values we need and transform them to local coords
    const Array4<const Real>& hydro_p4 = all_prim[state_idx];
    for (size_t n = 0; n < +HydroDef::PrimIdx::NUM; ++n) {
        cell_hydro_state[n] = hydro_p4(i, j, k, n);
    }
    transform_global2local(cell_hydro_state, wall_coord, HydroState::prim_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    std::copy(cell_hydro_state.begin(), cell_hydro_state.end(), wall_state.begin());

    // HACK!
    // adjust the proportion of velocity reversal and hence flow into/out of the wall
    // based on the local x-electric field strength

    const Real local_Dx = cell_em_state[+FieldDef::ConsIdx::Dx];
    const Real local_Dy = cell_em_state[+FieldDef::ConsIdx::Dy];
    const Real local_Dz = cell_em_state[+FieldDef::ConsIdx::Dz];

    const Real local_D = sqrt(local_Dx * local_Dx + local_Dy * local_Dy + local_Dz * local_Dz);

    wall_state[+HydroDef::PrimIdx::Xvel] *= -1 * fabs(local_Dx) / local_D;

    Real shk = 1.0;  // consider there to be a shock present if we have a switching flux
    flux_solver->solve(cell_hydro_state, wall_state, normal_flux, &shk);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, HydroState::cons_vector_idx);

    // split it up into components
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int n = 0; n < +HydroDef::ConsIdx::NUM; ++n) {
            F[d][n] = wall_coord[0][d] * normal_flux[n];
        }
    }

    return;
}

//-----------------------------------------------------------------------------

#endif
