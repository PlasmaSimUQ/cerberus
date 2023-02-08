#ifdef AMREX_USE_EB

#include "MFP_mhd_bc.H"
#include "MFP_mhd.H"
#include "MFP_transforms.H"

//-----------------------------------------------------------------------------

MHDBoundaryEB::MHDBoundaryEB(){}
MHDBoundaryEB::~MHDBoundaryEB(){}

//-----------------------------------------------------------------------------

void MHDState::set_eb_bc(const sol::table &bc_def)
{

    std::string bc_type = bc_def.get<std::string>("type");

    MHDRiemannSolver* flux = static_cast<MHDRiemannSolver*>(flux_solver.get());

    if (bc_type == MHDSlipWall::tag) {
        eb_bcs.push_back(std::unique_ptr<MHDBoundaryEB>(new MHDSlipWall(global_idx, flux)));
    } else if (bc_type == DirichletWallMHD::tag) {
        eb_bcs.push_back(std::unique_ptr<MHDBoundaryEB>(new DirichletWallMHD(global_idx, flux, bc_def)));
    } else {
        Abort("Requested EB bc of type '" + bc_type + "' which is not compatible with state '" + name + "'");
    }
}

//-----------------------------------------------------------------------------

std::string DirichletWallMHD::tag = "dirichlet";

DirichletWallMHD::DirichletWallMHD(){}
DirichletWallMHD::~DirichletWallMHD(){}

DirichletWallMHD::DirichletWallMHD(int idx, MHDRiemannSolver *flux,
                             const sol::table &bc_def)
{
    BL_PROFILE("DirichletWall::DirichletWall");
    state_idx = idx;
    flux_solver = flux;

    // grab the wall state from the lua definition
    // only get stuff that is defined, anything else will
    // be filled in from the fluid side of the wall
    int i = 0;
    for (const auto &name : MHDState::prim_names) {
        sol::object val = bc_def[name];
        if (val.valid()) {
            wall_value.push_back({i, val.as<Real>()});
        }
        ++i;
    }

    cell_state.resize(+MHDDef::PrimIdx::NUM);
    wall_state.resize(+MHDDef::PrimIdx::NUM);
    normal_flux.resize(+MHDDef::ConsIdx::NUM);

}

void DirichletWallMHD::solve(Array<Array<Real,3>,3> &wall_coord,
                             Array<Real,AMREX_SPACEDIM> wall_centre,
                             const Vector<Array4<const Real>>& all_prim,
                             const int i, const int j, const int k, const Real *dx,
                             Array<Vector<Real>,AMREX_SPACEDIM> &F)
{
    BL_PROFILE("DirichletWall::solve");
    //
    // get the inviscid flux
    //

    const Array4<const Real>& p4 = all_prim[state_idx];

    // grab the values we need
    for (size_t n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
        cell_state[n] = p4(i,j,k,n);
    }

    transform_global2local(cell_state, wall_coord, MHDState::prim_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    std::copy(cell_state.begin(), cell_state.end(), wall_state.begin());

    for (const auto& pair : wall_value) {
        wall_state[pair.first] = pair.second;
    }

    flux_solver->solve(cell_state, wall_state, normal_flux, nullptr);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, MHDState::cons_vector_idx);

    // split it up into components
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        for (int n=0; n<normal_flux.size(); ++n) {
            F[d][n] = wall_coord[0][d]*normal_flux[n];
        }
    }

}

//-----------------------------------------------------------------------------

// TODO: handle magnetic field

std::string MHDSlipWall::tag = "slip_wall";

MHDSlipWall::MHDSlipWall(){}
MHDSlipWall::~MHDSlipWall(){}

MHDSlipWall::MHDSlipWall(int idx, MHDRiemannSolver *flux)
{
    state_idx = idx;
    flux_solver = flux;
    cell_state.resize(+MHDDef::PrimIdx::NUM);
    wall_state.resize(+MHDDef::PrimIdx::NUM);
    normal_flux.resize(+MHDDef::ConsIdx::NUM);
}

void MHDSlipWall::solve(Array<Array<Real,3>,3> &wall_coord,
                        Array<Real,AMREX_SPACEDIM> wall_centre,
                        const Vector<Array4<const Real>>& all_prim,
                        const int i, const int j, const int k, const Real *dx,
                        Array<Vector<Real>,AMREX_SPACEDIM> &F)
{
    BL_PROFILE("MHDSlipWall::solve");
    //
    // get the inviscid flux
    //

    const Array4<const Real>& p4 = all_prim[state_idx];

    // grab the values we need
    for (size_t n=0; n<+MHDDef::PrimIdx::NUM; ++n) {
        cell_state[n] = p4(i,j,k,n);
    }

    transform_global2local(cell_state, wall_coord, MHDState::prim_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    std::copy(cell_state.begin(), cell_state.end(), wall_state.begin());

    wall_state[+MHDDef::PrimIdx::Xvel] *= -1;

    Real shk = 1.0; // consider there to be a shock present if we have a switching flux

    flux_solver->solve(cell_state, wall_state, normal_flux, &shk);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, MHDState::cons_vector_idx);

    // split it up into components
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        for (int n=0; n<+MHDDef::ConsIdx::NUM; ++n) {
            F[d][n] = wall_coord[0][d]*normal_flux[n];
        }
    }
    return;
}

#endif
