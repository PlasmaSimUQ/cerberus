#ifdef AMREX_USE_EB

#include "MFP_mhd_bc.H"
#include "MFP_mhd.H"
#include "MFP_transforms.H"

//-----------------------------------------------------------------------------

MHDBoundaryEB::MHDBoundaryEB(){}
MHDBoundaryEB::~MHDBoundaryEB(){}

//-----------------------------------------------------------------------------

std::string DirichletWallMHD::tag = "dirichlet";

DirichletWallMHD::DirichletWallMHD(){}
DirichletWallMHD::~DirichletWallMHD(){}

DirichletWallMHD::DirichletWallMHD(MHDRiemannSolver *flux,
                             const sol::table &bc_def)
{
    BL_PROFILE("DirichletWall::DirichletWall");
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
}

void DirichletWallMHD::solve(Array<Array<Real,3>,3> &wall_coord,
                          Array<Real,AMREX_SPACEDIM> wall_centre,
                          Array<Real,+MHDDef::PrimIdx::NUM> &cell_state,
                          Array4<const Real> const &prim4,
                          const int i, const int j, const int k, const Real *dx,
                          Array<Vector<Real>,AMREX_SPACEDIM> &F) const
{
    BL_PROFILE("DirichletWall::solve");
    //
    // get the inviscid flux
    //

    transform_global2local(cell_state, wall_coord, MHDState::prim_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    Array<Real,+MHDDef::PrimIdx::NUM> W = cell_state;

    for (const auto& pair : wall_value) {
        W[pair.first] = pair.second;
    }

    Array<Real,+MHDDef::ConsIdx::NUM> normal_flux;
    flux_solver->solve(cell_state, W, normal_flux, nullptr);

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

std::string MHDSlipWall::tag = "slip_wall";

MHDSlipWall::MHDSlipWall(){}
MHDSlipWall::~MHDSlipWall(){}

MHDSlipWall::MHDSlipWall(MHDRiemannSolver *flux)
{
    flux_solver = flux;
}

void MHDSlipWall::solve(Array<Array<Real,3>,3> &wall_coord,
                          Array<Real,AMREX_SPACEDIM> wall_centre,
                          Array<Real,+MHDDef::PrimIdx::NUM> &cell_state,
                          Array4<const Real> const &prim4,
                          const int i, const int j, const int k, const Real *dx,
                          Array<Vector<Real>,AMREX_SPACEDIM> &F) const
{
    BL_PROFILE("MHDSlipWall::solve");
    //
    // get the inviscid flux
    //

    transform_global2local(cell_state, wall_coord, MHDState::cons_vector_idx);

    // fabricate a state for inside the wall based on the provided state
    Array<Real,+MHDDef::PrimIdx::NUM> W = cell_state;

    W[+MHDDef::PrimIdx::Xvel] *= -1;

    Array<Real,+MHDDef::ConsIdx::NUM> normal_flux;

    Real shk = 1.0; // consider there to be a shock present if we have a switching flux
    flux_solver->solve(cell_state, W, normal_flux, &shk);

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
