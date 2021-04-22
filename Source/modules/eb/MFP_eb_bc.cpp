#ifdef AMREX_USE_EB
#include "MFP_eb_bc.H"
#include "MFP_transforms.H"
#include "AMReX_BLProfiler.H"

BoundaryEB::BoundaryEB(){};
BoundaryEB::~BoundaryEB(){};

const Vector<int> BoundaryEB::get_slopes() const
{
    return {};
}

Vector<Real> BoundaryEB::local2global(const Array<Array<Real,3>,3> &wall_coord) const
{
    return {};
}

const bool BoundaryEB::need_slopes() const
{
    return false;
}

//-----------------------------------------------------------------------------

std::string DirichletWall::tag = "dirichlet";

DirichletWall::DirichletWall(){}
DirichletWall::~DirichletWall(){}

DirichletWall::DirichletWall(RiemannSolver *flux,
                             const Vector<std::string> &names,
                             const Vector<int> &vec_idx,
                             const sol::table &bc_def)
{
    BL_PROFILE("DirichletWall::DirichletWall");
    flux_solver = flux;

    // grab the wall state from the lua definition
    // only get stuff that is defined, anything else will
    // be filled in from the fluid side of the wall
    int i = 0;
    for (const auto &name : names) {
        sol::object val = bc_def[name];
        if (val.valid()) {
            wall_value.push_back({i, val.as<Real>()});
        }
        ++i;
    }

    vector_idx = vec_idx;

}

void DirichletWall::solve(Array<Array<Real,3>,3> &wall_coord,
                          Vector<Real> &state,
                          Vector<Real> &normal_slope,
                          Array<Vector<Real>, AMREX_SPACEDIM> &F,
                          const Real *dx) const
{
    BL_PROFILE("DirichletWall::solve");
    //
    // get the inviscid flux
    //

    transform_global2local(state, wall_coord, vector_idx);

    // fabricate a state for inside the wall based on the provided state
    Vector<Real> W = state;

    for (const auto& pair : wall_value) {
        W[pair.first] = pair.second;
    }

    Vector<Real> normal_flux(flux_solver->get_flux_size());
    flux_solver->solve(state, W, normal_flux, nullptr);

    // convert back to global coordinate system
    transform_local2global(normal_flux, wall_coord, vector_idx);

    // split it up into components
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        for (int n=0; n<normal_flux.size(); ++n) {
            F[d][n] = wall_coord[0][d]*normal_flux[n];
        }
    }

}

//-----------------------------------------------------------------------------

#endif
