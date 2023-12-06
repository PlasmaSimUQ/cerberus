#ifndef MFP_HYDRORIEMANN_H
#define MFP_HYDRORIEMANN_H

#include "MFP_eulerian_riemann.H"
#include "MFP_factory.H"
#include "MFP_hydro_defs.H"
#include "MFP_utility.H"
#include "sol.hpp"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

using namespace amrex;

//================================================================================

class HydroRiemannSolver : public RiemannSolver
{
  public:
    HydroRiemannSolver();
    virtual ~HydroRiemannSolver();
};

template <typename D>
std::unique_ptr<HydroRiemannSolver> HydroRiemannSolverBuilder(const sol::table& def)
{
    if (def["flux"] == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

ClassFactory<HydroRiemannSolver>& GetHydroRiemannSolverFactory();

#endif  // MFP_FLUX_H
