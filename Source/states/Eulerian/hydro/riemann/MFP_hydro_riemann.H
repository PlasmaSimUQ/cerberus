#ifndef MFP_HYDRORIEMANN_H
#define MFP_HYDRORIEMANN_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "sol.hpp"
#include "MFP_factory.H"
#include "MFP_hydro_defs.H"
#include "MFP_utility.H"

using namespace amrex;

//================================================================================

class HydroRiemannSolver
{
public:
    HydroRiemannSolver();
    virtual ~HydroRiemannSolver();

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) {}

    virtual std::string get_tag() const {return "";}
    virtual bool requires_shock_detector() const {return false;}
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

ClassFactory<HydroRiemannSolver> &GetHydroRiemannSolverFactory();

#endif // MFP_FLUX_H
