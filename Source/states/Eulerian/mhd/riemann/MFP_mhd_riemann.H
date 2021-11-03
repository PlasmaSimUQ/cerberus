#ifndef MFP_MHDRIEMANN_H
#define MFP_MHDRIEMANN_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "sol.hpp"
#include "MFP_factory.H"
#include "MFP_mhd_defs.H"
#include "MFP_utility.H"

using namespace amrex;

//================================================================================

class MHDRiemannSolver
{
public:
    MHDRiemannSolver();
    virtual ~MHDRiemannSolver();

    virtual void solve(Array<Real,+MHDDef::PrimIdx::NUM> &L,
                       Array<Real,+MHDDef::PrimIdx::NUM> &R,
                       Array<Real,+MHDDef::ConsIdx::NUM> &F,
                       Real* shk) {}

    virtual std::string get_tag() const {return "";}
    virtual bool requires_shock_detector() const {return false;}

    void update_max_speed(const Real s) {max_speed = s;}

    Real max_speed;
};

template <typename D>
std::unique_ptr<MHDRiemannSolver> MHDRiemannSolverBuilder(const sol::table& def)
{
    if (def["flux"] == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

ClassFactory<MHDRiemannSolver> &GetMHDRiemannSolverFactory();

#endif // MFP_MHDRIEMANN_H
