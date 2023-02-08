#ifndef MFP_MHDRIEMANN_H
#define MFP_MHDRIEMANN_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "sol.hpp"
#include "MFP_factory.H"
#include "MFP_mhd_defs.H"
#include "MFP_utility.H"
#include "MFP_eulerian_riemann.H"

using namespace amrex;

//================================================================================

class MHDRiemannSolver : public RiemannSolver
{
public:
    MHDRiemannSolver();
    virtual ~MHDRiemannSolver();
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
