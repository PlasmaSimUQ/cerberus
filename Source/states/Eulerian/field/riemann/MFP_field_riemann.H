#ifndef MFP_FLUX_H
#define MFP_FLUX_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "sol.hpp"
#include "MFP_factory.H"
#include "MFP_field_defs.H"
#include "MFP_eulerian_riemann.H"

using namespace amrex;

//================================================================================

class FieldRiemannSolver : public RiemannSolver
{
public:
    FieldRiemannSolver();
    virtual ~FieldRiemannSolver();

    int idx;

};

template <typename D>
std::unique_ptr<FieldRiemannSolver> FieldRiemannSolverBuilder(const sol::table& def)
{
    int global_idx = def["global_idx"];

    if ((def["flux"] == D::tag) && (D::valid_state(global_idx))) {
        return std::unique_ptr<D>(new D(global_idx));
    } else {
        return nullptr;
    }
}

ClassFactory<FieldRiemannSolver> &GetFieldRiemannSolverFactory();

#endif // MFP_FLUX_H
