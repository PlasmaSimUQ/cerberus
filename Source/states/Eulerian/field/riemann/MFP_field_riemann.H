#ifndef MFP_FLUX_H
#define MFP_FLUX_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "sol.hpp"
#include "MFP_factory.H"
#include "MFP_field_defs.H"

using namespace amrex;

//================================================================================

class FieldRiemannSolver
{
public:
    FieldRiemannSolver();
    virtual ~FieldRiemannSolver();

    int idx;

    virtual void solve(Array<Real,+FieldDef::ConsIdx::NUM> &L,
                       Array<Real,+FieldDef::ConsIdx::NUM> &R,
                       Array<Real,+FieldDef::ConsIdx::NUM> &F) const = 0;

    virtual const std::string& get_tag() const = 0;
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
