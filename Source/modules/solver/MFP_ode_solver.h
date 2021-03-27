#ifndef MFP_ODE_H
#define MFP_ODE_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <Core>

#include "sol.hpp"
#include "MFP_factory.H"


using namespace amrex;

class ODESystem;


// ====================================================================================

// NEED TO IMPLEMENT A STIFF SOLVER

enum class SolveType : int {
    Null=-1,
    Explicit,
    Implicit,
    CVODE,
    NUM
};

class SolveODE
{
public:
    SolveODE();
    ~SolveODE();
    virtual int init_data(Vector<Real> *y,
                          Vector<Real> *yout);
    virtual int solve(Real x, Real y, Real z, Real t0, Real t1, int depth);
    virtual bool valid_solution();
    virtual void clear();
    virtual SolveType get_type() const {return SolveType::Null;}
    virtual bool has_freq() const {return false;}
    virtual std::string get_tag() const = 0;

    ODESystem* parent;
    Vector<Real> *y0, *y1, *extra;
    int verbosity = 0;
};

template <typename D>
std::unique_ptr<SolveODE> SolveODEBuilder(const sol::table& def)
{
    if (def["solver"] == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

PhysicsFactory<SolveODE>& GetSolveODEFactory();

#endif // MFP_ODE_H

