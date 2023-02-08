#ifndef RIEMANNSOLVER_H
#define RIEMANNSOLVER_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "sol.hpp"
#include "json.hpp"
#include "MFP_factory.H"
#include "MFP_hydro_defs.H"
#include "MFP_utility.H"

using namespace amrex;

class RiemannSolver
{
public:
    RiemannSolver();
    ~RiemannSolver();

    Real max_speed;

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) {}

    virtual const bool need_max_speed() const {return false;}
    virtual void update_max_speed(const Real speed) {}
    virtual const std::string get_tag() const {return "";}
    virtual int get_n_flux() const = 0;
    virtual bool requires_shock_detector() const {return false;}
};



//================================================================================

class State;


class ShockDetector
{
public:
    ShockDetector();
    ~ShockDetector();
    virtual Real solve(Vector<Real> &L,Vector<Real> &R) const {return 0.0;}

    int idx;
    State* istate;

    virtual const std::string get_tag() const {return "";}
    virtual void write_info(nlohmann::json& js) const {return;}
};


#endif // RIEMANNSOLVER_H
