#ifndef SOLVEEXPLICIT_H
#define SOLVEEXPLICIT_H

#include "MFP_ode_solver.h"

/*
 * Forward Euler timestep update
 */

class SolveRK : public SolveODE
{
public:
    SolveRK();
    SolveRK(const sol::table& def);
    ~SolveRK();
    virtual void euler_step_solve(Real x, Real y, Real z, Real t0, Real dt, Vector<Real>& y0, Vector<Real> &y1) = 0;
    virtual int solve(Real x, Real y, Real z, Real t0, Real t1, int depth=0) override;

    int n_sub = 2;
    int max_depth = 4;
    int verbosity = 0;
    int order = 1;
};

class SolveExplicit : public SolveRK
{
public:
    SolveExplicit();
    SolveExplicit(const sol::table& def);
    ~SolveExplicit();
    virtual void euler_step_solve(Real x, Real y, Real z, Real t0, Real dt, Vector<Real>& y0, Vector<Real> &y1) override;
    virtual SolveType get_type() const override {return SolveType::Explicit;}
    virtual bool has_freq() const override {return true;}
    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;
};


class SolveImplicit : public SolveRK
{
public:
    SolveImplicit();
    SolveImplicit(const sol::table& def);
    ~SolveImplicit();
    virtual void euler_step_solve(Real x, Real y, Real z, Real t0, Real dt, Vector<Real>& y0, Vector<Real> &y1) override;
    virtual SolveType get_type() const override {return SolveType::Implicit;}
    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;
    static Eigen::ColPivHouseholderQR<Eigen::MatrixXd> linear_solve;
};

#endif // SOLVEEXPLICIT_H
