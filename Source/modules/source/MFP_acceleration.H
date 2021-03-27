#ifndef ACCELERATION_H
#define ACCELERATION_H

#include "MFP_source.H"

#include "sol.hpp"

class Acceleration : public SourceTerm
{
public:
    Acceleration();
    Acceleration(const sol::table& def);
    ~Acceleration();

    Vector<Real> get_acc(Real x, Real y, Real z, Real t) const;

    static Vector<dual> accelerate(const Vector<dual> &y0,
                            const Vector<OffsetIndex> &offsets,
                            const Vector<Real> &a);

    virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const override;
    virtual int fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const override;

    static bool valid_state(const int global_idx);
    static bool valid_solver(const int solve_idx);

    virtual SourceType get_type() const override {return SourceType::Acceleration;}
    virtual const std::string& get_tag() const override {return tag;}

    static std::string tag;
    static bool registered;

    Array<Optional3D1VFunction, AMREX_SPACEDIM> acc; // the acceleration functions
};

#endif // ACCELERATION_H
