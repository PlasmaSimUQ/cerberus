#ifndef REACTIONS_H
#define REACTIONS_H

#include "MFP_source.H"

#include "MFP_interp.H"


class Ionisation : public SourceTerm
{
public:
    Ionisation();
    Ionisation(const sol::table &def);
    ~Ionisation();

    static Vector<Real> ionisation(const Vector<Real> &y0, const Vector<OffsetIndex> &offsets, const Interp2D &dens_interp, const Interp2D &mom_interp, const Interp2D &nrg_interp, const Real &bind_nrg);

    virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const override;
    virtual int fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const override;

    static bool valid_state(const int global_idx);
    static bool valid_solver(const int solve_idx);

    virtual SourceType get_type() const override {return SourceType::Ionisation;}
    virtual const std::string& get_tag() const override {return tag;}

    static std::string tag;
    static bool registered;

    Interp2D dens_interp;
    Interp2D mom_interp;
    Interp2D nrg_interp;
    Real bind_nrg;
};

class Recombination : public SourceTerm
{
public:
    Recombination();
    Recombination(const sol::table &def);
    ~Recombination();

    static Vector<Real> recombination(const Vector<Real> &y0, const Vector<OffsetIndex> &offsets, const Interp2D &dens_interp, const Interp2D &R_interp, const Interp2D &K_interp, const Interp2D &W_interp, const Interp2D &J_interp, const Real &bind_nrg);

    virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const override;
    virtual int fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const override;

    static bool valid_state(const int global_idx);
    static bool valid_solver(const int solve_idx);

    virtual SourceType get_type() const override {return SourceType::Recombination;}
    virtual const std::string& get_tag() const override {return tag;}

    static std::string tag;
    static bool registered;

    Interp2D dens_interp;
    Interp2D R_interp;
    Interp2D K_interp;
    Interp2D W_interp;
    Interp2D J_interp;
    Real bind_nrg;
};


#endif
