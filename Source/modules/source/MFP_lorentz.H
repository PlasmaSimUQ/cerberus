#ifndef LORENTZ_H
#define LORENTZ_H

#include "MFP_plasma5.H"
#include "MFP_source.H"

//-----------------------------------------------------------------------------

class Lorentz : public Plasma5
{
  public:
    Lorentz();
    Lorentz(const sol::table& def);
    ~Lorentz();
    static Vector<Real> source(const Vector<Real>& y, const Vector<OffsetIndex>& apply);
    virtual int fun_rhs(Real x,
                        Real y,
                        Real z,
                        Real t,
                        Vector<Real>& y0,
                        Vector<Real>& ydot,
                        Real dt = 0) const override;
    virtual int
      fun_jac(Real x, Real y, Real z, Real t, Vector<Real>& y0, Vector<Real>& J) const override;

    virtual bool has_freq() const override { return true; }
    virtual Real get_max_freq(Vector<Real>& y) const override;

    static bool valid_state(const int global_idx) { return Plasma5::valid_state(global_idx); }
    static bool valid_solver(const int solve_idx) { return Plasma5::valid_state(solve_idx); }

    virtual SourceType get_type() const override { return SourceType::Lorentz; }
    virtual const std::string& get_tag() const override { return tag; }

    static std::string tag;
    static bool registered;

    Vector<int> linked_electron_idx, linked_ion_idx;
    int linked_EM_idx;
};

#endif  // PLASMA5_H
