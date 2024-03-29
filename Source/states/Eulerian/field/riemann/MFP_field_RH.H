#ifndef FIELDRH_H
#define FIELDRH_H

#include "MFP_field_riemann.H"

class FieldRH : public FieldRiemannSolver
{
  public:
    FieldRH();
    FieldRH(const int i);

    virtual void solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk) override;

    static bool valid_state(const int idx);
    virtual const std::string get_tag() const override { return tag; }
    virtual int get_n_flux() const override { return +FieldDef::ConsIdx::NUM; }
    static std::string tag;
    static bool registered;

    Real c0, ch, ch2;
};

#endif  // FIELDRH_H
