#ifndef FIELDHLLE_H
#define FIELDHLLE_H

#include "MFP_field_riemann.H"

class FieldHLLE : public FieldRiemannSolver
{
  public:
    FieldHLLE();
    FieldHLLE(const int i);

    virtual void solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk) override;

    static bool valid_state(const int idx);
    virtual const std::string get_tag() const override { return tag; }
    virtual int get_n_flux() const override { return +FieldDef::ConsIdx::NUM; }
    static std::string tag;
    static bool registered;

    Real c0, c2, ch, ch2, cc;
};

#endif  // FIELDHLLE_H
