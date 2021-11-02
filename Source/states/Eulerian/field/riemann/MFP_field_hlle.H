#ifndef FIELDHLLE_H
#define FIELDHLLE_H

#include "MFP_field_riemann.H"

class FieldHLLE : public FieldRiemannSolver
{
public:

    FieldHLLE();
    FieldHLLE(const int i);

    virtual void solve(Array<Real,+FieldDef::ConsIdx::NUM> &L,
                       Array<Real,+FieldDef::ConsIdx::NUM> &R,
                       Array<Real,+FieldDef::ConsIdx::NUM> &F) const override;

    static bool valid_state(const int idx);
    virtual const std::string& get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

    Real c0, c2, ch, ch2, cc;
};

#endif // FIELDHLLE_H
