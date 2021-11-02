#ifndef FIELDRH_H
#define FIELDRH_H

#include "MFP_field_riemann.H"

class FieldRH : public FieldRiemannSolver
{
public:

    FieldRH();
    FieldRH( const int i);

    virtual void solve(Array<Real,+FieldDef::ConsIdx::NUM> &L,
                       Array<Real,+FieldDef::ConsIdx::NUM> &R,
                       Array<Real,+FieldDef::ConsIdx::NUM> &F) const override;

    static bool valid_state(const int idx);
    virtual const std::string& get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

    Real c0, ch, ch2;
};

#endif // FIELDRH_H
