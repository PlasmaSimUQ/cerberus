#ifndef MHDHLLE_H
#define MHDHLLE_H

#include "MFP_mhd_riemann.H"
#include "MFP_mhd_defs.H"

class MHDHLLE : public MHDRiemannSolver
{
public:
    MHDHLLE();
    MHDHLLE(const sol::table& def);

    virtual void solve(Array<Real,+MHDDef::PrimIdx::NUM> &L,
                       Array<Real,+MHDDef::PrimIdx::NUM> &R,
                       Array<Real,+MHDDef::ConsIdx::NUM> &F,
                       Real* shk) override;

    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

    Real gamma;
    Real div_transport;

};

#endif // MHDHLLE_H
