#ifndef HYDROHLLE_H
#define HYDROHLLE_H

#include "MFP_hydro_riemann.H"
#include "MFP_hydro_defs.H"

class HydroHLLE : public HydroRiemannSolver
{
public:
    HydroHLLE();
    HydroHLLE(const sol::table& def);

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) override;

    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

    Vector<Real> fvL, svL;
    Vector<Real> fvR, svR;
    Vector<Real> trL, trR;

};

#endif // HYDROHLLE_H
