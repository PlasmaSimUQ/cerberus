#ifndef HYDROHLLC_H
#define HYDROHLLC_H

#include "MFP_hydro_riemann.H"

class HydroHLLC : public HydroRiemannSolver
{
public:
    HydroHLLC();
    HydroHLLC(const sol::table& def);

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) override;

    static bool valid_state(const int idx);
    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

    Vector<Real> svLs, fvL, svL;
    Vector<Real> svRs, fvR, svR;
    Vector<Real> trL, trR;

};

#endif // HYDROHLLC_H
