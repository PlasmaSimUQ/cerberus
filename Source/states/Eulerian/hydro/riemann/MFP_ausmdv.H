#ifndef HYDROAUSMDV_H
#define HYDROAUSMDV_H

#include "MFP_hydro_riemann.H"

class HydroAUSMDV : public HydroRiemannSolver
{
public:
    HydroAUSMDV();
    HydroAUSMDV(const sol::table& def);

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) override;

    virtual std::string get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

    Vector<Real> trL;
    Vector<Real> trR;
};

#endif // HYDROAUSMDV_H
