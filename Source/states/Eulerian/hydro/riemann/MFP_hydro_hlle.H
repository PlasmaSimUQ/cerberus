#ifndef HYDROHLLE_H
#define HYDROHLLE_H

#include "MFP_hydro_defs.H"
#include "MFP_hydro_riemann.H"

class HydroHLLE : public HydroRiemannSolver
{
  public:
    HydroHLLE();
    HydroHLLE(const sol::table& def);

    virtual void solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk) override;

    virtual const std::string get_tag() const override { return tag; }
    static std::string tag;

    virtual int get_n_flux() const override { return n_flux; }
    int n_flux;

    static bool registered;

    Vector<Real> fvL, svL;
    Vector<Real> fvR, svR;
    Vector<Real> trL, trR;
};

#endif  // HYDROHLLE_H
