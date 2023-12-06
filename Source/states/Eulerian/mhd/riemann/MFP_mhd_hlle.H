#ifndef MHDHLLE_H
#define MHDHLLE_H

#include "MFP_mhd_defs.H"
#include "MFP_mhd_riemann.H"

class MHDHLLE : public MHDRiemannSolver
{
  public:
    MHDHLLE();
    MHDHLLE(const sol::table& def);

    virtual void solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk) override;

    virtual const std::string get_tag() const override { return tag; }
    static std::string tag;

    virtual int get_n_flux() const override { return +MHDDef::ConsIdx::NUM; }
    virtual const bool need_max_speed() const override { return div_transport > 0.0; }

    static bool registered;

    Real gamma;
    Real div_transport;
};

#endif  // MHDHLLE_H
