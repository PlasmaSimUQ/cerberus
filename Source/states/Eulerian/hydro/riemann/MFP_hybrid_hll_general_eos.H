#ifndef HYDROHYBRIDHLLGENERALEOS_H
#define HYDROHYBRIDHLLGENERALEOS_H

// General-EOS shock-switching hybrid (companion to MFP_hybrid_hll.cpp).
//
// Runs the general-EOS HLLC in smooth flow (sharp contacts/shears) and falls
// back to the general-EOS HLLE across strong shocks (robust, no star state),
// blending linearly by the shock weight *shk produced by a shock detector.
// This is the general-EOS answer to the HLLC starred-energy blow-up
// (p/(rho*(S_L-u_L)) as the tabulated sound speed -> 0 in a cold hyper-Mach
// rarefaction): the detector routes those faces through the star-less HLLE.
//
// Differences from the ideal HydroHybridHLL:
//   * set_gas() forwards the non-owning gas pointer to BOTH sub-solvers.
//   * the transition branch blends the FULL flux vector (cons + tracers),
//     n_flux long, so passive scalars are advected consistently with the hydro
//     state — the ideal hybrid blends only ConsIdx::NUM and leaves the tracer
//     slots of F stale (see doc/general_eos_hybrid_riemann_plan.md).

#include "MFP_hllc_general_eos.H"
#include "MFP_hlle_general_eos.H"
#include "MFP_hydro_defs.H"
#include "MFP_hydro_riemann.H"

class HydroGas;

class HydroHybridHLLGeneralEOS : public HydroRiemannSolver
{
  public:
    HydroHybridHLLGeneralEOS();
    HydroHybridHLLGeneralEOS(const sol::table& def);

    HydroHLLCGeneralEOS hllc;
    HydroHLLEGeneralEOS hlle;

    virtual void solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk) override;

    virtual bool requires_shock_detector() const { return true; }

    void set_gas(const HydroGas* g)
    {
        gas = g;
        hllc.set_gas(g);
        hlle.set_gas(g);
    }

    // fallbacks accrue inside the embedded HLLC (the only sub-solver with a
    // fallback path); surface its per-guard counts for the per-block diagnostic
    virtual RiemannFallbackCounts get_and_reset_fallback_counts() override
    {
        return hllc.get_and_reset_fallback_counts();
    }

    // used to evaluate the table ONCE per side in the blend branch, then reused
    // by both sub-solvers
    const HydroGas* gas = nullptr;

    virtual const std::string get_tag() const override { return tag; }
    static std::string tag;

    virtual int get_n_flux() const override { return n_flux; }
    int n_flux;

    static bool registered;

    Vector<Real> F_hlle, F_hllc;
};

#endif  // HYDROHYBRIDHLLGENERALEOS_H
