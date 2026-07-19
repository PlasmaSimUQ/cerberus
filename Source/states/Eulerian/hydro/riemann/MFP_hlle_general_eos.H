#ifndef HYDROHLLEGENERALEOS_H
#define HYDROHLLEGENERALEOS_H

// General-EOS HLLE Riemann solver (companion to MFP_hllc_general_eos.cpp).
//
// Structurally a copy of MFP_hydro_hlle.cpp with the two closure swaps that
// apply to a solver with NO star state (HLLE has no two-shock q-factor, so the
// third swap used by the general-EOS HLLC does not appear here):
//
//   1. face total energy:  p/(gam-1) + ke   ->  rho * gas->e(Q) + ke
//   2. face sound speed:   sqrt(gam*p/rho)  ->  gas->a(Q)
//      (1+2 come from ONE combined gas evaluation per face — for a tabulated
//      gas that is a single table inversion)
//
// The HLL average and the outermost signal-speed estimates
//   sL = min(uL-aL, uR-aR),  sR = max(uL+aL, uR+aR)
// are closure-independent and identical to the ideal HLLE. On a gamma-law gas
// the HydroGas base defaults reproduce sqrt(gam*p/rho) and p/((gam-1)rho), so
// this solver agrees with HydroHLLE to round-off — the gate that isolates
// solver bugs from EOS-table bugs.
//
// FPE safety: HLLE performs no division by an internal wave speed (unlike
// HLLC's starred energy p/(rho*(S_L-u_L))), so there is no a->0 singularity to
// guard; faces arriving here are already hull-floored (rho, p > 0) and table
// sound speeds are floored non-negative inside the evaluator.

#include "MFP_hydro_defs.H"
#include "MFP_hydro_riemann.H"

class HydroGas;

class HydroHLLEGeneralEOS : public HydroRiemannSolver
{
  public:
    HydroHLLEGeneralEOS();
    HydroHLLEGeneralEOS(const sol::table& def);

    virtual void solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk) override;

    // eval-reuse overload: consume face internal-energy/sound-speed pairs the
    // caller has already computed (HLLC fallback, hybrid blend) instead of
    // re-inverting the table for L and R — one table inversion is the dominant
    // cost of the tabulated path, so this avoids doubling it.
    void solve(Vector<Real>& L,
               Vector<Real>& R,
               Vector<Real>& F,
               Real* shk,
               Real eL,
               Real aL,
               Real eR,
               Real aR);

    void set_gas(const HydroGas* g) { gas = g; }

    // faces on which the wave fan had to be widened to keep the HLL average
    // finite; read-and-reset, merged into the owner's counts.
    //
    // SPECULATIVE: this has read zero in every run to date. The guard was added
    // on a hypothesis that measurement later refuted — see the long note at the
    // widening site in MFP_hlle_general_eos.cpp before relying on it.
    long get_and_reset_fan_widened()
    {
        const long n = n_fan_widened;
        n_fan_widened = 0;
        return n;
    }

    long n_fan_widened = 0;

    virtual const std::string get_tag() const override { return tag; }
    static std::string tag;

    virtual int get_n_flux() const override { return n_flux; }
    int n_flux;

    static bool registered;

    const HydroGas* gas = nullptr;

    Vector<Real> fvL, svL;
    Vector<Real> fvR, svR;
    Vector<Real> trL, trR;
};

#endif  // HYDROHLLEGENERALEOS_H
