#ifndef HYDROHLLCGENERALEOS_H
#define HYDROHLLCGENERALEOS_H

// General-EOS HLLC Riemann solver (plan W10, Stage 5 — STAGE5.md; Athena++
// hllc.cpp GENERAL_EOS pattern, plan D2).
//
// Structurally a copy of MFP_hllc.cpp with exactly three ingredient swaps —
// the star-state algebra (S_star, the starred states, the flux assembly,
// tracers) is EOS-independent Rankine–Hugoniot bookkeeping and is identical:
//
//   1. face total energy:  p/(gam-1) + ke   ->  rho * gas->e(Q) + ke
//   2. face sound speed:   sqrt(gam*p/rho)  ->  gas->a(Q)
//      (1+2 come from ONE combined gas evaluation per face — for a
//      tabulated gas that is a single table inversion, the G8 budget)
//   3. two-shock q-factor: gamma slot       ->  local acoustic gamma
//      Gamma_1 = a^2*rho/p at the face, reusing the sound speed from 2
//      (no extra table calls; STAGE5.md D-c recorded alternative, adopted
//      after the G8 measurement — the PVRS-p* variant was physics-identical
//      within every gate but cost ~1 extra inversion per face)
//
// The gas model is a non-owning pointer injected by HydroState::set_flux()
// (construction order is safe: the gas is built before the flux solver).
// On a gamma-law gas the HydroGas base defaults reproduce MFP_hllc.cpp's
// algebra exactly, so this solver agrees with HLLC to round-off — the gate
// that isolates solver bugs from EOS bugs (STAGE5.md G1).
//
// FPE safety (traps are enabled in the EOS test cases): faces arriving here
// have been hull-floored by apply_prim_floor (rho, p > 0); table sound
// speeds are floored non-negative inside eval_rt; the q-factor divisions are
// guarded by the p_star > p_face branch condition and g_loc > 0.

#include "MFP_hydro_riemann.H"

class HydroGas;

class HydroHLLCGeneralEOS : public HydroRiemannSolver
{
  public:
    HydroHLLCGeneralEOS();
    HydroHLLCGeneralEOS(const sol::table& def);

    virtual void solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk) override;

    void set_gas(const HydroGas* g) { gas = g; }

    virtual const std::string get_tag() const override { return tag; }
    static std::string tag;

    virtual int get_n_flux() const override { return n_flux; }
    int n_flux;

    static bool registered;

    const HydroGas* gas = nullptr;

    Vector<Real> svLs, fvL, svL;
    Vector<Real> svRs, fvR, svR;
    Vector<Real> trL, trR;
};

#endif  // HYDROHLLCGENERALEOS_H
