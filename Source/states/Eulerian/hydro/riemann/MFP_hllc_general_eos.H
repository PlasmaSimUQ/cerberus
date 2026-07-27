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

#include "MFP_hlle_general_eos.H"
#include "MFP_hydro_riemann.H"

class HydroGas;

class HydroHLLCGeneralEOS : public HydroRiemannSolver
{
  public:
    HydroHLLCGeneralEOS();
    HydroHLLCGeneralEOS(const sol::table& def);

    virtual void solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk) override;

    // eval-reuse overload: consume face internal-energy/sound-speed pairs the
    // caller has already computed (the hybrid blend), avoiding a redundant
    // table inversion per side; the plain solve() above computes them then
    // delegates here.
    void solve(Vector<Real>& L,
               Vector<Real>& R,
               Vector<Real>& F,
               Real* shk,
               Real eL,
               Real aL,
               Real eR,
               Real aR);

    // forward the table pointer to the embedded HLLE used as the degeneracy
    // fallback (Section 2 of doc/general_eos_hllc_fallback_plan.md)
    void set_gas(const HydroGas* g)
    {
        gas = g;
        hlle.set_gas(g);
    }

    // diagnostic accumulator: faces that tripped a degeneracy guard and fell
    // back to HLLE since the last report, split by guard (read-and-reset)
    virtual RiemannFallbackCounts get_and_reset_fallback_counts() override
    {
        RiemannFallbackCounts n = n_fallback;
        n_fallback = RiemannFallbackCounts();
        // fold in the embedded fallback's own distress signal
        n.fan = hlle.get_and_reset_fan_widened();
        return n;
    }

    virtual const std::string get_tag() const override { return tag; }
    static std::string tag;

    virtual int get_n_flux() const override { return n_flux; }
    int n_flux;

    static bool registered;

    const HydroGas* gas = nullptr;

    // star-less fallback for degenerate faces; armed via set_gas
    HydroHLLEGeneralEOS hlle;

    // relative floor for the degeneracy guards (Lua key 'fallback_eps')
    Real eps_rel = 1.0e-6;

    // W29 per-guard switches (Lua keys 'fallback_guard_{eval,sound,wave}',
    // default true = established behaviour). A disabled guard runs the
    // star-state algebra unprotected on its mode — deliberate, see the ctor.
    bool guard_eval = true;
    bool guard_sound = true;
    bool guard_wave = true;

    // per-rank, per-guard counts of fallback substitutions (plain increments in
    // solve; never reduced inside the per-face path)
    RiemannFallbackCounts n_fallback;

    Vector<Real> svLs, fvL, svL;
    Vector<Real> svRs, fvR, svR;
    Vector<Real> trL, trR;
};

#endif  // HYDROHLLCGENERALEOS_H
