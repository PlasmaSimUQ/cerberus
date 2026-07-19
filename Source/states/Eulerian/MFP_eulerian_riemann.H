#ifndef RIEMANNSOLVER_H
#define RIEMANNSOLVER_H

#include "MFP_factory.H"
#include "MFP_hydro_defs.H"
#include "MFP_utility.H"
#include "json.hpp"
#include "sol.hpp"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

using namespace amrex;

// Per-guard tally of degeneracy fallbacks (general-EOS HLLC -> HLLE). Split by
// guard so the *mechanism* is attributable: a run dominated by 'sound' is
// hitting the cold table hull, whereas 'wave' points at the wave-speed
// estimate. See Section 4 of doc/general_eos_hllc_fallback_plan.md.
//
// Compiled in only under MFP_SOLVER_DIAG (follows DEBUG; see the Makefile) —
// the increments sit in the per-face inner loop. The type itself is always
// defined so the accessor signatures do not change between builds.
struct RiemannFallbackCounts {
    long eval = 0;   // guard 1: non-finite or negative EOS evaluation
    long sound = 0;  // guard 2: relative sound-speed floor (cold hull)
    long wave = 0;   // guard 4: wave ordering / relative gap floor

    // NOTE: no 'pstar' counter — the p_star <= 0 guard was removed as incorrect
    // (it is not a degeneracy for HLLC). See MFP_hllc_general_eos.cpp.

    // NOT an HLLC guard: faces on which the HLLE *fallback* had to widen a
    // collapsed wave fan (sR - sL -> 0) to keep its own HLL average finite.
    // Counted separately because it measures the fallback failing, not the
    // primary solver.
    long fan = 0;

    long total() const { return eval + sound + wave; }

    RiemannFallbackCounts& operator+=(const RiemannFallbackCounts& o)
    {
        eval += o.eval;
        sound += o.sound;
        wave += o.wave;
        fan += o.fan;
        return *this;
    }
};

class RiemannSolver
{
  public:
    RiemannSolver();
    ~RiemannSolver();

    Real max_speed;

    virtual void solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk) {}

    virtual const bool need_max_speed() const { return false; }
    virtual void update_max_speed(const Real speed) {}
    virtual const std::string get_tag() const { return ""; }
    virtual int get_n_flux() const = 0;
    virtual bool requires_shock_detector() const { return false; }

    // diagnostic: faces on which the solver substituted a robust fallback flux
    // since the last call, split by which guard tripped; reads and resets the
    // accumulators. Mirrors EulerianState::get_and_reset_thermo_clamps. Zeroed
    // by default for solvers with no fallback path.
    //
    // Deliberately the ONLY reader: a second read-and-reset accessor (e.g. a
    // total-only one) would race with this one, whichever ran first zeroing the
    // counts the other reports.
    virtual RiemannFallbackCounts get_and_reset_fallback_counts()
    {
        return RiemannFallbackCounts();
    }
};

//================================================================================

class State;

class ShockDetector
{
  public:
    ShockDetector();
    ~ShockDetector();
    virtual Real solve(Vector<Real>& L, Vector<Real>& R) const { return 0.0; }

    int idx;
    State* istate;

    virtual const std::string get_tag() const { return ""; }
    virtual void write_info(nlohmann::json& js) const { return; }
};

#endif  // RIEMANNSOLVER_H
