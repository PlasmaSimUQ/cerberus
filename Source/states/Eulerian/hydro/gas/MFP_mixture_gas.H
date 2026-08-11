#ifndef MFP_MIXTURE_GAS_H
#define MFP_MIXTURE_GAS_H

// Composition-weighted mixture EOS gas model (plan W24-W26, Stage 9 —
// doc/eos_mixture_dalton_plan.md; design decisions in
// doc/eos_wide_range_and_mixtures.md section 5).
//
// A HydroGas backend that binds ONE conditioned .eostab table PER COMPONENT
// and closes the state under a mixing rule, with the advected tracer mass
// fractions alpha_k as the composition. Dalton (partial pressure, shared
// volume + temperature) is the v1 bring-up rule; Amagat (additive volume,
// plan W27) is the follow-on accurate mode for condensed mixtures and adds a
// branch behind `rule` with the same driver signatures. Configured per state:
//
//   gas = {
//       type        = 'tabulated_mixture',
//       mixing_rule = 'dalton',            -- 'amagat': switch plumbed (AM1),
//                                          --   drivers land with AM2 —
//                                          --   doc/eos_amagat_plan.md
//       components  = {
//           { table = 'heavy.eostab', name = 'heavy', mass = 1.0, charge = 0.0 },
//           { table = 'light.eostab', name = 'light', mass = 0.5, charge = 0.0 },
//       },
//       ttol = 1e-10, max_newton = 100,    -- inversion controls (optional)
//       pure_tol = 1e-10, drop_tol = 1e-10,
//   }
//   value = { ..., alpha = alpha_heavy }   -- N-1 tracers; last = 1 - sum
//
// Component ORDER defines the alpha mapping: components 1..N-1 (Lua order)
// bind to tracer slots 0..N-2 and the last component carries the derived
// fraction 1 - sum(alpha).
//
// Dalton closure: rho_k = alpha_k rho (shared volume), one shared T, and
//   p = sum_k p_k(rho_k, T),  e = sum_k alpha_k e_k(rho_k, T).
// Both inverse drivers (given e or p, recover T) are 1-D guarded
// Newton/bisection solves in the shared T: the slopes sum_k alpha_k cv_k and
// sum_k dpdT_k are positive by table conditioning, so the residual-converged
// solve is unconditionally convergent on the bracket. The bracket is the
// intersection of the retained components' T-hulls; the ctor aborts if the
// full-component intersection is empty (Dalton needs one shared T inside
// every hull).
//
// KNOWN LIMITATION (stated up front, plan section 3): Dalton puts component
// k at alpha_k*rho in the full volume, so a dilute CONDENSED component lands
// far below its table hull — constant (legitimate, per-component-tallied)
// clamps and an unphysical partial pressure. That is intrinsic to the
// shared-volume rule, which is exact only for ideal-gas components. Amagat
// (W27) is the production closure for condensed heavy+light pairings; the
// per-component clamp counters are the signal to prioritise it.
//
// The frozen mixture sound speed (composition fixed across the wave) is
//   a^2 = sum_k alpha_k dpdrho_k + T (sum_k dpdT_k)^2 / (rho^2 sum_k alpha_k cv_k)
// — all existing EosEval columns, no new table data (plan W26). The state
// runs under the Stage-5 HLLC_general_eos solver via get_face_eval_from_prim
// (needs_general_eos_solver() below); the solver needs zero edits.
//
// NOT compatible with MHD states or the Braginskii action (config-time
// guards, plan section 8).

#include "MFP_eos_table.H"
#include "MFP_hydro_gas.H"

class MixtureEOS : public HydroGas
{
  public:
    MixtureEOS();
    MixtureEOS(const int global_idx, const sol::table& def);
    ~MixtureEOS();

    virtual Real get_gamma_from_cons(const Vector<Real>& U,
                                     const int density_idx = +HydroDef::ConsIdx::Density,
                                     const int tracer_idx = +HydroDef::ConsIdx::NUM) const override;
    virtual Real get_cp_from_cons(const Vector<Real>& U,
                                  const int density_idx = +HydroDef::ConsIdx::Density,
                                  const int tracer_idx = +HydroDef::ConsIdx::NUM) const override;

    virtual Real get_gamma_from_prim(const Vector<Real>& Q,
                                     const int idx = +HydroDef::PrimIdx::NUM) const override;
    virtual Real get_cp_from_prim(const Vector<Real>& Q,
                                  const int tracer_idx = +HydroDef::PrimIdx::NUM) const override;

    virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const override;
    virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const override;
    virtual void define_rho_p_T(Vector<Real>& Q) const override;

    virtual Real get_temperature_from_cons(const Vector<Real>& U) const override;
    virtual RealArray get_speed_from_cons(const Vector<Real>& U) const override;
    virtual RealArray get_speed_from_prim(const Vector<Real>& Q) const override;

    // face evaluations for HLLC_general_eos — one rp mixture solve answers
    // both quantities (the Stage-5 G8 wall-time discipline)
    virtual Real get_internal_energy_from_prim(const Vector<Real>& Q) const override;
    virtual Real get_sound_speed_from_prim_rp(const Vector<Real>& Q) const override;
    virtual void get_face_eval_from_prim(const Vector<Real>& Q,
                                         Real& e,
                                         Real& a) const override;

    virtual std::string get_tag() const override { return tag; }
    virtual void write_info(nlohmann::json& js) const override;

    virtual bool needs_general_eos_solver() const override { return true; }

    // conservative scalar hull floors (ctor-computed, plan section 5.6); the
    // real guard is the per-component clamp inside the drivers
    virtual int apply_prim_floor(Vector<Real>& Q) const override;

    virtual long get_and_reset_hull_clamps() const override
    {
        long n = n_hull_clamps;
        for (long& nk : n_hull_clamps_k) {
            n += nk;
            nk = 0;
        }
        n_hull_clamps = 0;
        return n;
    }

    virtual long get_and_reset_alpha_fixes() const override
    {
        const long n = n_alpha_fixes;
        n_alpha_fixes = 0;
        return n;
    }

    // TEMPORARY DIAGNOSTIC (W31 investigation) -- see the counter block below
    virtual std::string get_and_reset_pin_report() const override
    {
        if (n_asm == 0) return {};
        std::ostringstream o;
        o << "asm=" << n_asm << " partial=" << n_asm_partial
          << " allpin=" << n_asm_allpin << " partial_ovf(rho*cs2>4)=" << n_asm_ovf
          << " p<=0=" << n_asm_nonpos_p << " max_rho_cs2 partial=" << max_rcs2_partial
          << " nopin=" << max_rcs2_nopin
          << " edge_flag1=" << n_edge_flag1 << " edge_flag2=" << n_edge_flag2
          << " cap_cs2=" << n_cap_cs2 << " cap_K=" << n_cap_K
          << " cap_dpdrho=" << n_cap_dpdrho;
        if (n_asm_partial > 0) {
            o << " | partial pin/live by component:";
            for (size_t k = 0; k < n_pin_partial_k.size(); ++k) {
                o << " " << k << ":" << n_pin_partial_k[k] << "/" << n_live_partial_k[k];
            }
            o << " | worst partial: rho=" << worst_rho << " p=" << worst_p
              << " T=" << worst_T << " w={";
            for (size_t k = 0; k < worst_w.size(); ++k) {
                o << (k ? "," : "") << worst_w[k];
            }
            o << "}";
        }
        n_asm = n_asm_partial = n_asm_allpin = n_asm_ovf = n_asm_nonpos_p = 0;
        n_edge_flag1 = n_edge_flag2 = n_cap_cs2 = n_cap_K = n_cap_dpdrho = 0;
        max_rcs2_partial = max_rcs2_nopin = 0.0;
        std::fill(n_pin_partial_k.begin(), n_pin_partial_k.end(), 0L);
        std::fill(n_live_partial_k.begin(), n_live_partial_k.end(), 0L);
        return o.str();
    }

    // B (T4b plan Part I): floor-time momentum repairs performed by
    // cons2prim (resurrections zeroed + soft hits capped)
    virtual long get_and_reset_mom_fixes() const override
    {
        const long n = n_mom_fixes;
        n_mom_fixes = 0;
        return n;
    }

    // true iff the LAST cons2prim call repaired the conserved vector in
    // place — calc_primitives keys the conserved-fab write-back on this
    // (the 2026-07-24 forensics found rho <= 0 persisted in checkpoints
    // when only the primitive copy was repaired)
    virtual bool cons_was_repaired() const override { return last_repair; }

    static bool registered;
    static std::string tag;

  private:
    enum class MixRule { Dalton, Amagat };

    // --- the mixing-rule drivers (the seam Amagat/W27 re-implements) ------
    // (rho, e_int, alpha) -> full mixture evaluation via the shared-T re
    // solve; alpha must be sanitized (size N, simplex). st.flag follows the
    // EosInvertStats semantics with max-aggregation across components.
    void eval_from_rho_e(Real rho,
                         Real e_int,
                         const Vector<Real>& alpha,
                         EosEval& ev,
                         EosInvertStats& st) const;
    // (rho, p, alpha) -> full mixture evaluation via the shared-T rp solve
    void eval_from_rho_p(Real rho,
                         Real p,
                         const Vector<Real>& alpha,
                         EosEval& ev,
                         EosInvertStats& st) const;
    // shared body of the two above: solve for T with target e (by_e) or p
    void dalton_solve_T(Real rho,
                        Real target,
                        bool by_e,
                        const Vector<Real>& alpha,
                        EosEval& ev,
                        EosInvertStats& st) const;

    // pure-cell short-circuit body shared by BOTH mixing rules (plan 7.1):
    // formula-identical to the single-table TabulatedEOS driver — the
    // mechanism behind every pure-cell round-off/bitwise gate
    void pure_cell_eval(int kp,
                        Real rho,
                        Real target,
                        bool by_e,
                        EosEval& ev,
                        EosInvertStats& st) const;

    // --- Amagat (W27; doc/eos_amagat_plan.md). AM2 landed the rp path;
    // AM3 the nested re-driver + the (rho,T) -> p forward. ---------------
    // Volume sum Sum_k w_k v_k(p, T) at fixed (p, T): each retained rho_k
    // from the guarded isotherm inversion (warm-started through m_rho_k;
    // an unattainable p pins rho_k at the hull edge — the pin mask goes to
    // m_pinned, the raw table flag to m_pin_flag, and pinned components
    // contribute zero derivative weight IN THIS SOLVE — amagat_p_solve's
    // unconditional convergence relies on that. This is a DIFFERENT
    // convention from amagat_assemble, which (pinned_compliance_fix_plan.md
    // Step 1) gives pinned components their edge compliance instead of
    // dropping them; the pin flag therefore carries two meanings depending
    // on which function reads it. Fills m_evk with the (p,T)-consistent
    // per-component evals and
    // returns the two slope sums (dv/dT)_p and (dv/dp)_T. `tally` controls
    // the accounting: iteration-interior evaluations run quiet (transient
    // pins while bracketing are not state pins — the same rule as the
    // Dalton rho-inversion); only the final at-solution evaluation tallies
    // per component and escalates the caller's flag.
    Real amagat_volume_sum(Real p,
                           Real T,
                           bool tally,
                           EosInvertStats& st,
                           Real& dvdT,
                           Real& dvdp) const;
    // inner solve of the volume constraint at fixed T (plan eq 4): p from
    // R_v(p) = Sum_k w_k v_k(p,T) - 1/rho, STRICTLY monotone decreasing in
    // p (dpdrho > 0 by conditioning) — unconditionally convergent on the
    // union bracket of the attainable per-component p-ranges (plan eq 7).
    // p_warm carries the previous outer iterate's solution as the seed
    // (in/out; < 0 = bracket midpoint). A no-root state (tension /
    // over-compression) clamps to the bracket end and flags.
    Real amagat_p_solve(Real rho, Real T, Real& p_warm, EosInvertStats& st) const;
    // driver twin of dalton_solve_T behind `rule`: same pure-cell head;
    // by_e = false -> guarded 1-D T-solve of the volume residual at given
    // p (AM2); by_e = true -> nested solve, outer 1-D in T on the energy
    // sum with the inner p-solve per iterate (AM3)
    void amagat_solve_T(Real rho,
                        Real target,
                        bool by_e,
                        const Vector<Real>& alpha,
                        EosEval& ev,
                        EosInvertStats& st) const;
    // mixture EosEval from the cached per-component evals at the solution:
    // alpha-weighted volume-derivative sums (plan section 3); reduces to
    // gamma*p/rho for a single ideal-gas component. Pinned components carry
    // their EDGE compliance (Remedy A, pinned_compliance_fix_plan.md Step
    // 1) rather than being dropped — the table's (dp/drho)_k etc. at the
    // clamped rho_k, using the same m_rho_k/m_evk the volume sum cached,
    // so the sums stay bounded and mutually consistent. The assembled
    // modulus is then capped from above by the stiffest retained
    // component's own modulus (Remedy B, Step 2; upper bound only — see
    // mixture_stiffness_bounds_and_pinning.md section 3 for why no lower
    // bound is sound here).
    void amagat_assemble(Real rho, Real p, Real T, EosEval& ev, EosInvertStats& st) const;
    // given (T, p, alpha) recover the total rho — define_rho_p_T's p+T
    // branch (initial-condition fills; not a hot path)
    Real invert_rho_from_p_T(Real T, Real p, const Vector<Real>& alpha, EosInvertStats& st) const;

    // clamp raw tracer output into the composition simplex: each alpha_k
    // into [0,1] then renormalise to sum 1. REQUIRED at every gas-model
    // entry because the base get_alpha_fractions_* helpers do no clamping
    // (verified) — reconstruction/advection overshoot passes straight
    // through and the derived last fraction goes negative when the stored
    // tracers sum past 1. Returns the number of repairs.
    int sanitize_alpha(Vector<Real>& alpha) const;

    // build the retained-component weights from sanitized alpha:
    // components with alpha_k <= drop_tol are dropped (weight 0) and the
    // rest renormalised (C0-smooth to O(drop_tol)); prevents rho_k -> 0
    // driving a component off-hull for a below-round-off contribution.
    // Returns the component index if one weight >= 1 - pure_tol (the
    // pure-cell short-circuit), else -1. Fills m_w / m_retained.
    int prepare_weights(const Vector<Real>& alpha) const;

    // clamp a partial density into component k's hull, tallying per
    // component and escalating the stats flag (clamping freezes rho_k, not
    // the T-dependence, so the shared-T solve stays monotone)
    Real clamp_rho_k(int k, Real rho_k, EosInvertStats& st) const;

    // per-call setup of the clamped partial densities m_rho_k AND the
    // dilute-limit pressure scales m_scale_k: below its density hull a
    // component is evaluated at the hull edge with p and dpdT scaled by
    // rho_k/rho_hull_min — the ideal-gas low-density limit (p linear in rho
    // at fixed T; e, cv rho-independent). Without this a dilute component
    // would contribute the FULL p(rho_hull_min, T) — a finite, alpha-
    // independent spurious pressure — and the drop_tol crossing would kink
    // by ~p(rho_hull_min,T)/p instead of O(drop_tol) (the M2 gate).
    // High-side clamps stay plain: compression beyond the table is
    // genuinely off-model.
    void setup_partial_densities(Real rho, EosInvertStats& st) const;

    // forward Dalton pressure at given (rho, T, weights): sum of the scaled
    // component pressures (define_rho_p_T's rho+T branch and the self-test)
    Real mix_p_rho_T(Real rho, Real T, EosInvertStats& st) const;

    // clamp the TOTAL density into the union of the component hulls before
    // a solve (the Stage-4 clamp-rho-first rule; the per-component clamps
    // then own the real work)
    Real clamp_rho_union(Real rho) const;

    void tally(const EosInvertStats& st) const
    {
        if (st.flag != 0) ++n_hull_clamps;
    }

    // one-zone self-test sweep (plan M2), run from the ctor when the gas
    // def carries `self_test = <n_sweep>`: driver round-trips on an
    // alpha-grid x hull-sweep, the frozen sound speed vs an isentropic
    // finite difference of the mixture pressure, pure-cell formula parity
    // with the single-table path, and a kink-free drop_tol crossing. Prints
    // one "MIXEOS-SELFTEST[...] PASS/FAIL ..." line per check (check.py
    // greps the run log).
    void run_self_test(int n_sweep) const;

    std::vector<EosTable> tables;    // owners: loaded + nondimensionalised in ctor
    std::vector<EosTableView> tvs;   // index-aligned non-owning views

    // per-component hull caches in code units (ctor)
    std::vector<Real> rho_hull_min_k, rho_hull_max_k, p_hull_min_k;
    std::vector<Real> T_min_k, T_max_k;

    // union scalars derived from the caches (ctor)
    Real rho_union_min = 0.0, rho_union_max = 0.0;  // min_k / max_k of the rho hulls
    Real rho_floor = 0.0, p_floor = 0.0;            // union (min_k) floors; rho_floor == rho_union_min
    Real floor_u_max = -1.0;  // B: |u| cap for soft floor hits (Lua `floor_u_max`; <= 0 disables)
    Real T_lo_all = 0.0, T_hi_all = 0.0;            // full-component T-hull intersection

    MixRule rule = MixRule::Dalton;
    Real ttol = 1.0e-10;
    int max_newton = 100;
    Real pure_tol = 1.0e-10;  // alpha_k >= 1 - pure_tol -> single-table path
    Real drop_tol = 1.0e-10;  // alpha_k <= drop_tol     -> dropped from the sums

    // mutable per-call scratch + accounting. Safe as plain members for the
    // same reason as TabulatedEOS::n_hull_clamps: every evaluation path is
    // const and Cerberus runs the cell loops serially (OpenMP is
    // force-disabled) — this avoids a heap allocation per face call.
    mutable Vector<Real> m_alpha;            // sanitized fractions (size N)
    mutable std::vector<Real> m_w;           // retained renormalised weights
    mutable std::vector<int> m_retained;     // indexes with m_w > 0
    mutable std::vector<Real> m_rho_k;       // clamped partial densities
    mutable std::vector<Real> m_scale_k;     // dilute-limit p/dpdT scales
    mutable std::vector<EosEval> m_evk;      // per-component evals (cached)
    mutable std::vector<char> m_pinned;      // Amagat: hull-pinned this eval
    mutable std::vector<int> m_pin_flag;     // Amagat: raw stk.flag per component
                                              // this eval (0/1/2); m_pinned is
                                              // just (flag != 0) — the W31
                                              // counters need the split between
                                              // flag==1 (true hull edge) and
                                              // flag==2 (Newton-exhausted
                                              // bracket midpoint)

    mutable long n_hull_clamps = 0;              // total-rho / aggregate clamps
    mutable std::vector<long> n_hull_clamps_k;   // per-component clamp accounting

    // ---- TEMPORARY DIAGNOSTIC (W31 investigation) ----------------------
    // amagat_assemble pin accounting. The question these answer: of the
    // assemblies that pin at all, how many are PARTIAL (some pinned, one or
    // more live) -- originally the only branch that truncated the dvdp sum
    // (pre pinned_compliance_fix_plan.md; Remedy A now retains the pinned
    // term instead) -- which component is the pinned one, and how stiff the
    // assembled mixture gets. rho*cs2 > 4 (= DBL_MAX*DBL_MIN) is the
    // threshold above which gam1 overflows once p reaches 0. Delete this
    // block and its use sites to revert.
    mutable long n_asm = 0;             // amagat_assemble calls
    mutable long n_asm_partial = 0;     // partially pinned (some retained, some not)
    mutable long n_asm_allpin = 0;      // all pinned
    mutable long n_asm_ovf = 0;         // partial AND rho*cs2 > 4
    mutable long n_asm_nonpos_p = 0;    // assembled with p <= 0
    mutable std::vector<long> n_pin_partial_k;  // per-component, partial only
    mutable std::vector<long> n_live_partial_k; // per-component, partial only
    // pinned_compliance_fix_plan.md Step 5: Remedy A / Remedy B accounting,
    // extending the block above rather than replacing it. n_edge_flag1/2
    // count every retained pinned component (partial or all-pin), split by
    // the §1 caveat (flag 1 = true hull edge, flag 2 = Newton-exhausted
    // midpoint) — both now feed the assembly via Remedy A, but they are not
    // the same table state. n_cap_* count Remedy B cap hits, separately per
    // cap, so a cap that fires on a supposedly-healthy cell is visible
    // rather than silent (a self-test gate on these, asserting zero on
    // unpinned sweeps, is Step 7 — out of scope here).
    mutable long n_edge_flag1 = 0;      // retained pin, flag==1 (hull edge)
    mutable long n_edge_flag2 = 0;      // retained pin, flag==2 (midpoint)
    mutable long n_cap_cs2 = 0;         // cap 1 hits: cs2 <= max_k cs_k^2
    mutable long n_cap_K = 0;           // cap 2 hits: rho*cs2 <= max_k K_k
    mutable long n_cap_dpdrho = 0;      // dpdrho <= max_k dpdrho_k hits
    mutable Real max_rcs2_partial = 0.0;
    mutable Real max_rcs2_nopin = 0.0;
    mutable std::vector<Real> worst_w;  // weights at the worst partial event
    mutable Real worst_p = 0.0, worst_rho = 0.0, worst_T = 0.0;
    mutable long n_alpha_fixes = 0;              // sanitize_alpha repairs
    mutable long n_mom_fixes = 0;                // B: floor-time momentum repairs
    mutable bool last_repair = false;            // B: last cons2prim repaired U in place
};

#endif  // MFP_MIXTURE_GAS_H
