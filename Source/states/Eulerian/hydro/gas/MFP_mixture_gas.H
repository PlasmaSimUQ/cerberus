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
//       mixing_rule = 'dalton',            -- 'amagat' aborts until plan W27
//       components  = {
//           { table = 'casing.eostab', name = 'casing', mass = 1.0, charge = 0.0 },
//           { table = 'fuel.eostab',   name = 'fuel',   mass = 0.5, charge = 0.0 },
//       },
//       ttol = 1e-10, max_newton = 100,    -- inversion controls (optional)
//       pure_tol = 1e-10, drop_tol = 1e-10,
//   }
//   value = { ..., alpha = alpha_casing }  -- N-1 tracers; last = 1 - sum
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
// (W27) is the production closure for condensed casing+fuel problems; the
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

    // union / conservative scalars derived from the caches (ctor)
    Real rho_union_min = 0.0, rho_union_max = 0.0;  // min_k / max_k of the rho hulls
    Real rho_floor = 0.0, p_floor = 0.0;            // conservative max_k floors (plan 5.6)
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

    mutable long n_hull_clamps = 0;              // total-rho / aggregate clamps
    mutable std::vector<long> n_hull_clamps_k;   // per-component clamp accounting
    mutable long n_alpha_fixes = 0;              // sanitize_alpha repairs
};

#endif  // MFP_MIXTURE_GAS_H
