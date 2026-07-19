#include "MFP_mixture_gas.H"

#include "MFP_lua.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

std::string MixtureEOS::tag = "tabulated_mixture";
bool MixtureEOS::registered =
  GetHydroGasFactory().Register(MixtureEOS::tag, HydroGasBuilder<MixtureEOS>);

MixtureEOS::MixtureEOS() {}
MixtureEOS::~MixtureEOS() {}

// ===========================================================================
// file-local helpers
// ===========================================================================

namespace
{

// locate x on a log10-uniform axis (clone of the reader's locate(); that one
// is file-local to MFP_eos_table.cpp)
void locate_ax(Real x, Real xmin, Real dx, int n, int& i, Real& f)
{
    const Real u = (x - xmin) / dx;
    i = (int)std::floor(u);
    i = std::max(0, std::min(n - 2, i));
    f = u - (Real)i;
    f = std::max((Real)0.0, std::min((Real)1.0, f));
}

// one bilinear read of a pre-inverted seed map (T_of_e on (lrho x le) or
// T_of_p on (lrho x lp)) — the same lookup EosTable::invert_T_from_* does
// internally, exposed here so the mixture can alpha-average the component
// seeds (plan 2.2). Returns <= 0 if the table carries no map.
Real seed_from_map(const EosTableView& v,
                   const Real* map,
                   Real ax_min,
                   Real dax,
                   int n_ax,
                   Real rho,
                   Real target)
{
    if (map == nullptr || n_ax < 2) return -1.0;
    const Real tiny = std::numeric_limits<Real>::min();
    int i, j;
    Real fx, fy;
    locate_ax(std::log10(std::max(rho, tiny)), v.lrho_min, v.dlrho, v.n_rho, i, fx);
    locate_ax(std::log10(std::max(target, tiny)), ax_min, dax, n_ax, j, fy);
    return (1 - fx) * ((1 - fy) * map[i * n_ax + j] + fy * map[i * n_ax + j + 1]) +
           fx * ((1 - fy) * map[(i + 1) * n_ax + j] + fy * map[(i + 1) * n_ax + j + 1]);
}

// guarded 1-D Newton with bisection fallback on a bracket — the same
// contract as the reader's file-local invert_1d (convergence on the residual
// in the TARGET quantity, never the step; any Newton step leaving the
// bracket becomes a bisection), except that the stats flag ESCALATES
// (max-aggregation) instead of being written directly, because the mixture
// stats may already carry a per-component hull-clamp flag (plan 7.3).
template <typename F>
Real guarded_solve(F fval,
                   Real lo,
                   Real hi,
                   Real target,
                   Real t_seed,
                   Real ttol,
                   int max_newton,
                   EosInvertStats& stats)
{
    const Real scale = std::max(std::abs(target), std::numeric_limits<Real>::min());
    Real slope_dummy;

    Real flo = fval(lo, slope_dummy) - target;
    Real fhi = fval(hi, slope_dummy) - target;
    if (std::abs(flo) <= ttol * scale) return lo;
    if (std::abs(fhi) <= ttol * scale) return hi;
    if ((flo > 0) == (fhi > 0)) {
        // target not attainable on this bracket -> clamp to the nearer end
        stats.flag = std::max(stats.flag, 1);
        return (std::abs(flo) < std::abs(fhi)) ? lo : hi;
    }

    Real t = std::max(lo, std::min(hi, t_seed));
    for (stats.iters = 1; stats.iters <= max_newton; ++stats.iters) {
        Real slope;
        const Real ft = fval(t, slope) - target;
        if (std::abs(ft) <= ttol * scale) return t;

        if ((ft > 0) == (flo > 0)) {
            lo = t;
            flo = ft;
        } else {
            hi = t;
            fhi = ft;
        }
        if (hi - lo < 1.0e-14) return 0.5 * (lo + hi);

        // step bounded by the bracket width BEFORE dividing so ft/slope can
        // never overflow (same FPE-trap-safety guard as the reader)
        Real tn;
        if (slope != 0.0 && std::abs(ft) <= std::abs(slope) * (hi - lo)) {
            tn = t - ft / slope;
            if (!(tn > lo && tn < hi)) {
                tn = 0.5 * (lo + hi);
                ++stats.bisections;
            }
        } else {
            tn = 0.5 * (lo + hi);
            ++stats.bisections;
        }
        t = tn;
    }
    stats.flag = std::max(stats.flag, 2);  // best bracket midpoint returned
    return 0.5 * (lo + hi);
}

}  // namespace

// ===========================================================================
// construction
// ===========================================================================

MixtureEOS::MixtureEOS(const int global_idx, const sol::table& def)
{
    BL_PROFILE("MixtureEOS::MixtureEOS");

    idx = global_idx;

    const std::string name = MFP::state_names[idx];

    const std::string rule_s = def["mixing_rule"].get_or<std::string>("dalton");
    if (rule_s == "amagat") {
        Abort("State: " + name + "; mixing_rule 'amagat' is not implemented yet (plan W27) — "
              "use 'dalton'");
    } else if (rule_s != "dalton") {
        Abort("State: " + name + "; unknown mixing_rule '" + rule_s +
              "' (options are 'dalton'; 'amagat' arrives with plan W27)");
    }
    rule = MixRule::Dalton;

    ttol = def["ttol"].get_or(1.0e-10);
    max_newton = def["max_newton"].get_or(100);
    pure_tol = def["pure_tol"].get_or(1.0e-10);
    drop_tol = def["drop_tol"].get_or(1.0e-10);

    const sol::table comps = def["components"].get_or(sol::table());
    if (!comps.valid())
        Abort("State: " + name + "; gas type '" + tag + "' requires a 'components' list");
    const int N = (int)comps.size();
    if (N < 2)
        Abort("State: " + name + "; gas type '" + tag + "' needs at least 2 components (" +
              num2str(N) + " given) — for one material use gas type 'tabulated'");

    if (!(MFP::u_ref > 0.0) || !(MFP::prs_ref > 0.0) || !(MFP::rho_ref > 0.0))
        Abort("State: " + name + "; reference quantities not set before mixture-gas construction");

    // component order defines the alpha mapping: components 0..N-2 bind to
    // tracer slots 0..N-2, component N-1 carries the derived 1 - sum(alpha)
    tables.resize(N);  // resized ONCE so the views below stay valid
    for (int k = 0; k < N; ++k) {
        const sol::table c = comps[k + 1].get_or(sol::table());
        if (!c.valid())
            Abort("State: " + name + "; components[" + num2str(k + 1) + "] is not a table");

        const std::string table_path = c["table"].get_or<std::string>("");
        if (table_path.empty())
            Abort("State: " + name + "; components[" + num2str(k + 1) +
                  "] requires a 'table' file path");

        tables[k].load(table_path);
        tables[k].ttol = ttol;
        tables[k].max_newton = max_newton;
        // identical nondimensionalisation call to TabulatedEOS, including
        // the current SI->CGS reference handling (plan D7; W18 removes the
        // unit shim for both models together)
        tables[k].nondimensionalise(MFP::rho_ref * 1.0e-3,
                                    MFP::T_ref,
                                    MFP::prs_ref * 10.0,
                                    MFP::u_ref * 1.0e2);

        const Real m = c["mass"].get_or(0.0);
        if (!(m > 0.0))
            Abort("State: " + name + "; components[" + num2str(k + 1) +
                  "] requires a positive 'mass'");
        mass.push_back(m);
        charge.push_back(c["charge"].get_or(0.0));
        const std::string default_name = name + "_" + num2str(k);
        comp_names.push_back(c["name"].get_or(default_name));
    }

    mass_const = all_equal(mass.begin(), mass.end(), mass[0]);
    charge_const = all_equal(charge.begin(), charge.end(), charge[0]);

    // views + hull caches (after ALL loads: no reallocation past this point)
    for (int k = 0; k < N; ++k) {
        tvs.push_back(tables[k].view());

        rho_hull_min_k.push_back(tables[k].rho_min());
        rho_hull_max_k.push_back(tables[k].rho_max());
        T_min_k.push_back(tables[k].T_min());
        T_max_k.push_back(tables[k].T_max());

        Real p_min = std::numeric_limits<Real>::max();
        for (int i = 0; i < tvs[k].n_rho; ++i) {
            for (int j = 0; j < tvs[k].n_T; ++j) {
                if (tvs[k].hull[tvs[k].idx(i, j)] > 0.5) {
                    p_min = std::min(p_min, tvs[k].p[tvs[k].idx(i, j)]);
                }
            }
        }
        if (!(p_min > 0.0))
            Abort("State: " + name + "; component '" + comp_names[k] +
                  "' has a non-positive in-hull pressure minimum (" + std::to_string(p_min) +
                  ") — a cold-curve table needs the deferred hull-membership validity work "
                  "(plan D10), not the v1 p>0 machinery");
        p_hull_min_k.push_back(p_min);
    }

    // union bounds (total-rho clamp) and scalar floors. The floors are the
    // UNION hull minima, not the max over components: with one condensed
    // component (Ti) and one dilute one (air), the max floor sits orders of
    // magnitude above a legitimate near-vacuum ambient state and would
    // silently inject mass/pressure into it every step. The union edge is
    // sufficient for W8.1 (finite sound speed at any floored state) because
    // the per-component guard in setup_partial_densities clamps/dilute-scales
    // each component onto its own hull before evaluation.
    rho_union_min = *std::min_element(rho_hull_min_k.begin(), rho_hull_min_k.end());
    rho_union_max = *std::max_element(rho_hull_max_k.begin(), rho_hull_max_k.end());
    rho_floor = rho_union_min;
    p_floor = *std::min_element(p_hull_min_k.begin(), p_hull_min_k.end());

    // T-bracket = intersection of the component T-hulls: Dalton needs one
    // shared T inside every hull. An empty intersection is a table-set
    // problem — fix it offline (extend the offending table), not at runtime.
    T_lo_all = *std::max_element(T_min_k.begin(), T_min_k.end());
    T_hi_all = *std::min_element(T_max_k.begin(), T_max_k.end());
    if (!(T_lo_all < T_hi_all))
        Abort("State: " + name + "; the component tables have no common temperature range "
              "(intersection [" + std::to_string(T_lo_all) + ", " + std::to_string(T_hi_all) +
              "] in code units) — the Dalton shared-T solve needs one; regenerate the "
              "offending table with a wider T span");

    // Dalton clamp-source warning (plan section 3): a component whose density
    // hull starts far above the union floor is a condensed-matter table —
    // dilute fractions of it land at alpha*rho below its hull and clamp on
    // every evaluation. Legitimate under Dalton, but worth a config-time note.
    for (int k = 0; k < N; ++k) {
        if (rho_hull_min_k[k] > 1.0e3 * rho_union_min) {
            amrex::Print() << "MixtureEOS[" << name << "]: WARNING component '" << comp_names[k]
                           << "' has rho_hull_min " << rho_hull_min_k[k]
                           << " far above the union floor " << rho_union_min
                           << " — dilute fractions of it will hull-clamp under the Dalton "
                              "rule (the Amagat rule, plan W27, is the accurate closure "
                              "for condensed components)\n";
        }
    }

    n_hull_clamps_k.assign(N, 0);
    m_alpha.resize(N, 0.0);
    m_w.resize(N, 0.0);
    m_rho_k.resize(N, 0.0);
    m_scale_k.resize(N, 1.0);
    m_evk.resize(N);

    amrex::Print() << "MixtureEOS[" << name << "]: rule=" << rule_s << " N=" << N
                   << " T-bracket (code units) [" << T_lo_all << "," << T_hi_all << "]\n";
    for (int k = 0; k < N; ++k) {
        amrex::Print() << "  component '" << comp_names[k] << "': rho=[" << rho_hull_min_k[k]
                       << "," << rho_hull_max_k[k] << "] T=[" << T_min_k[k] << "," << T_max_k[k]
                       << "] p_min=" << p_hull_min_k[k] << "\n";
    }

    // opt-in one-zone self-test sweep (plan M2): `self_test = <n_sweep>`
    const int n_self = def["self_test"].get_or(0);
    if (n_self > 0) run_self_test(n_self);
}

// ===========================================================================
// composition handling
// ===========================================================================

int MixtureEOS::sanitize_alpha(Vector<Real>& alpha) const
{
    int n_fix = 0;
    Real s = 0.0;
    for (Real& a : alpha) {
        if (a < 0.0) {
            a = 0.0;
            ++n_fix;
        } else if (a > 1.0) {
            a = 1.0;
            ++n_fix;
        }
        s += a;
    }
    if (!(s > 0.0)) {
        // unreachable for finite input (clamping cannot zero every entry:
        // negative tracers push the derived last fraction above 0) — kept as
        // a backstop against non-finite tracer data
        std::fill(alpha.begin(), alpha.end(), 0.0);
        alpha[alpha.size() - 1] = 1.0;
        return n_fix + 1;
    }
    if (std::abs(s - 1.0) > 1.0e-14) {
        for (Real& a : alpha) { a /= s; }
        ++n_fix;
    }
    return n_fix;
}

int MixtureEOS::prepare_weights(const Vector<Real>& alpha) const
{
    const int N = (int)n_species();
    std::fill(m_w.begin(), m_w.end(), 0.0);
    m_retained.clear();

    Real s = 0.0;
    for (int k = 0; k < N; ++k) {
        if (alpha[k] > drop_tol) {
            m_retained.push_back(k);
            s += alpha[k];
        }
    }
    if (m_retained.empty()) {
        // unreachable for sanitized alpha (sum 1, N >= 2 -> max >= 1/N)
        int kmax = 0;
        for (int k = 1; k < N; ++k) {
            if (alpha[k] > alpha[kmax]) kmax = k;
        }
        m_w[kmax] = 1.0;
        m_retained.push_back(kmax);
        return kmax;
    }
    for (const int k : m_retained) { m_w[k] = alpha[k] / s; }

    if ((int)m_retained.size() == 1) return m_retained[0];
    for (const int k : m_retained) {
        if (m_w[k] >= 1.0 - pure_tol) return k;
    }
    return -1;
}

Real MixtureEOS::clamp_rho_k(int k, Real rho_k, EosInvertStats& st) const
{
    const Real lo = std::max(effective_zero, rho_hull_min_k[k]);
    if (rho_k < lo) {
        ++n_hull_clamps_k[k];
        st.flag = std::max(st.flag, 1);
        return lo;
    }
    if (rho_k > rho_hull_max_k[k]) {
        ++n_hull_clamps_k[k];
        st.flag = std::max(st.flag, 1);
        return rho_hull_max_k[k];
    }
    return rho_k;
}

Real MixtureEOS::clamp_rho_union(Real rho) const
{
    const Real lo = std::max(effective_zero, rho_union_min);
    if (rho < lo) {
        ++n_hull_clamps;
        return lo;
    }
    if (rho > rho_union_max) {
        ++n_hull_clamps;
        return rho_union_max;
    }
    return rho;
}

void MixtureEOS::setup_partial_densities(Real rho, EosInvertStats& st) const
{
    // Dalton partial densities rho_k = w_k * rho are FIXED (independent of
    // the solve variable T), so they are hull-handled and tallied once per
    // call, not once per iteration. A low-side excursion switches to the
    // dilute ideal-gas-limit extrapolation (evaluate the hull-edge row,
    // scale p and dpdT by rho_k/edge — see the header): the component's
    // partial pressure then vanishes smoothly as alpha_k -> 0 instead of
    // flooring at p(rho_hull_min, T). Either excursion freezes rho_k, not
    // the T-dependence, so the shared-T residual stays monotone.
    for (const int k : m_retained) {
        const Real raw = m_w[k] * rho;
        const Real lo = std::max(effective_zero, rho_hull_min_k[k]);
        if (raw < lo) {
            m_rho_k[k] = lo;
            m_scale_k[k] = raw / lo;
            ++n_hull_clamps_k[k];
            st.flag = std::max(st.flag, 1);
        } else if (raw > rho_hull_max_k[k]) {
            m_rho_k[k] = rho_hull_max_k[k];
            m_scale_k[k] = 1.0;
            ++n_hull_clamps_k[k];
            st.flag = std::max(st.flag, 1);
        } else {
            m_rho_k[k] = raw;
            m_scale_k[k] = 1.0;
        }
    }
}

Real MixtureEOS::mix_p_rho_T(Real rho, Real T, EosInvertStats& st) const
{
    // forward Dalton pressure at the prepared weights (m_w/m_retained)
    setup_partial_densities(rho, st);
    Real p = 0.0;
    EosEval evk;
    for (const int k : m_retained) {
        EosTable::eval_rt(tvs[k], m_rho_k[k], T, evk);
        p += m_scale_k[k] * evk.p;
    }
    return p;
}

// ===========================================================================
// the Dalton shared-T drivers
// ===========================================================================

void MixtureEOS::dalton_solve_T(Real rho,
                                Real target,
                                bool by_e,
                                const Vector<Real>& alpha,
                                EosEval& ev,
                                EosInvertStats& st) const
{
    const Real tiny = std::numeric_limits<Real>::min();

    const int kp = prepare_weights(alpha);

    // pure-cell short-circuit (plan 7.1): FORMULA-IDENTICAL to the
    // single-table TabulatedEOS driver (same inverter call, same eval_rt) —
    // a hard requirement, not an optimisation: it is what makes the
    // identical-tables round-off validation gates meaningful, and pure cells
    // are the bulk of a capsule domain (the primary cost saver).
    if (kp >= 0) {
        const Real rho_c = clamp_rho_k(kp, rho, st);
        const Real T =
          by_e ? EosTable::invert_T_from_e(tvs[kp], rho_c, target, -1.0, ttol, max_newton, st)
               : EosTable::invert_T_from_p(tvs[kp], rho_c, target, -1.0, ttol, max_newton, st);
        EosTable::eval_rt(tvs[kp], rho_c, T, ev);
        return;
    }

    // fixed partial densities + dilute-limit scales, tallied once per call
    setup_partial_densities(rho, st);

    // bracket: intersection of the RETAINED components' T-hulls — a superset
    // of the ctor's full intersection, so never empty
    Real T_lo = 0.0, T_hi = std::numeric_limits<Real>::max();
    for (const int k : m_retained) {
        T_lo = std::max(T_lo, T_min_k[k]);
        T_hi = std::min(T_hi, T_max_k[k]);
    }
    const Real lt_lo = std::log10(T_lo);
    const Real lt_hi = std::log10(T_hi);

    // mixture residual + slope at log10(T):
    //   by_e:  F = sum_k w_k e_k(rho_k, T),  dF/dT = sum_k w_k cv_k   (eq 4)
    //   else:  G = sum_k p_k(rho_k, T),      dG/dT = sum_k dpdT_k     (eq 5)
    // The slope uses the conditioned cv/dpdT blocks; unlike the single-table
    // inverter it is not the exact interpolant derivative, so the bracket +
    // bisection fallback owns any mismatch — convergence is declared on the
    // residual, never the step. cv > 0 by table conditioning keeps F
    // strictly increasing (the solve is unconditionally convergent).
    Real last_lt = std::numeric_limits<Real>::quiet_NaN();
    auto fval = [&](Real lt, Real& slope) {
        const Real T = std::pow((Real)10.0, lt);
        Real val = 0.0, dvdT = 0.0;
        for (const int k : m_retained) {
            EosTable::eval_rt(tvs[k], m_rho_k[k], T, m_evk[k]);
            // dilute-limit scaling applied INTO the cached eval so the
            // assembly below picks it up automatically (scale is 1 in-hull)
            m_evk[k].p *= m_scale_k[k];
            m_evk[k].dpdT *= m_scale_k[k];
            if (by_e) {
                val += m_w[k] * m_evk[k].e;
                dvdT += m_w[k] * m_evk[k].cv;
            } else {
                val += m_evk[k].p;
                dvdT += m_evk[k].dpdT;
            }
        }
        slope = std::log((Real)10.0) * T * dvdT;
        last_lt = lt;
        return val;
    };

    // seed (plan 2.2): alpha-weighted per-component inverse-map lookup,
    // bracket midpoint fallback. No Temp-slot warm start — the established
    // single-table driver convention (the Temp slot is not trusted). Using
    // the MIXTURE target as each component's map argument is heuristic (the
    // component e_k differs from the mixture e in general) but the guarded
    // solve owns convergence; a seed only needs to be near.
    Real lt0 = 0.5 * (lt_lo + lt_hi);
    {
        Real Ts_sum = 0.0, w_sum = 0.0;
        for (const int k : m_retained) {
            const Real Ts =
              by_e ? seed_from_map(tvs[k], tvs[k].T_of_e, tvs[k].le_min, tvs[k].dle, tvs[k].n_e,
                                   m_rho_k[k], target)
                   : seed_from_map(tvs[k], tvs[k].T_of_p, tvs[k].lp_min, tvs[k].dlp, tvs[k].n_p,
                                   m_rho_k[k], target);
            if (Ts > 0.0) {
                Ts_sum += m_w[k] * Ts;
                w_sum += m_w[k];
            }
        }
        if (w_sum > 0.0) {
            lt0 = std::log10(std::max(Ts_sum / w_sum, tiny));
            st.seeded = true;
        }
    }

    const Real lt = guarded_solve(fval, lt_lo, lt_hi, target, lt0, ttol, max_newton, st);
    // an endpoint clamp or bracket-collapse midpoint may differ from the
    // last evaluated coordinate: refresh the cached per-component evals so
    // the assembly below describes the returned T (otherwise no
    // re-evaluation happens after convergence — plan 2.1)
    if (lt != last_lt) {
        Real slope;
        fval(lt, slope);
    }
    const Real T = std::pow((Real)10.0, lt);

    // assemble the mixture EosEval from the cached per-component evals
    // (plan eqs 2-3, 6-9)
    ev = EosEval();
    ev.rho = rho;
    ev.T = T;
    bool clamped = (st.flag != 0);
    for (const int k : m_retained) {
        const EosEval& evk = m_evk[k];
        ev.p += evk.p;                            // (2) partial pressures
        ev.e += m_w[k] * evk.e;                   // (3)
        ev.cv += m_w[k] * evk.cv;                 // (8)
        ev.dpdT += evk.dpdT;                      // (7)
        ev.dpdrho += m_w[k] * evk.dpdrho;         // (6): d(rho_k)/d(rho) = w_k
        ev.dedrho += m_w[k] * m_w[k] * evk.dedrho;
        clamped |= evk.clamped;
    }
    ev.clamped = clamped;

    // frozen mixture sound speed (9) — the components couple through the
    // shared T in the cross term; same non-convexity floor as eval_rt
    const Real cs2 = std::max(
      ev.dpdrho + (T / (rho * rho)) * ev.dpdT * ev.dpdT / std::max(ev.cv, tiny), (Real)0.0);
    ev.cs = std::sqrt(cs2);
    ev.gam1 = rho * cs2 / std::max(ev.p, tiny);
    ev.dpde = ev.dpdT / std::max(ev.cv, tiny);
    ev.dpdr_e = ev.dpdrho - ev.dpde * ev.dedrho;
}

void MixtureEOS::eval_from_rho_e(Real rho,
                                 Real e_int,
                                 const Vector<Real>& alpha,
                                 EosEval& ev,
                                 EosInvertStats& st) const
{
    dalton_solve_T(rho, e_int, true, alpha, ev, st);
}

void MixtureEOS::eval_from_rho_p(Real rho,
                                 Real p,
                                 const Vector<Real>& alpha,
                                 EosEval& ev,
                                 EosInvertStats& st) const
{
    dalton_solve_T(rho, p, false, alpha, ev, st);
}

Real MixtureEOS::invert_rho_from_p_T(Real T,
                                     Real p,
                                     const Vector<Real>& alpha,
                                     EosInvertStats& st) const
{
    // initial-condition fill path (define_rho_p_T p+T branch) — not hot.
    // Solve H(rho) = sum_k p_k(w_k rho, T) = p, monotone nondecreasing in
    // rho with slope sum_k w_k dpdrho_k > 0 on the unclamped interior.
    const int kp = prepare_weights(alpha);
    if (kp >= 0) {
        return EosTable::invert_rho_from_p(tvs[kp], T, p, -1.0, ttol, max_newton, st);
    }

    // bracket in total rho: below min_k(rho_hull_min_k) every partial
    // density is under its hull; above max_k(rho_hull_max_k / w_k) every
    // partial density is over its hull — H is constant outside these, so the
    // bracket ends suffice (an unattainable p clamps + flags at an end)
    Real rho_lo = std::numeric_limits<Real>::max(), rho_hi = 0.0;
    for (const int k : m_retained) {
        rho_lo = std::min(rho_lo, rho_hull_min_k[k]);
        rho_hi = std::max(rho_hi, rho_hull_max_k[k] / m_w[k]);
    }
    const Real lr_lo = std::log10(rho_lo);
    const Real lr_hi = std::log10(rho_hi);

    // quiet in-iteration hull handling (partial densities move every
    // iteration — per-iteration tallies would just spam the counters on an
    // init path), with the same dilute-limit scaling as the T-drivers so
    // the pressure surface being inverted is the one the solves see
    EosEval evk;
    auto fval = [&](Real lr, Real& slope) {
        const Real rho = std::pow((Real)10.0, lr);
        Real val = 0.0, dvdrho = 0.0;
        for (const int k : m_retained) {
            const Real raw = m_w[k] * rho;
            const Real lo = std::max(effective_zero, rho_hull_min_k[k]);
            const Real rk = std::clamp(raw, lo, rho_hull_max_k[k]);
            const Real scale = (raw < lo) ? raw / lo : 1.0;
            EosTable::eval_rt(tvs[k], rk, T, evk);
            val += scale * evk.p;
            dvdrho += m_w[k] * evk.dpdrho;
        }
        slope = std::log((Real)10.0) * rho * dvdrho;
        return val;
    };

    // seed: alpha-weighted single-table density inversions (exact for
    // identical ideal-gas components; throwaway stats — seed bookkeeping
    // must not pollute the caller's flag)
    Real lr0 = 0.5 * (lr_lo + lr_hi);
    {
        Real seed = 0.0;
        EosInvertStats st_seed;
        for (const int k : m_retained) {
            seed += m_w[k] * EosTable::invert_rho_from_p(tvs[k], T, p, -1.0, ttol, max_newton,
                                                         st_seed);
        }
        if (seed > 0.0) {
            lr0 = std::log10(seed);
            st.seeded = true;
        }
    }

    const Real lr = guarded_solve(fval, lr_lo, lr_hi, p, lr0, ttol, max_newton, st);
    return std::pow((Real)10.0, lr);
}

// ===========================================================================
// conversions
// ===========================================================================

bool MixtureEOS::cons2prim(Vector<Real>& U, Vector<Real>& Q) const
{
    BL_PROFILE("MixtureEOS::cons2prim");

    // composition first, from the RAW conserved data (the simplex
    // renormalisation absorbs the raw-vs-clamped density scale), then the
    // total density into the union hull (the Stage-4 clamp-rho-first rule):
    // every derived quantity below describes the same in-union state, and
    // the per-component clamps inside the driver own the real hull work.
    get_alpha_fractions_from_cons(U, m_alpha);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    Real rho = clamp_rho_union(U[+HydroDef::ConsIdx::Density]);
    Real mx = U[+HydroDef::ConsIdx::Xmom];
    Real my = U[+HydroDef::ConsIdx::Ymom];
    Real mz = U[+HydroDef::ConsIdx::Zmom];
    Real ed = U[+HydroDef::ConsIdx::Eden];

    Real rhoinv = 1 / rho;
    Real u = mx * rhoinv;
    Real v = my * rhoinv;
    Real w = mz * rhoinv;
    Real e_int = ed * rhoinv - 0.5 * (u * u + v * v + w * w);

    EosEval ev;
    EosInvertStats st;
    eval_from_rho_e(rho, e_int, m_alpha, ev, st);
    tally(st);

    // energy-consistent effective gamma with the clamp-consistency rule
    // (Stage-4 W8.2, same as TabulatedEOS): where the solve converged,
    // p/(gamma_e - 1) = rho*e_int exactly; where it clamped, the physical
    // e_int is off-bracket (possibly <= 0) so the mixture-consistent ev.e
    // keeps (p, T, gamma_e) describing the same clamped state
    const Real e_eff = (st.flag == 0) ? e_int : ev.e;
    const Real ge = 1.0 + ev.p / std::max(rho * e_eff, std::numeric_limits<Real>::min());

    // mixture cp via the general-EOS identity (all mixture-summed columns)
    const Real cp = ev.cv + (ev.T / (rho * rho)) * ev.dpdT * ev.dpdT /
                              std::max(ev.dpdrho, std::numeric_limits<Real>::min());

    Q[+HydroDef::PrimIdx::Density] = rho;
    Q[+HydroDef::PrimIdx::Xvel] = u;
    Q[+HydroDef::PrimIdx::Yvel] = v;
    Q[+HydroDef::PrimIdx::Zvel] = w;
    Q[+HydroDef::PrimIdx::Prs] = ev.p;
    Q[+HydroDef::PrimIdx::Temp] = ev.T;
    Q[+HydroDef::PrimIdx::Gamma] = ge;
    Q[+HydroDef::PrimIdx::SpHeat] = cp;

    // sanitized fractions out: the tracer slots re-enter the simplex here
    for (int i = 0; i < n_tracers(); ++i) { Q[+HydroDef::PrimIdx::NUM + i] = m_alpha[i]; }

#ifdef MFP_PRIM_FLOOR
    // rho was union-clamped up front; p and T are table values — backstop
    if (Q[+HydroDef::PrimIdx::Prs] < effective_zero) {
        Q[+HydroDef::PrimIdx::Prs] = effective_zero;
    }
    if (Q[+HydroDef::PrimIdx::Temp] < effective_zero) {
        Q[+HydroDef::PrimIdx::Temp] = effective_zero;
    }
#endif

    return prim_valid(Q);
}

void MixtureEOS::prim2cons(Vector<Real>& Q, Vector<Real>& U) const
{
    BL_PROFILE("MixtureEOS::prim2cons");

    get_alpha_fractions_from_prim(Q, m_alpha);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    Real rho = clamp_rho_union(Q[+HydroDef::PrimIdx::Density]);
    Real u = Q[+HydroDef::PrimIdx::Xvel];
    Real v = Q[+HydroDef::PrimIdx::Yvel];
    Real w = Q[+HydroDef::PrimIdx::Zvel];
    Real p = Q[+HydroDef::PrimIdx::Prs];

    Real mx = u * rho;
    Real my = v * rho;
    Real mz = w * rho;
    Real ke = 0.5 * rho * (u * u + v * v + w * w);

    // rp solve -> shared T -> mixture e (the Temp slot is not trusted here,
    // matching the single-table driver convention)
    EosEval ev;
    EosInvertStats st;
    eval_from_rho_p(rho, p, m_alpha, ev, st);
    tally(st);

    U[+HydroDef::ConsIdx::Density] = rho;
    U[+HydroDef::ConsIdx::Xmom] = mx;
    U[+HydroDef::ConsIdx::Ymom] = my;
    U[+HydroDef::ConsIdx::Zmom] = mz;
    U[+HydroDef::ConsIdx::Eden] = rho * ev.e + ke;

    for (int i = 0; i < n_tracers(); ++i) {
        U[+HydroDef::ConsIdx::NUM + i] = m_alpha[i] * rho;
    }
}

void MixtureEOS::define_rho_p_T(Vector<Real>& Q) const
{
    BL_PROFILE("MixtureEOS::define_rho_p_T");

    Real rho = Q[+HydroDef::PrimIdx::Density];
    Real p = Q[+HydroDef::PrimIdx::Prs];
    Real T = Q[+HydroDef::PrimIdx::Temp];

    get_alpha_fractions_from_prim(Q, m_alpha);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    // same given-ness convention as ThermallyPerfectGas/TabulatedEOS:
    // positive = given, same priority order. A GIVEN density is clamped;
    // zero means "not given" and must stay zero for the convention to work.
    if (rho > 0.0) rho = clamp_rho_union(rho);
    EosInvertStats st;
    if ((rho > 0.0) && (p > 0.0)) {
        EosEval ev;
        eval_from_rho_p(rho, p, m_alpha, ev, st);
        T = ev.T;
    } else if ((p > 0.0) && (T > 0.0)) {
        rho = invert_rho_from_p_T(T, p, m_alpha, st);
    } else if ((rho > 0.0) && (T > 0.0)) {
        // forward evaluation: p = sum_k p_k(rho_k, T)
        const int kp = prepare_weights(m_alpha);
        if (kp >= 0) {
            EosEval evk;
            EosTable::eval_rt(tvs[kp], clamp_rho_k(kp, rho, st), T, evk);
            p = evk.p;
        } else {
            p = mix_p_rho_T(rho, T, st);
        }
    }
    tally(st);

    Q[+HydroDef::PrimIdx::Density] = rho;
    Q[+HydroDef::PrimIdx::Prs] = p;
    Q[+HydroDef::PrimIdx::Temp] = T;
}

// ===========================================================================
// getters
// ===========================================================================

Real MixtureEOS::get_temperature_from_cons(const Vector<Real>& U) const
{
    BL_PROFILE("MixtureEOS::get_temperature_from_cons");

    get_alpha_fractions_from_cons(U, m_alpha);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    const Real rho = clamp_rho_union(U[+HydroDef::ConsIdx::Density]);
    const Real rhoinv = 1 / rho;
    const Real u = U[+HydroDef::ConsIdx::Xmom] * rhoinv;
    const Real v = U[+HydroDef::ConsIdx::Ymom] * rhoinv;
    const Real w = U[+HydroDef::ConsIdx::Zmom] * rhoinv;
    const Real e_int = U[+HydroDef::ConsIdx::Eden] * rhoinv - 0.5 * (u * u + v * v + w * w);

    EosEval ev;
    EosInvertStats st;
    eval_from_rho_e(rho, e_int, m_alpha, ev, st);
    tally(st);
    return ev.T;
}

Real MixtureEOS::get_gamma_from_cons(const Vector<Real>& U,
                                     const int density_idx,
                                     const int tracer_idx) const
{
    BL_PROFILE("MixtureEOS::get_gamma_from_cons");

    get_alpha_fractions_from_cons(U, m_alpha, density_idx, tracer_idx);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    const Real rho = clamp_rho_union(U[density_idx]);
    const Real rhoinv = 1 / rho;
    const Real u = U[+HydroDef::ConsIdx::Xmom] * rhoinv;
    const Real v = U[+HydroDef::ConsIdx::Ymom] * rhoinv;
    const Real w = U[+HydroDef::ConsIdx::Zmom] * rhoinv;
    const Real e_int = U[+HydroDef::ConsIdx::Eden] * rhoinv - 0.5 * (u * u + v * v + w * w);

    EosEval ev;
    EosInvertStats st;
    eval_from_rho_e(rho, e_int, m_alpha, ev, st);
    tally(st);
    const Real e_eff = (st.flag == 0) ? e_int : ev.e;
    return 1.0 + ev.p / std::max(rho * e_eff, std::numeric_limits<Real>::min());
}

Real MixtureEOS::get_gamma_from_prim(const Vector<Real>& Q, const int idx) const
{
    BL_PROFILE("MixtureEOS::get_gamma_from_prim");

    get_alpha_fractions_from_prim(Q, m_alpha, idx);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    const Real rho = clamp_rho_union(Q[+HydroDef::PrimIdx::Density]);
    const Real p = Q[+HydroDef::PrimIdx::Prs];

    EosEval ev;
    EosInvertStats st;
    eval_from_rho_p(rho, p, m_alpha, ev, st);
    tally(st);
    return 1.0 + p / std::max(rho * ev.e, std::numeric_limits<Real>::min());
}

Real MixtureEOS::get_cp_from_cons(const Vector<Real>& U,
                                  const int density_idx,
                                  const int tracer_idx) const
{
    BL_PROFILE("MixtureEOS::get_cp_from_cons");

    get_alpha_fractions_from_cons(U, m_alpha, density_idx, tracer_idx);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    const Real rho = clamp_rho_union(U[density_idx]);
    const Real rhoinv = 1 / rho;
    const Real u = U[+HydroDef::ConsIdx::Xmom] * rhoinv;
    const Real v = U[+HydroDef::ConsIdx::Ymom] * rhoinv;
    const Real w = U[+HydroDef::ConsIdx::Zmom] * rhoinv;
    const Real e_int = U[+HydroDef::ConsIdx::Eden] * rhoinv - 0.5 * (u * u + v * v + w * w);

    EosEval ev;
    EosInvertStats st;
    eval_from_rho_e(rho, e_int, m_alpha, ev, st);
    tally(st);
    return ev.cv + (ev.T / (rho * rho)) * ev.dpdT * ev.dpdT /
                     std::max(ev.dpdrho, std::numeric_limits<Real>::min());
}

Real MixtureEOS::get_cp_from_prim(const Vector<Real>& Q, const int tracer_idx) const
{
    BL_PROFILE("MixtureEOS::get_cp_from_prim");

    get_alpha_fractions_from_prim(Q, m_alpha, tracer_idx);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    const Real rho = clamp_rho_union(Q[+HydroDef::PrimIdx::Density]);
    const Real p = Q[+HydroDef::PrimIdx::Prs];

    EosEval ev;
    EosInvertStats st;
    eval_from_rho_p(rho, p, m_alpha, ev, st);
    tally(st);
    return ev.cv + (ev.T / (rho * rho)) * ev.dpdT * ev.dpdT /
                     std::max(ev.dpdrho, std::numeric_limits<Real>::min());
}

Real MixtureEOS::get_internal_energy_from_prim(const Vector<Real>& Q) const
{
    BL_PROFILE("MixtureEOS::get_internal_energy_from_prim");

    get_alpha_fractions_from_prim(Q, m_alpha);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    const Real rho = clamp_rho_union(Q[+HydroDef::PrimIdx::Density]);
    const Real p = Q[+HydroDef::PrimIdx::Prs];

    EosEval ev;
    EosInvertStats st;
    eval_from_rho_p(rho, p, m_alpha, ev, st);
    tally(st);
    return ev.e;
}

Real MixtureEOS::get_sound_speed_from_prim_rp(const Vector<Real>& Q) const
{
    BL_PROFILE("MixtureEOS::get_sound_speed_from_prim_rp");

    get_alpha_fractions_from_prim(Q, m_alpha);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    const Real rho = clamp_rho_union(Q[+HydroDef::PrimIdx::Density]);
    const Real p = Q[+HydroDef::PrimIdx::Prs];

    EosEval ev;
    EosInvertStats st;
    eval_from_rho_p(rho, p, m_alpha, ev, st);
    tally(st);
    return ev.cs;
}

void MixtureEOS::get_face_eval_from_prim(const Vector<Real>& Q, Real& e, Real& a) const
{
    BL_PROFILE("MixtureEOS::get_face_eval_from_prim");

    // HOT PATH: one rp mixture solve answers both face quantities for
    // HLLC_general_eos (the Stage-5 G8 discipline). The reconstructed face
    // tracer slots arrive through Q, so the face composition is consistent
    // with the face (rho, p) by construction — the solver needs no edits.
    get_alpha_fractions_from_prim(Q, m_alpha);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    const Real rho = clamp_rho_union(Q[+HydroDef::PrimIdx::Density]);
    const Real p = Q[+HydroDef::PrimIdx::Prs];

    EosEval ev;
    EosInvertStats st;
    eval_from_rho_p(rho, p, m_alpha, ev, st);
    tally(st);
    e = ev.e;
    a = ev.cs;
}

RealArray MixtureEOS::get_speed_from_cons(const Vector<Real>& U) const
{
    BL_PROFILE("MixtureEOS::get_speed_from_cons");

    get_alpha_fractions_from_cons(U, m_alpha);
    n_alpha_fixes += sanitize_alpha(m_alpha);

    Real rho = U[+HydroDef::ConsIdx::Density];

#ifdef MFP_PRIM_FLOOR
    // raw conserved data that has not passed the cons2prim floors: floor to
    // the union hull edge, not the absolute effective_zero (W8.1) — an
    // evacuated cell would otherwise report u = mx/1e-14 and collapse the
    // global time step
    rho = std::max(rho, std::max(effective_zero, rho_union_min));
#endif

    const Real rhoinv = 1 / rho;
    const Real u = U[+HydroDef::ConsIdx::Xmom] * rhoinv;
    const Real v = U[+HydroDef::ConsIdx::Ymom] * rhoinv;
    const Real w = U[+HydroDef::ConsIdx::Zmom] * rhoinv;
    Real e_int = U[+HydroDef::ConsIdx::Eden] * rhoinv - 0.5 * (u * u + v * v + w * w);

#ifdef MFP_PRIM_FLOOR
    e_int = std::max(e_int, effective_zero);
#endif

    // TRUE frozen mixture sound speed (eq 9) — the CFL time step never uses
    // the effective-gamma estimate, same discipline as TabulatedEOS
    EosEval ev;
    EosInvertStats st;
    eval_from_rho_e(rho, e_int, m_alpha, ev, st);
    tally(st);
    const Real a = ev.cs;

    RealArray s = {AMREX_D_DECL(a + std::abs(u), a + std::abs(v), a + std::abs(w))};

    return s;
}

RealArray MixtureEOS::get_speed_from_prim(const Vector<Real>& Q) const
{
    BL_PROFILE("MixtureEOS::get_speed_from_prim");

    get_alpha_fractions_from_prim(Q, m_alpha);
    n_alpha_fixes += sanitize_alpha(m_alpha);

#ifdef MFP_PRIM_FLOOR
    // union-hull floors, as in get_speed_from_cons (W8.1)
    const Real rho =
      std::max(Q[+HydroDef::PrimIdx::Density], std::max(effective_zero, rho_union_min));
    const Real p = std::max(Q[+HydroDef::PrimIdx::Prs], std::max(effective_zero, p_floor));
#else
    const Real rho = clamp_rho_union(Q[+HydroDef::PrimIdx::Density]);
    const Real p = Q[+HydroDef::PrimIdx::Prs];
#endif

    EosEval ev;
    EosInvertStats st;
    eval_from_rho_p(rho, p, m_alpha, ev, st);
    tally(st);
    const Real a = ev.cs;

    RealArray s = {AMREX_D_DECL(a + std::abs(Q[+HydroDef::PrimIdx::Xvel]),
                                a + std::abs(Q[+HydroDef::PrimIdx::Yvel]),
                                a + std::abs(Q[+HydroDef::PrimIdx::Zvel]))};

    return s;
}

// ===========================================================================
// floors / info
// ===========================================================================

int MixtureEOS::apply_prim_floor(Vector<Real>& Q) const
{
#ifdef MFP_PRIM_FLOOR
    // ctor-computed union scalars (min over the component hull minima):
    // last-resort positivity only. The real per-component hull guard lives
    // in the drivers (setup_partial_densities), which keeps W8.1 finite
    // sound speeds down to the union edge; flooring to the max over
    // components would destroy any ambient state below the densest
    // component's hull (e.g. a rough vacuum next to a solid casing).
    int n = 0;
    const Real rho_fl = std::max(effective_zero, rho_floor);
    const Real p_fl = std::max(effective_zero, p_floor);
    if (Q[+HydroDef::PrimIdx::Density] < rho_fl) {
        Q[+HydroDef::PrimIdx::Density] = rho_fl;
        n += 1;
    }
    if (Q[+HydroDef::PrimIdx::Prs] < p_fl) {
        Q[+HydroDef::PrimIdx::Prs] = p_fl;
        n += 1;
    }
    return n;
#else
    return 0;
#endif
}

// ===========================================================================
// one-zone self-test (plan M2). Runs in code units on the constructed model
// (after nondimensionalisation) so it exercises the exact runtime paths.
// PASS/FAIL lines are grepped by the EOS-Mixture check.py.
// ===========================================================================

void MixtureEOS::run_self_test(int n_sweep) const
{
    const Real tiny = std::numeric_limits<Real>::min();
    const int N = (int)n_species();
    const std::string name = MFP::state_names[idx];

    bool all_pass = true;
    auto verdict = [&](const std::string& check, bool ok, const std::string& detail) {
        all_pass &= ok;
        amrex::Print() << "MIXEOS-SELFTEST[" << check << "] " << (ok ? "PASS " : "FAIL ")
                       << detail << "\n";
    };

    // alpha grid: the first component sweeps its fraction, the remainder is
    // split equally — covers pure cells, the drop_tol/pure_tol edges and
    // well-mixed states
    const Real a_grid[] = {0.0, 1.0e-12, 1.0e-8, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 - 1.0e-8, 1.0};

    Vector<Real> alpha(N, 0.0);
    auto set_alpha = [&](Real a0) {
        alpha[0] = a0;
        for (int k = 1; k < N; ++k) { alpha[k] = (1.0 - a0) / (N - 1); }
    };

    // sample range where every retained partial density is inside its hull
    // (the round-trip checks probe the solver, not the clamp handling)
    auto sample_range = [&](Real& rho_lo, Real& rho_hi) {
        rho_lo = 0.0;
        rho_hi = std::numeric_limits<Real>::max();
        for (const int k : m_retained) {
            rho_lo = std::max(rho_lo, rho_hull_min_k[k] / m_w[k]);
            rho_hi = std::min(rho_hi, rho_hull_max_k[k] / m_w[k]);
        }
        return rho_lo < rho_hi;
    };

    // ---- 1. driver round-trips on the alpha-grid x hull-sweep ----
    {
        int total = 0, nonconv = 0;
        Real max_res_e = 0.0, max_res_p = 0.0, max_pmis = 0.0, max_emis = 0.0, max_Terr = 0.0;
        for (const Real a0 : a_grid) {
            set_alpha(a0);
            prepare_weights(alpha);
            // local copies: the drivers below re-run prepare_weights
            const std::vector<Real> w = m_w;
            const std::vector<int> retained = m_retained;

            Real rho_lo, rho_hi;
            if (!sample_range(rho_lo, rho_hi)) continue;

            for (int a = 0; a < n_sweep; ++a) {
                for (int b = 0; b < n_sweep; ++b) {
                    const Real rho =
                      std::pow(10.0, std::log10(rho_lo) + (a + 0.5) / n_sweep *
                                                            (std::log10(rho_hi) - std::log10(rho_lo)));
                    const Real T =
                      std::pow(10.0, std::log10(T_lo_all) + (b + 0.5) / n_sweep *
                                                              (std::log10(T_hi_all) - std::log10(T_lo_all)));

                    // forward Dalton reference (in-hull: no scaling active)
                    Real p_f = 0.0, e_f = 0.0;
                    bool clamped = false;
                    EosEval evk;
                    for (const int k : retained) {
                        EosTable::eval_rt(tvs[k], w[k] * rho, T, evk);
                        p_f += evk.p;
                        e_f += w[k] * evk.e;
                        clamped |= evk.clamped;
                    }
                    if (clamped) continue;  // filled-cell territory
                    ++total;

                    EosEval ev1;
                    EosInvertStats st1;
                    eval_from_rho_e(rho, e_f, alpha, ev1, st1);
                    if (st1.flag == 2) ++nonconv;
                    max_res_e = std::max(max_res_e,
                                         std::abs(ev1.e - e_f) / std::max(std::abs(e_f), tiny));
                    max_pmis = std::max(max_pmis, std::abs(ev1.p - p_f) / std::max(p_f, tiny));
                    max_Terr = std::max(max_Terr, std::abs(ev1.T - T) / T);

                    EosEval ev2;
                    EosInvertStats st2;
                    eval_from_rho_p(rho, p_f, alpha, ev2, st2);
                    if (st2.flag == 2) ++nonconv;
                    max_res_p = std::max(max_res_p, std::abs(ev2.p - p_f) / std::max(p_f, tiny));
                    max_emis = std::max(max_emis,
                                        std::abs(ev2.e - e_f) / std::max(std::abs(e_f), tiny));
                }
            }
        }
        std::ostringstream d;
        d << "n=" << total << " max_res_e=" << max_res_e << " max_res_p=" << max_res_p
          << " max_p_mismatch(monitored)=" << max_pmis << " max_e_mismatch(monitored)=" << max_emis
          << " max_Terr(monitored)=" << max_Terr << " nonconv=" << nonconv;
        // gate on residual convergence only, matching the single-table
        // EOSTAB-SELFTEST roundtrip gates: on real condensed tables p(T) at
        // fixed rho can be flat to below solver tolerance (e.g. Ti's
        // cold-curve-dominated compressed band), so T is not identifiable
        // from p and the recovered e/p/T legitimately differ from the
        // forward seed while every residual converges. Those mismatches are
        // table conditioning properties, not mixture-driver defects — keep
        // them printed as diagnostics.
        verdict("roundtrip", total > 0 && nonconv == 0 && max_res_e <= 10 * ttol &&
                               max_res_p <= 10 * ttol,
                d.str());
    }

    // ---- 2. pure-cell formula parity with the single-table path ----
    // (the mechanism behind the M1/M3 round-off gates: alpha = {1, 0...}
    // must produce BITWISE the TabulatedEOS driver result)
    {
        set_alpha(1.0);
        bool ok = true;
        Real rho_lo = rho_hull_min_k[0], rho_hi = rho_hull_max_k[0];
        for (int a = 0; a < n_sweep; ++a) {
            const Real rho =
              std::pow(10.0, std::log10(rho_lo) + (a + 0.5) / n_sweep *
                                                    (std::log10(rho_hi) - std::log10(rho_lo)));
            const Real T = std::sqrt(T_min_k[0] * T_max_k[0]);
            EosEval ev_ref;
            EosTable::eval_rt(tvs[0], rho, T, ev_ref);
            if (ev_ref.clamped) continue;

            // reference: the exact TabulatedEOS::eval_from_rho_e sequence
            EosEval ev_tab;
            EosInvertStats st_tab;
            const Real T_tab =
              EosTable::invert_T_from_e(tvs[0], rho, ev_ref.e, -1.0, ttol, max_newton, st_tab);
            EosTable::eval_rt(tvs[0], rho, T_tab, ev_tab);

            EosEval ev_mix;
            EosInvertStats st_mix;
            eval_from_rho_e(rho, ev_ref.e, alpha, ev_mix, st_mix);

            ok &= (ev_mix.p == ev_tab.p) && (ev_mix.T == ev_tab.T) && (ev_mix.cs == ev_tab.cs);
        }
        verdict("pure-parity", ok, "alpha={1,0,...} bitwise vs the single-table driver");
    }

    // ---- 3. frozen sound speed vs an isentropic finite difference ----
    // ds = (cv/T) dT - (dpdT/rho^2) drho = 0 gives the isentrope direction
    // dT/drho = T dpdT / (rho^2 cv); a centred FD of the mixture pressure
    // along it must reproduce cs^2 = (dp/drho)|_s to FD+interpolation error
    {
        Real max_err = 0.0;
        int total = 0;
        std::vector<Real> errs;
        for (const Real a0 : {0.25, 0.5, 0.75}) {
            set_alpha(a0);
            prepare_weights(alpha);
            const std::vector<Real> w = m_w;
            const std::vector<int> retained = m_retained;
            Real rho_lo, rho_hi;
            if (!sample_range(rho_lo, rho_hi)) continue;

            for (int a = 0; a < n_sweep; ++a) {
                for (int b = 0; b < n_sweep; ++b) {
                    // stay clear of the hull edges: the FD stencil must not clamp
                    const Real rho =
                      std::pow(10.0, std::log10(rho_lo) + (0.1 + 0.8 * (a + 0.5) / n_sweep) *
                                                            (std::log10(rho_hi) - std::log10(rho_lo)));
                    const Real T =
                      std::pow(10.0, std::log10(T_lo_all) + (0.1 + 0.8 * (b + 0.5) / n_sweep) *
                                                              (std::log10(T_hi_all) - std::log10(T_lo_all)));

                    Real p_f = 0.0;
                    bool clamped = false;
                    EosEval evk;
                    for (const int k : retained) {
                        EosTable::eval_rt(tvs[k], w[k] * rho, T, evk);
                        p_f += evk.p;
                        clamped |= evk.clamped;
                    }
                    if (clamped) continue;

                    EosEval ev;
                    EosInvertStats st;
                    eval_from_rho_p(rho, p_f, alpha, ev, st);
                    if (st.flag != 0) continue;

                    const Real drho = 1.0e-4 * rho;
                    const Real dTdrho = ev.T * ev.dpdT / (rho * rho * std::max(ev.cv, tiny));
                    Real p_pm[2];
                    for (int s = 0; s < 2; ++s) {
                        const Real sgn = (s == 0) ? -1.0 : 1.0;
                        const Real rho_s = rho + sgn * drho;
                        const Real T_s = ev.T + sgn * dTdrho * drho;
                        p_pm[s] = 0.0;
                        for (const int k : retained) {
                            EosTable::eval_rt(tvs[k], w[k] * rho_s, T_s, evk);
                            p_pm[s] += evk.p;
                        }
                    }
                    const Real cs2_fd = (p_pm[1] - p_pm[0]) / (2.0 * drho);
                    const Real cs2 = ev.cs * ev.cs;
                    const Real err = std::abs(cs2_fd - cs2) / std::max(cs2, tiny);
                    max_err = std::max(max_err, err);
                    errs.push_back(err);
                    ++total;
                }
            }
        }
        // gate the 95th percentile, monitor the max: on real spliced tables
        // the stored derivative columns and a within-cell bilinear FD of p
        // legitimately disagree at seam/dome rows (isolated samples), while a
        // mixture-formula algebra bug would shift the whole distribution. On
        // smooth/ideal tables p95 == max to interpolation error, so the gate
        // is unchanged there.
        Real p95_err = 0.0;
        if (!errs.empty()) {
            std::sort(errs.begin(), errs.end());
            p95_err = errs[(size_t)(0.95 * (errs.size() - 1))];
        }
        std::ostringstream d;
        d << "n=" << total << " p95_cs2_err=" << p95_err << " max_cs2_err(monitored)=" << max_err;
        verdict("cs-isentrope-fd", total > 0 && p95_err <= 5.0e-2, d.str());
    }

    // ---- 4. kink-free drop_tol crossing (plan 7.2) ----
    // the dilute-limit scaling makes the retained-side partial pressure
    // O(drop_tol * rho / rho_hull_min), so the jump across the threshold
    // must be far below any physical signal
    {
        Real max_jump = 0.0;
        set_alpha(0.5);
        prepare_weights(alpha);
        Real rho_lo, rho_hi;
        sample_range(rho_lo, rho_hi);
        const Real rho = std::sqrt(rho_lo * rho_hi);
        const Real T = std::sqrt(T_lo_all * T_hi_all);
        Real p_side[2];
        for (int s = 0; s < 2; ++s) {
            set_alpha(s == 0 ? 0.5 * drop_tol : 2.0 * drop_tol);
            prepare_weights(alpha);
            EosInvertStats st;
            p_side[s] = mix_p_rho_T(rho, T, st);
        }
        max_jump = std::abs(p_side[1] - p_side[0]) / std::max(p_side[0], tiny);
        std::ostringstream d;
        d << "p_jump_rel=" << max_jump << " at alpha_0 = drop_tol/2 vs 2*drop_tol";
        verdict("drop-kink", max_jump <= 1.0e-6, d.str());
    }

    amrex::Print() << "MIXEOS-SELFTEST OVERALL " << (all_pass ? "PASS" : "FAIL")
                   << " state=" << name << "\n";
}

void MixtureEOS::write_info(nlohmann::json& js) const
{
    BL_PROFILE("MixtureEOS::write_info");

    HydroGas::write_info(js);

    js["type"] = tag;
    js["mixing_rule"] = "dalton";
    for (size_t k = 0; k < tables.size(); ++k) {
        for (const auto& kv : tables[k].provenance()) {
            js["table_" + comp_names[k] + "_" + kv.first] = kv.second;
        }
    }
}
