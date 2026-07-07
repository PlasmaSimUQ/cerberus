#ifndef MFP_EOS_TABLE_H
#define MFP_EOS_TABLE_H

// Tabulated-EOS core (plan W4, Stage 2 — doc/eos_implementation_plan.md).
//
// Reads the canonical `.eostab` format (spec frozen 2026-07-07, mirrored in
// Exec/python_analysis/eos_table_prep.py and Exec/testing/EOS-Table/README.md)
// and provides bilinear evaluation on the (log10 rho, log10 T) surface plus
// guarded inversions (given e or p, recover T; given p and T, recover rho).
//
// Deliberately has NO dependency on the HydroGas interface — the gas-model
// wrapper arrives in Stage 3 (W6). CPU-only (plan D6): storage is plain
// std::vector; evaluation goes through the non-owning POD view below, which
// is the calling convention a future GPU port needs (init/evaluate split).

#include "sol.hpp"

#include <AMReX_REAL.H>

#include <map>
#include <string>
#include <vector>

// Non-owning, passed-by-value view of a loaded table (POD: plain-old-data,
// safe to copy into any evaluation context). Axes are reconstructed from
// (min, spacing, count) — never stored as arrays, matching the file format.
struct EosTableView {
    int n_rho = 0, n_T = 0;
    amrex::Real lrho_min = 0, dlrho = 0;  // log10-uniform density axis
    amrex::Real lT_min = 0, dlT = 0;      // log10-uniform temperature axis

    const amrex::Real* p = nullptr;       // pressure
    const amrex::Real* e = nullptr;       // specific internal energy
    const amrex::Real* dpdT = nullptr;    // smoothed derivative tables --
    const amrex::Real* dpdrho = nullptr;  //  outputs only, never used as
    const amrex::Real* cv = nullptr;      //  the Newton slope (plan D5)
    const amrex::Real* dedrho = nullptr;

    // pre-inverted seed maps (Athena++-informed, plan D5/D8); null if the
    // file omits them (they are optional blocks)
    const amrex::Real* T_of_e = nullptr;  // T(rho, e) on (lrho x le)
    const amrex::Real* T_of_p = nullptr;  // T(rho, p) on (lrho x lp)
    amrex::Real le_min = 0, dle = 0;
    int n_e = 0;
    amrex::Real lp_min = 0, dlp = 0;
    int n_p = 0;

    // REQUIRED block: 1.0 = inside the source-data hull, 0.0 = cell filled
    // by the conditioning tool (the FPEOS source is ragged)
    const amrex::Real* hull = nullptr;

    int idx(int i, int j) const { return i * n_T + j; }  // i = rho, j = T

    // true if every corner of cell (i, j) carries source data
    bool cell_in_hull(int i, int j) const
    {
        return hull[idx(i, j)] > 0.5 && hull[idx(i + 1, j)] > 0.5 &&
               hull[idx(i, j + 1)] > 0.5 && hull[idx(i + 1, j + 1)] > 0.5;
    }
};

// Outputs of one (rho, T) evaluation
struct EosEval {
    amrex::Real rho = 0, T = 0, p = 0, e = 0;
    amrex::Real cs = 0, gam1 = 0, dpde = 0, dpdr_e = 0;  // Riemann contract
    amrex::Real cv = 0, dpdT = 0, dpdrho = 0, dedrho = 0;
    bool clamped = false;  // query was outside the grid or in a filled cell
};

// Bookkeeping for one inversion call
struct EosInvertStats {
    int iters = 0;       // Newton/bisection iterations taken
    int bisections = 0;  // how many of those were bisection fallbacks
    bool seeded = false; // start point came from an inverse map
    // 0 = converged; 1 = target outside the table's range at this line
    // (clamped to the nearer end); 2 = max_newton exhausted (best value
    // returned) — never aborts (plan D5)
    int flag = 0;
};

class EosTable
{
  public:
    EosTable() = default;

    // Parse + hard-validate a .eostab file (dimensional, as stored).
    // Malformed input -> amrex::Abort with context. Every rank reads (D9).
    void load(const std::string& path);

    // Divide once by the Cerberus reference quantities (plan D7). Must be
    // called exactly once, after load. The self-test never calls this
    // (dimensional mode) — Stage 3 wires MFP::*_ref in here.
    void nondimensionalise(amrex::Real rho_ref,
                           amrex::Real T_ref,
                           amrex::Real prs_ref,
                           amrex::Real u_ref);

    EosTableView view() const;

    // --- evaluation / inversion (static: operate on the view only) ---

    // bilinear interpolation in (log10 rho, log10 T), hard hull clamping,
    // then the Riemann-contract identities (non-convexity floored: cs^2>=0)
    static void eval_rt(const EosTableView& v,
                        amrex::Real rho,
                        amrex::Real T,
                        EosEval& out);

    // Given rho and a target e (or p), recover T. Inverse-map seed + Newton
    // whose slope is the analytic derivative of the bilinear interpolant
    // itself; hull-bracketed with bisection fallback; convergence is
    // declared on the e (resp. p) residual, NOT the T increment.
    static amrex::Real invert_T_from_e(const EosTableView& v,
                                       amrex::Real rho,
                                       amrex::Real e_target,
                                       amrex::Real T_guess,
                                       amrex::Real ttol,
                                       int max_newton,
                                       EosInvertStats& stats);
    static amrex::Real invert_T_from_p(const EosTableView& v,
                                       amrex::Real rho,
                                       amrex::Real p_target,
                                       amrex::Real T_guess,
                                       amrex::Real ttol,
                                       int max_newton,
                                       EosInvertStats& stats);

    // Given T and a target p, recover rho (tp mode; initial-condition fills)
    static amrex::Real invert_rho_from_p(const EosTableView& v,
                                         amrex::Real T,
                                         amrex::Real p_target,
                                         amrex::Real rho_guess,
                                         amrex::Real ttol,
                                         int max_newton,
                                         EosInvertStats& stats);

    // --- hull accessors (in current units: dimensional until nondim) ---
    amrex::Real rho_min() const;
    amrex::Real rho_max() const;
    amrex::Real T_min() const;
    amrex::Real T_max() const;

    const std::map<std::string, std::string>& provenance() const { return prov; }
    bool is_nondimensional() const { return nondim_done; }

    // Registers the debug-only Lua hook `eos_table_self_test(path, n_sweep)`
    // (SDF self-test pattern; compiled to a no-op in non-DEBUG builds)
    static void register_with_lua(sol::state& lua);

#ifdef AMREX_DEBUG
    // One-zone self-test harness (Stage-2 W5): prints one
    // "EOSTAB-SELFTEST[...] PASS/FAIL ..." line per check for check.py.
    static bool self_test(const std::string& path, int n_sweep);
#endif

    // inversion controls (Lua-configurable via the Stage-3 gas model)
    amrex::Real ttol = 1.0e-10;
    int max_newton = 100;

  private:
    int n_rho = 0, n_T = 0, n_e = 0, n_p = 0;
    amrex::Real lrho_min = 0, lrho_max = 0, lT_min = 0, lT_max = 0;
    amrex::Real le_min = 0, le_max = 0, lp_min = 0, lp_max = 0;

    std::vector<amrex::Real> p_v, e_v, dpdT_v, dpdrho_v, cv_v, dedrho_v;
    std::vector<amrex::Real> T_of_e_v, T_of_p_v, hull_v;

    std::map<std::string, std::string> prov;
    bool loaded = false;
    bool nondim_done = false;
};

#endif  // MFP_EOS_TABLE_H
