#include "MFP_eos_table.H"

#include <AMReX.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

using amrex::Real;

// ===========================================================================
// helpers
// ===========================================================================

namespace
{

// locate x on a log10-uniform axis: cell index i in [0, n-2], fraction f in
// [0, 1]; clamped reports whether x fell outside the axis range
void locate(Real x, Real xmin, Real dx, int n, int& i, Real& f, bool& clamped)
{
    const Real u = (x - xmin) / dx;
    if (u < 0.0 || u > (Real)(n - 1)) clamped = true;
    i = (int)std::floor(u);
    i = std::max(0, std::min(n - 2, i));
    f = u - (Real)i;
    f = std::max((Real)0.0, std::min((Real)1.0, f));
}

Real bilin(const Real* a, const EosTableView& v, int i, int j, Real fx, Real fy)
{
    const Real a00 = a[v.idx(i, j)], a01 = a[v.idx(i, j + 1)];
    const Real a10 = a[v.idx(i + 1, j)], a11 = a[v.idx(i + 1, j + 1)];
    return (1 - fx) * ((1 - fy) * a00 + fy * a01) + fx * ((1 - fy) * a10 + fy * a11);
}

std::string trim(const std::string& s)
{
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// "key=value" tokens on a header line -> value for the requested key
std::string kv_token(const std::string& line, const std::string& key)
{
    std::istringstream iss(line);
    std::string tok;
    while (iss >> tok) {
        if (tok.rfind(key + "=", 0) == 0) return tok.substr(key.size() + 1);
    }
    return "";
}

// shared 1-D root find along one axis of the bilinear surface.
// value(t) and dvalue/dt come from the interpolant itself (plan D5); the
// bracket [lo, hi] is the axis range; any Newton step leaving the bracket
// becomes a bisection step; convergence is on the residual in the target
// quantity. `fval` evaluates the interpolated block at axis coordinate t and
// also returns the in-cell slope d(value)/dt.
template <typename F>
Real invert_1d(F fval,
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
        // target not attainable on this line -> clamp to the nearer end
        stats.flag = 1;
        return (std::abs(flo) < std::abs(fhi)) ? lo : hi;
    }

    Real t = std::max(lo, std::min(hi, t_seed));
    for (stats.iters = 1; stats.iters <= max_newton; ++stats.iters) {
        Real slope;
        const Real ft = fval(t, slope) - target;
        if (std::abs(ft) <= ttol * scale) return t;

        // maintain the bracket
        if ((ft > 0) == (flo > 0)) {
            lo = t;
            flo = ft;
        } else {
            hi = t;
            fhi = ft;
        }
        if (hi - lo < 1.0e-14) return 0.5 * (lo + hi);

        // Newton candidate; bisection if the slope is unusable or the step
        // leaves the bracket. The |ft| <= |slope|*(hi-lo) guard bounds the
        // step by the bracket width BEFORE dividing, so ft/slope can never
        // overflow (matters with amrex.fpe_trap_* enabled: near filled-cell
        // boundaries the interpolant slope underflows to ~0).
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
    stats.flag = 2;  // not converged; return the best bracket midpoint
    return 0.5 * (lo + hi);
}

}  // namespace

// ===========================================================================
// evaluation
// ===========================================================================

void EosTable::eval_rt(const EosTableView& v, Real rho, Real T, EosEval& out)
{
    const Real tiny = std::numeric_limits<Real>::min();
    const Real lr = std::log10(std::max(rho, tiny));
    const Real lt = std::log10(std::max(T, tiny));

    bool clamped = false;
    int i, j;
    Real fx, fy;
    locate(lr, v.lrho_min, v.dlrho, v.n_rho, i, fx, clamped);
    locate(lt, v.lT_min, v.dlT, v.n_T, j, fy, clamped);

    out.rho = rho;
    out.T = T;
    out.p = bilin(v.p, v, i, j, fx, fy);
    out.e = bilin(v.e, v, i, j, fx, fy);
    out.dpdT = bilin(v.dpdT, v, i, j, fx, fy);
    out.dpdrho = bilin(v.dpdrho, v, i, j, fx, fy);
    out.cv = bilin(v.cv, v, i, j, fx, fy);
    out.dedrho = bilin(v.dedrho, v, i, j, fx, fy);

    // filled (out-of-source-hull) cells count as clamped territory
    if (!v.cell_in_hull(i, j)) clamped = true;
    out.clamped = clamped;

    // Riemann-contract identities with the non-convexity floor (plan D5)
    out.dpde = out.dpdT / out.cv;
    out.dpdr_e = out.dpdrho - out.dpde * out.dedrho;
    const Real cs2 =
        std::max(out.dpdrho + (T / (rho * rho)) * out.dpdT * out.dpdT / out.cv,
                 (Real)0.0);
    out.cs = std::sqrt(cs2);
    out.gam1 = rho * cs2 / std::max(out.p, tiny);
}

// ===========================================================================
// inversions
// ===========================================================================

Real EosTable::invert_T_from_e(const EosTableView& v,
                               Real rho,
                               Real e_target,
                               Real T_guess,
                               Real ttol,
                               int max_newton,
                               EosInvertStats& stats)
{
    const Real tiny = std::numeric_limits<Real>::min();
    const Real lr = std::log10(std::max(rho, tiny));
    bool cl = false;
    int i;
    Real fx;
    locate(lr, v.lrho_min, v.dlrho, v.n_rho, i, fx, cl);

    // value + in-cell slope of the interpolated e along log10 T at fixed rho
    auto fval = [&](Real lt, Real& slope) {
        bool c2 = false;
        int j;
        Real fy;
        locate(lt, v.lT_min, v.dlT, v.n_T, j, fy, c2);
        const Real ea = (1 - fx) * v.e[v.idx(i, j)] + fx * v.e[v.idx(i + 1, j)];
        const Real eb =
            (1 - fx) * v.e[v.idx(i, j + 1)] + fx * v.e[v.idx(i + 1, j + 1)];
        slope = (eb - ea) / v.dlT;
        return (1 - fy) * ea + fy * eb;
    };

    // seed: inverse map lookup (one bilinear read), else caller guess, else
    // the axis midpoint (plan D8)
    const Real lo = v.lT_min, hi = v.lT_min + v.dlT * (v.n_T - 1);
    Real lt0 = 0.5 * (lo + hi);
    if (v.T_of_e != nullptr) {
        bool c2 = false;
        int je;
        Real fe;
        const Real le = std::log10(std::max(e_target, tiny));
        locate(le, v.le_min, v.dle, v.n_e, je, fe, c2);
        const Real Ts = (1 - fx) * ((1 - fe) * v.T_of_e[i * v.n_e + je] +
                                    fe * v.T_of_e[i * v.n_e + je + 1]) +
                        fx * ((1 - fe) * v.T_of_e[(i + 1) * v.n_e + je] +
                              fe * v.T_of_e[(i + 1) * v.n_e + je + 1]);
        if (Ts > 0.0) {
            lt0 = std::log10(Ts);
            stats.seeded = true;
        }
    } else if (T_guess > 0.0) {
        lt0 = std::log10(T_guess);
    }

    const Real lt = invert_1d(fval, lo, hi, e_target, lt0, ttol, max_newton, stats);
    return std::pow((Real)10.0, lt);
}

Real EosTable::invert_T_from_p(const EosTableView& v,
                               Real rho,
                               Real p_target,
                               Real T_guess,
                               Real ttol,
                               int max_newton,
                               EosInvertStats& stats)
{
    const Real tiny = std::numeric_limits<Real>::min();
    const Real lr = std::log10(std::max(rho, tiny));
    bool cl = false;
    int i;
    Real fx;
    locate(lr, v.lrho_min, v.dlrho, v.n_rho, i, fx, cl);

    auto fval = [&](Real lt, Real& slope) {
        bool c2 = false;
        int j;
        Real fy;
        locate(lt, v.lT_min, v.dlT, v.n_T, j, fy, c2);
        const Real pa = (1 - fx) * v.p[v.idx(i, j)] + fx * v.p[v.idx(i + 1, j)];
        const Real pb =
            (1 - fx) * v.p[v.idx(i, j + 1)] + fx * v.p[v.idx(i + 1, j + 1)];
        slope = (pb - pa) / v.dlT;
        return (1 - fy) * pa + fy * pb;
    };

    const Real lo = v.lT_min, hi = v.lT_min + v.dlT * (v.n_T - 1);
    Real lt0 = 0.5 * (lo + hi);
    if (v.T_of_p != nullptr) {
        bool c2 = false;
        int jp;
        Real fp;
        const Real lp = std::log10(std::max(p_target, tiny));
        locate(lp, v.lp_min, v.dlp, v.n_p, jp, fp, c2);
        const Real Ts = (1 - fx) * ((1 - fp) * v.T_of_p[i * v.n_p + jp] +
                                    fp * v.T_of_p[i * v.n_p + jp + 1]) +
                        fx * ((1 - fp) * v.T_of_p[(i + 1) * v.n_p + jp] +
                              fp * v.T_of_p[(i + 1) * v.n_p + jp + 1]);
        if (Ts > 0.0) {
            lt0 = std::log10(Ts);
            stats.seeded = true;
        }
    } else if (T_guess > 0.0) {
        lt0 = std::log10(T_guess);
    }

    const Real lt = invert_1d(fval, lo, hi, p_target, lt0, ttol, max_newton, stats);
    return std::pow((Real)10.0, lt);
}

Real EosTable::invert_rho_from_p(const EosTableView& v,
                                 Real T,
                                 Real p_target,
                                 Real rho_guess,
                                 Real ttol,
                                 int max_newton,
                                 EosInvertStats& stats)
{
    const Real tiny = std::numeric_limits<Real>::min();
    const Real lt = std::log10(std::max(T, tiny));
    bool cl = false;
    int j;
    Real fy;
    locate(lt, v.lT_min, v.dlT, v.n_T, j, fy, cl);

    auto fval = [&](Real lr, Real& slope) {
        bool c2 = false;
        int i;
        Real fx;
        locate(lr, v.lrho_min, v.dlrho, v.n_rho, i, fx, c2);
        const Real pa = (1 - fy) * v.p[v.idx(i, j)] + fy * v.p[v.idx(i, j + 1)];
        const Real pb =
            (1 - fy) * v.p[v.idx(i + 1, j)] + fy * v.p[v.idx(i + 1, j + 1)];
        slope = (pb - pa) / v.dlrho;
        return (1 - fx) * pa + fx * pb;
    };

    const Real lo = v.lrho_min, hi = v.lrho_min + v.dlrho * (v.n_rho - 1);
    Real lr0 = (rho_guess > 0.0) ? std::log10(rho_guess) : 0.5 * (lo + hi);

    const Real lr = invert_1d(fval, lo, hi, p_target, lr0, ttol, max_newton, stats);
    return std::pow((Real)10.0, lr);
}

// ===========================================================================
// loading / units
// ===========================================================================

void EosTable::load(const std::string& path)
{
    std::ifstream f(path);
    if (!f) amrex::Abort("EosTable::load: cannot open '" + path + "'");

    std::string line;
    std::getline(f, line);
    if (trim(line).rfind("EOSTAB 1", 0) != 0)
        amrex::Abort("EosTable::load: bad magic in '" + path +
                     "' (expected 'EOSTAB 1', got '" + trim(line) + "')");

    // ---- header ----
    std::string block_name;
    while (std::getline(f, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        line = trim(line);
        if (line.empty()) continue;
        if (line.rfind("block:", 0) == 0) {
            block_name = trim(line.substr(6));
            break;
        }
        const auto colon = line.find(':');
        if (colon == std::string::npos)
            amrex::Abort("EosTable::load: malformed header line '" + line + "'");
        const std::string key = trim(line.substr(0, colon));
        const std::string val = trim(line.substr(colon + 1));
        if (key == "grid") {
            n_rho = std::stoi(kv_token(val, "n_rho"));
            n_T = std::stoi(kv_token(val, "n_T"));
        } else if (key == "lrho") {
            std::istringstream iss(val);
            iss >> lrho_min >> lrho_max;
        } else if (key == "lT") {
            std::istringstream iss(val);
            iss >> lT_min >> lT_max;
        } else if (key == "le") {
            std::istringstream iss(val);
            iss >> le_min >> le_max;
            n_e = std::stoi(kv_token(val, "n_e"));
        } else if (key == "lp") {
            std::istringstream iss(val);
            iss >> lp_min >> lp_max;
            n_p = std::stoi(kv_token(val, "n_p"));
        } else {
            prov[key] = val;
        }
    }
    if (n_rho < 2 || n_T < 2)
        amrex::Abort("EosTable::load: grid dims must be > 1 in '" + path + "'");
    if (!(lrho_max > lrho_min) || !(lT_max > lT_min))
        amrex::Abort("EosTable::load: bad axis ranges in '" + path + "'");

    // ---- blocks ----
    const size_t fwd = (size_t)n_rho * n_T;
    while (!block_name.empty()) {
        std::vector<Real>* dst = nullptr;
        size_t count = fwd;
        if (block_name == "p") dst = &p_v;
        else if (block_name == "e") dst = &e_v;
        else if (block_name == "dpdT") dst = &dpdT_v;
        else if (block_name == "dpdrho") dst = &dpdrho_v;
        else if (block_name == "cv") dst = &cv_v;
        else if (block_name == "dedrho") dst = &dedrho_v;
        else if (block_name == "hull") dst = &hull_v;
        else if (block_name == "T_of_e") {
            dst = &T_of_e_v;
            count = (size_t)n_rho * n_e;
        } else if (block_name == "T_of_p") {
            dst = &T_of_p_v;
            count = (size_t)n_rho * n_p;
        } else {
            amrex::Abort("EosTable::load: unknown block '" + block_name +
                         "' in '" + path + "'");
        }
        if (count == 0)
            amrex::Abort("EosTable::load: block '" + block_name +
                         "' has no axis declared in '" + path + "'");
        dst->resize(count);
        for (size_t k = 0; k < count; ++k) {
            if (!(f >> (*dst)[k]))
                amrex::Abort("EosTable::load: block '" + block_name +
                             "' too short in '" + path + "'");
        }
        // advance to the next block header (or EOF)
        block_name.clear();
        while (std::getline(f, line)) {
            line = trim(line);
            if (line.empty()) continue;
            if (line.rfind("block:", 0) == 0) {
                block_name = trim(line.substr(6));
                break;
            }
            amrex::Abort("EosTable::load: unexpected content '" + line +
                         "' in '" + path + "'");
        }
    }

    // ---- validation ----
    auto require = [&](const std::vector<Real>& vv, const char* nm) {
        if (vv.size() != fwd)
            amrex::Abort("EosTable::load: required block '" + std::string(nm) +
                         "' missing or short in '" + path + "'");
        for (const Real x : vv) {
            if (!std::isfinite(x))
                amrex::Abort("EosTable::load: non-finite value in block '" +
                             std::string(nm) + "' of '" + path + "'");
        }
    };
    require(p_v, "p");
    require(e_v, "e");
    require(dpdT_v, "dpdT");
    require(dpdrho_v, "dpdrho");
    require(cv_v, "cv");
    require(dedrho_v, "dedrho");
    require(hull_v, "hull");
    for (const Real x : cv_v) {
        if (!(x > 0.0))
            amrex::Abort("EosTable::load: cv <= 0 found in '" + path + "'");
    }
    bool any_hull = false;
    for (const Real x : hull_v) {
        if (x != 0.0 && x != 1.0)
            amrex::Abort("EosTable::load: hull values must be 0 or 1 in '" +
                         path + "'");
        any_hull |= (x == 1.0);
    }
    if (!any_hull) amrex::Abort("EosTable::load: empty hull in '" + path + "'");
    if (!T_of_e_v.empty()) {
        for (const Real x : T_of_e_v) {
            if (!std::isfinite(x))
                amrex::Abort("EosTable::load: non-finite T_of_e in '" + path + "'");
        }
    }
    if (!T_of_p_v.empty()) {
        for (const Real x : T_of_p_v) {
            if (!std::isfinite(x))
                amrex::Abort("EosTable::load: non-finite T_of_p in '" + path + "'");
        }
    }

    loaded = true;

    // load report: dims + the provenance lines that matter for consistency
    amrex::Print() << "EosTable: loaded '" << path << "' (" << n_rho << " x "
                   << n_T << ", inverse maps: "
                   << (T_of_e_v.empty() ? "no" : "yes") << ")\n";
    for (const char* k : {"material", "source", "e_shift", "conditioning"}) {
        const auto it = prov.find(k);
        if (it != prov.end())
            amrex::Print() << "  " << k << ": " << it->second << "\n";
    }
}

void EosTable::nondimensionalise(Real rho_ref, Real T_ref, Real prs_ref, Real u_ref)
{
    if (!loaded) amrex::Abort("EosTable::nondimensionalise before load");
    if (nondim_done) amrex::Abort("EosTable::nondimensionalise called twice");
    const Real u2 = u_ref * u_ref;

    for (Real& x : p_v) x /= prs_ref;
    for (Real& x : e_v) x /= u2;
    for (Real& x : dpdT_v) x *= T_ref / prs_ref;
    for (Real& x : dpdrho_v) x *= rho_ref / prs_ref;
    for (Real& x : cv_v) x *= T_ref / u2;
    for (Real& x : dedrho_v) x *= rho_ref / u2;
    for (Real& x : T_of_e_v) x /= T_ref;
    for (Real& x : T_of_p_v) x /= T_ref;

    const Real dlr = std::log10(rho_ref), dlt = std::log10(T_ref);
    lrho_min -= dlr;
    lrho_max -= dlr;
    lT_min -= dlt;
    lT_max -= dlt;
    le_min -= std::log10(u2);
    le_max -= std::log10(u2);
    lp_min -= std::log10(prs_ref);
    lp_max -= std::log10(prs_ref);

    nondim_done = true;
}

EosTableView EosTable::view() const
{
    if (!loaded) amrex::Abort("EosTable::view before load");
    EosTableView v;
    v.n_rho = n_rho;
    v.n_T = n_T;
    v.lrho_min = lrho_min;
    v.dlrho = (lrho_max - lrho_min) / (n_rho - 1);
    v.lT_min = lT_min;
    v.dlT = (lT_max - lT_min) / (n_T - 1);
    v.p = p_v.data();
    v.e = e_v.data();
    v.dpdT = dpdT_v.data();
    v.dpdrho = dpdrho_v.data();
    v.cv = cv_v.data();
    v.dedrho = dedrho_v.data();
    v.hull = hull_v.data();
    if (!T_of_e_v.empty()) {
        v.T_of_e = T_of_e_v.data();
        v.le_min = le_min;
        v.dle = (le_max - le_min) / (n_e - 1);
        v.n_e = n_e;
    }
    if (!T_of_p_v.empty()) {
        v.T_of_p = T_of_p_v.data();
        v.lp_min = lp_min;
        v.dlp = (lp_max - lp_min) / (n_p - 1);
        v.n_p = n_p;
    }
    return v;
}

Real EosTable::rho_min() const { return std::pow(10.0, lrho_min); }
Real EosTable::rho_max() const { return std::pow(10.0, lrho_max); }
Real EosTable::T_min() const { return std::pow(10.0, lT_min); }
Real EosTable::T_max() const { return std::pow(10.0, lT_max); }

// ===========================================================================
// Lua registration + debug self-test (SDF hook pattern:
// MFP_ebgeometry_nodeshared.cpp — free function, #ifdef AMREX_DEBUG)
// ===========================================================================

void EosTable::register_with_lua(sol::state& lua)
{
#ifdef AMREX_DEBUG
    lua.set_function("eos_table_self_test", &EosTable::self_test);
#else
    amrex::ignore_unused(lua);
#endif
}

#ifdef AMREX_DEBUG

namespace
{
// simple accumulator for iteration statistics
struct SweepStats {
    std::vector<int> iters;
    int bisections = 0, nonconv = 0, seeded = 0, total = 0;
    Real max_res = 0.0, max_Terr = 0.0;
    void add(const EosInvertStats& s, Real res, Real Terr)
    {
        ++total;
        iters.push_back(s.iters);
        bisections += s.bisections;
        if (s.flag == 2) ++nonconv;
        if (s.seeded) ++seeded;
        max_res = std::max(max_res, res);
        max_Terr = std::max(max_Terr, Terr);
    }
    int pct(Real q) const
    {
        if (iters.empty()) return 0;
        auto v = iters;
        std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, (size_t)(q * v.size()))];
    }
    int imax() const { return iters.empty() ? 0 : *std::max_element(iters.begin(), iters.end()); }
};
}  // namespace

bool EosTable::self_test(const std::string& path, int n_sweep)
{
    // Runs in DIMENSIONAL mode (no nondimensionalise): the hook fires during
    // Lua config execution, before reference quantities are guaranteed set
    // (review amendment; Stage-3 exercises the nondimensional path instead).
    bool all_pass = true;
    auto verdict = [&](const std::string& name, bool ok, const std::string& detail) {
        all_pass &= ok;
        amrex::Print() << "EOSTAB-SELFTEST[" << name << "] "
                       << (ok ? "PASS " : "FAIL ") << detail << "\n";
    };

    EosTable tab;
    tab.load(path);
    const EosTableView v = tab.view();
    const Real ttol = tab.ttol;
    const int max_newton = tab.max_newton;

    // ---- 1. reader integrity ----
    {
        int in_hull = 0;
        for (const Real x : tab.hull_v) in_hull += (x > 0.5);
        std::ostringstream d;
        d << "dims=" << v.n_rho << "x" << v.n_T << " hull_frac="
          << (Real)in_hull / (v.n_rho * v.n_T)
          << " inv_maps=" << (v.T_of_e ? "yes" : "no");
        auto mm = [&](const std::vector<Real>& b, const char* nm) {
            const auto lohi = std::minmax_element(b.begin(), b.end());
            d << " " << nm << "=[" << *lohi.first << "," << *lohi.second << "]";
        };
        mm(tab.p_v, "p");
        mm(tab.e_v, "e");
        mm(tab.cv_v, "cv");
        verdict("reader", true, d.str());
    }

    // sample coordinates: cell-interior points covering the grid rectangle
    auto sample = [&](int k, int n, Real lo, Real d_ax, int n_ax) {
        const Real u = (k + 0.5) / n;  // in (0, 1)
        return lo + u * d_ax * (n_ax - 1);
    };

    // ---- 2. round trips on hull cells (rt->re->rt and rt->rp->rt) ----
    {
        SweepStats se, sp;
        for (int a = 0; a < n_sweep; ++a) {
            for (int b = 0; b < n_sweep; ++b) {
                const Real lr = sample(a, n_sweep, v.lrho_min, v.dlrho, v.n_rho);
                const Real lt = sample(b, n_sweep, v.lT_min, v.dlT, v.n_T);
                const Real rho = std::pow(10.0, lr), T = std::pow(10.0, lt);
                EosEval ev;
                eval_rt(v, rho, T, ev);
                if (ev.clamped) continue;  // filled cells -> check 4 territory

                EosInvertStats st1;
                const Real T1 = invert_T_from_e(v, rho, ev.e, -1.0, ttol, max_newton, st1);
                EosEval e1;
                eval_rt(v, rho, T1, e1);
                se.add(st1, std::abs(e1.e - ev.e) / std::max(std::abs(ev.e), 1e-300),
                       std::abs(T1 - T) / T);

                EosInvertStats st2;
                const Real T2 = invert_T_from_p(v, rho, ev.p, -1.0, ttol, max_newton, st2);
                EosEval e2;
                eval_rt(v, rho, T2, e2);
                sp.add(st2, std::abs(e2.p - ev.p) / std::max(std::abs(ev.p), 1e-300),
                       std::abs(T2 - T) / T);
            }
        }
        for (auto* s : {&se, &sp}) {
            const bool is_e = (s == &se);
            std::ostringstream d;
            d << "n=" << s->total << " max_res=" << s->max_res
              << " max_Terr=" << s->max_Terr << " iters_max=" << s->imax()
              << " iters_p99=" << s->pct(0.99) << " bisect=" << s->bisections
              << " seedrate=" << (s->total ? (Real)s->seeded / s->total : 0)
              << " nonconv=" << s->nonconv;
            verdict(is_e ? "roundtrip-e" : "roundtrip-p",
                    s->total > 0 && s->max_res <= ttol * 10 && s->nonconv == 0,
                    d.str());
        }
    }

    // ---- 3. derivative identities ----
    {
        // (a) recombination consistency of eval outputs (exact identities)
        // (b) smoothed derivative blocks vs finite differences of the value
        //     surface across each cell (agreement to interpolation order)
        Real max_id = 0.0, max_fd = 0.0;
        for (int a = 0; a < n_sweep; ++a) {
            for (int b = 0; b < n_sweep; ++b) {
                const Real lr = sample(a, n_sweep, v.lrho_min, v.dlrho, v.n_rho);
                const Real lt = sample(b, n_sweep, v.lT_min, v.dlT, v.n_T);
                const Real rho = std::pow(10.0, lr), T = std::pow(10.0, lt);
                EosEval ev;
                eval_rt(v, rho, T, ev);
                if (ev.clamped) continue;

                const Real cs2 = ev.cs * ev.cs;
                const Real cs2_id =
                    ev.dpdrho + (T / (rho * rho)) * ev.dpdT * ev.dpdT / ev.cv;
                max_id = std::max(max_id, std::abs(cs2 - std::max(cs2_id, (Real)0.0)) /
                                              std::max(cs2, (Real)1e-300));
                max_id = std::max(max_id,
                                  std::abs(ev.gam1 - rho * cs2 / ev.p) /
                                      std::max(std::abs(ev.gam1), (Real)1e-300));
                max_id = std::max(max_id,
                                  std::abs(ev.dpde - ev.dpdT / ev.cv) /
                                      std::max(std::abs(ev.dpde), (Real)1e-300));

                // FD across the containing cell in LINEAR T and rho
                bool cdum = false;
                int i, j;
                Real fx, fy;
                locate(lr, v.lrho_min, v.dlrho, v.n_rho, i, fx, cdum);
                locate(lt, v.lT_min, v.dlT, v.n_T, j, fy, cdum);
                const Real T1 = std::pow(10.0, v.lT_min + j * v.dlT);
                const Real T2 = std::pow(10.0, v.lT_min + (j + 1) * v.dlT);
                const Real r1 = std::pow(10.0, v.lrho_min + i * v.dlrho);
                const Real r2 = std::pow(10.0, v.lrho_min + (i + 1) * v.dlrho);
                EosEval c00, c01, c10, c11;
                eval_rt(v, r1, T1, c00);
                eval_rt(v, r1, T2, c01);
                eval_rt(v, r2, T1, c10);
                eval_rt(v, r2, T2, c11);
                const Real cv_fd =
                    0.5 * ((c01.e - c00.e) + (c11.e - c10.e)) / (T2 - T1);
                const Real dpdT_fd =
                    0.5 * ((c01.p - c00.p) + (c11.p - c10.p)) / (T2 - T1);
                const Real dpdr_fd =
                    0.5 * ((c10.p - c00.p) + (c11.p - c01.p)) / (r2 - r1);
                // compare the cell-centred secants against the derivative
                // blocks at the CELL CENTRE (mean of the four corners) —
                // comparing against an arbitrary in-cell sample point would
                // measure half-a-cell of variation, not conditioning error
                const Real cv_c = 0.25 * (c00.cv + c01.cv + c10.cv + c11.cv);
                const Real dpdT_c = 0.25 * (c00.dpdT + c01.dpdT + c10.dpdT + c11.dpdT);
                const Real dpdr_c =
                    0.25 * (c00.dpdrho + c01.dpdrho + c10.dpdrho + c11.dpdrho);
                max_fd = std::max(max_fd, std::abs(cv_fd - cv_c) /
                                              std::max(std::abs(cv_c), (Real)1e-300));
                max_fd = std::max(max_fd, std::abs(dpdT_fd - dpdT_c) /
                                              std::max(std::abs(dpdT_c), (Real)1e-300));
                max_fd = std::max(max_fd, std::abs(dpdr_fd - dpdr_c) /
                                              std::max(std::abs(dpdr_c), (Real)1e-300));
            }
        }
        std::ostringstream d;
        d << "max_identity=" << max_id << " max_fd_vs_block=" << max_fd;
        verdict("identities", max_id <= 1e-12, d.str());
        verdict("fd-vs-blocks", max_fd <= 5e-2, d.str());
    }

    // ---- 4. hull / out-of-range behaviour ----
    {
        bool ok = true;
        std::ostringstream d;
        const Real rmid = std::pow(10.0, 0.5 * (v.lrho_min + v.lrho_min + v.dlrho * (v.n_rho - 1)));
        const Real tmid = std::pow(10.0, 0.5 * (v.lT_min + v.lT_min + v.dlT * (v.n_T - 1)));
        const Real probes[4][2] = {{tab.rho_min() * 1e-3, tmid},
                                   {tab.rho_max() * 1e3, tmid},
                                   {rmid, tab.T_min() * 1e-3},
                                   {rmid, tab.T_max() * 1e3}};
        for (const auto& pr : probes) {
            EosEval ev;
            eval_rt(v, pr[0], pr[1], ev);
            const bool fin = std::isfinite(ev.p) && std::isfinite(ev.e) &&
                             std::isfinite(ev.cs) && std::isfinite(ev.T);
            ok &= ev.clamped && fin;
        }
        // a filled cell, if the table has any
        int filled = -1;
        for (int i = 0; i + 1 < v.n_rho && filled < 0; ++i) {
            for (int j = 0; j + 1 < v.n_T && filled < 0; ++j) {
                if (!v.cell_in_hull(i, j)) filled = v.idx(i, j);
            }
        }
        if (filled >= 0) {
            const int i = filled / v.n_T, j = filled % v.n_T;
            EosEval ev;
            eval_rt(v, std::pow(10.0, v.lrho_min + (i + 0.5) * v.dlrho),
                    std::pow(10.0, v.lT_min + (j + 0.5) * v.dlT), ev);
            ok &= ev.clamped && std::isfinite(ev.p) && std::isfinite(ev.cs);
            d << "filled_cell_probe=yes ";
        } else {
            d << "filled_cell_probe=none ";
        }
        d << "rect_probes=4";
        verdict("hull", ok, d.str());
    }

    // ---- 5. degenerate-corner stress (high rho, low T) ----
    {
        SweepStats s;
        const int n = std::max(8, n_sweep / 2);
        for (int a = 0; a < n; ++a) {
            for (int b = 0; b < n; ++b) {
                // top 15% of the density axis, bottom 15% of the T axis
                const Real lr = v.lrho_min + v.dlrho * (v.n_rho - 1) *
                                                 (0.85 + 0.15 * (a + 0.5) / n);
                const Real lt = v.lT_min + v.dlT * (v.n_T - 1) * (0.15 * (b + 0.5) / n);
                const Real rho = std::pow(10.0, lr), T = std::pow(10.0, lt);
                EosEval ev;
                eval_rt(v, rho, T, ev);
                if (ev.clamped) continue;
                EosInvertStats st;
                const Real T1 = invert_T_from_e(v, rho, ev.e, -1.0, ttol, max_newton, st);
                EosEval e1;
                eval_rt(v, rho, T1, e1);
                s.add(st, std::abs(e1.e - ev.e) / std::max(std::abs(ev.e), 1e-300),
                      std::abs(T1 - T) / T);
            }
        }
        std::ostringstream d;
        d << "n=" << s.total << " max_eres=" << s.max_res
          << " max_Terr(monitored)=" << s.max_Terr << " bisect=" << s.bisections
          << " nonconv=" << s.nonconv;
        // pass on the e-residual; T error is monitored only (ill-posed where
        // cv -> 0, plan amendment #1). A fully-filled corner (n=0) passes.
        verdict("corner", s.total == 0 || (s.max_res <= ttol * 10 && s.nonconv == 0),
                d.str());
    }

    amrex::Print() << "EOSTAB-SELFTEST OVERALL " << (all_pass ? "PASS" : "FAIL")
                   << " table=" << path << "\n";
    return all_pass;
}

#endif  // AMREX_DEBUG
