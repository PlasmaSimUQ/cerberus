#ifndef MFP_RIEMANN_STATE_PROBE_H
#define MFP_RIEMANN_STATE_PROBE_H

// ===========================================================================
// TEMPORARY DIAGNOSTIC -- delete once the l>=2 flight failure is identified.
// See rho0_roomtemp/FAILURE_HUNT_PLAN.md. Enable with -DMFP_RIEMANN_PROBE via
// Make.local (NOT on the make command line -- a command-line CXXFLAGS+= or
// USERSuffix= REPLACES the makefile's value, which both de-optimises the build
// and suppresses the .PARTICLES suffix, so the binary is neither comparable
// nor correctly named).
//
// Both observed FPE sites are OVERFLOWS on inputs that are already bad:
//   MFP_hlle_general_eos.cpp  fvR[Eden] = uR*(nrgR + pR)
//   MFP_hlle_general_eos.cpp  the HLL average -- its (sR - sL) divide is
//                             already guarded and that counter has never fired
// so the stack address says nothing about provenance. This reports the FACE
// and the CELL instead.
// ===========================================================================

#ifdef MFP_RIEMANN_PROBE

    #include "MFP_utility.H"

    #include <AMReX_Array4.H>
    #include <AMReX_Geometry.H>
    #include <AMReX_ParallelDescriptor.H>
    #include <AMReX_Print.H>
    #include <AMReX_REAL.H>

    #include <cmath>
    #include <set>
    #include <string>
    #include <tuple>

namespace mfp_riemann_probe
{

// PER-SITE budgets. A single shared counter is wrong: the calc_fluxes probe and
// the two solver-entry probes would each burn reports on the same face, and the
// ~180 floor-repaired faces that precede the fatal one would exhaust a shared
// budget before the interesting event.
constexpr int max_reports_per_site = 60;

// Corruption tiers. Operating point: rho_amb 1.1805e-4, rho_casing 2.465,
// p ~1.2e-7, u_slug 4.92 (100 km/s) or 0.0492 (1 km/s).
constexpr amrex::Real u_absurd = 1.0e3;    // |u| this large is m/rho at a floored rho
constexpr amrex::Real pe_absurd = 1.0e12;  // p or e far past any representable shock
constexpr amrex::Real rho_absurd = 1.0e3;  // ~400x solid casing -> corruption
                                           // (WITHOUT this tier the trap can fire
                                           //  via a huge rho while u, p, e all
                                           //  look sane and NOTHING is reported)

inline int& counter(int site)
{
    static int n[4] = {0, 0, 0, 0};
    return n[site];
}

// Dedupe: with CTU the flux box is grown transversally, so the same face is
// visited from more than one box and "first event" is not otherwise unique.
inline std::set<std::tuple<int, int, int, int, int>>& seen()
{
    static std::set<std::tuple<int, int, int, int, int>> s;
    return s;
}

// Classify on PRIMITIVES ONLY. e and a are not available at the calc_fluxes
// site (they are computed inside solve() by get_face_eval_from_prim), and
// calling the EOS here could trap inside the probe and would perturb the
// inversion tallies.
inline std::string classify(amrex::Real rho, amrex::Real u, amrex::Real v,
                            amrex::Real w, amrex::Real p)
{
    if (!std::isfinite(rho) || !std::isfinite(u) || !std::isfinite(v) ||
        !std::isfinite(w) || !std::isfinite(p))
        return "NON-FINITE (nan/inf) -> arrived corrupt";
    if (std::abs(rho) > rho_absurd) return "HUGE density -> corruption upstream";
    if (std::abs(u) > u_absurd || std::abs(v) > u_absurd || std::abs(w) > u_absurd)
        return "HUGE |velocity| -> m/rho at a floored density";
    if (std::abs(p) > pe_absurd) return "HUGE pressure -> EOS inversion runaway";
    if (rho <= 0.0) return "NON-POSITIVE density";
    if (p < 0.0) return "NEGATIVE pressure";
    return "";
}

// Would the HLLE energy flux overflow? uR*(nrg + p) with nrg = rho*e + KE.
// This tier is UNCAPPED: with traps on it can fire at most once, so it cannot
// flood, and it guarantees the fatal face is reported even if the budget above
// is already spent.
inline bool trap_imminent(amrex::Real rho, amrex::Real u, amrex::Real v,
                          amrex::Real w, amrex::Real p)
{
    const amrex::Real ke = 0.5 * std::abs(rho) * (u * u + v * v + w * w);
    const amrex::Real mag = std::abs(u) * (ke + std::abs(p));
    return !std::isfinite(mag) || mag > 1.0e300;
}

inline void dump(const char* tag, const amrex::Vector<amrex::Real>& S, int nprim)
{
    amrex::AllPrint() << "      " << tag << " prim:";
    for (int n = 0; n < nprim; ++n) amrex::AllPrint() << " [" << n << "]=" << S[n];
    amrex::AllPrint() << "\n";
}

// The full report. Called from EulerianState::calc_fluxes with the RAW gathered
// face states -- BEFORE apply_prim_floor and BEFORE transform_global2local, so
// the rho/p pathology is still visible and u/v/w are still in global axes.
template <typename A4>
void probe_face(int site, const char* where, int pass, int level, int grid, int dir,
                int i, int j, int k, const amrex::Real* dx,
                const amrex::Vector<amrex::Real>& L, const amrex::Vector<amrex::Real>& R,
                int nprim, const A4& cons, int ncons, const int* index)
{
    const amrex::Real rL = L[+HydroDef::PrimIdx::Density], uL = L[+HydroDef::PrimIdx::Xvel];
    const amrex::Real vL = L[+HydroDef::PrimIdx::Yvel], wL = L[+HydroDef::PrimIdx::Zvel];
    const amrex::Real pL = L[+HydroDef::PrimIdx::Prs];
    const amrex::Real rR = R[+HydroDef::PrimIdx::Density], uR = R[+HydroDef::PrimIdx::Xvel];
    const amrex::Real vR = R[+HydroDef::PrimIdx::Yvel], wR = R[+HydroDef::PrimIdx::Zvel];
    const amrex::Real pR = R[+HydroDef::PrimIdx::Prs];

    std::string why = classify(rL, uL, vL, wL, pL);
    const char* side = "L";
    if (why.empty()) {
        why = classify(rR, uR, vR, wR, pR);
        side = "R";
    }
    const bool imminent = trap_imminent(rL, uL, vL, wL, pL) || trap_imminent(rR, uR, vR, wR, pR);
    if (why.empty() && !imminent) return;
    if (imminent && why.empty()) why = "TRAP IMMINENT: u*(KE+p) would overflow";

    if (!imminent) {
        if (counter(site) >= max_reports_per_site) return;
        if (!seen().insert(std::make_tuple(level, dir, i, j, k)).second) return;
        ++counter(site);
    }

    const amrex::Real* plo = amrex::DefaultGeometry().ProbLo();
    // cell centres of the two cells sharing this face (the face itself sits at
    // plo[dir] + i*dx[dir] in direction dir)
    const amrex::Real xR = plo[0] + (i + 0.5) * dx[0];
    const amrex::Real yR = plo[1] + (j + 0.5) * dx[1];
    const amrex::Real zR = plo[2] + (k + 0.5) * dx[2];

    amrex::AllPrint()
      << "\n[RIEMANN-PROBE] rank " << amrex::ParallelDescriptor::MyProc() << "  " << where
      << "  CTU-pass=" << pass << (imminent ? "  *** TRAP IMMINENT ***" : "") << "\n"
      << "    " << why << "   (bad side: " << side << ")\n"
      << "    level=" << level << " grid=" << grid << " dir=" << dir << "  face(i,j,k)=(" << i
      << "," << j << "," << k << ")\n"
      << "    R-cell centre  x=" << xR << " y=" << yR << " z=" << zR << "   dx=(" << dx[0] << ","
      << dx[1] << "," << dx[2] << ")\n"
      << "    L-cell index   (" << i - index[0] << "," << j - index[1] << "," << k - index[2]
      << ")\n";
    dump("L", L, nprim);
    dump("R", R, nprim);
    amrex::AllPrint() << "      L cons:";
    for (int n = 0; n < ncons; ++n)
        amrex::AllPrint() << " [" << n << "]=" << cons(i - index[0], j - index[1], k - index[2], n);
    amrex::AllPrint() << "\n      R cons:";
    for (int n = 0; n < ncons; ++n) amrex::AllPrint() << " [" << n << "]=" << cons(i, j, k, n);
    amrex::AllPrint() << "\n";
}

}  // namespace mfp_riemann_probe

    #define MFP_PROBE_FACE_AT_CELL(site, where, pass, lev, grid, dir, i, j, k, dx, L, R, nprim, \
                                   cons, ncons, index)                                          \
        mfp_riemann_probe::probe_face(site, where, pass, lev, grid, dir, i, j, k, dx, L, R,      \
                                      nprim, cons, ncons, index)

#else
    #define MFP_PROBE_FACE_AT_CELL(site, where, pass, lev, grid, dir, i, j, k, dx, L, R, nprim, \
                                   cons, ncons, index)                                          \
        ((void)0)
#endif  // MFP_RIEMANN_PROBE
#endif  // MFP_RIEMANN_STATE_PROBE_H
