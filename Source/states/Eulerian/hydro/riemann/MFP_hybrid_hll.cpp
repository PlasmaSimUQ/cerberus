#include "MFP_hybrid_hll.H"

#include "MFP.H"
#include "MFP_state.H"
#include "MFP_utility.H"

//================================================================================

std::string HydroHybridHLL::tag = "HLLE/HLLC";
bool HydroHybridHLL::registered = GetHydroRiemannSolverFactory().Register(
  HydroHybridHLL::tag, HydroRiemannSolverBuilder<HydroHybridHLL>);

HydroHybridHLL::HydroHybridHLL() {}
HydroHybridHLL::HydroHybridHLL(const sol::table& def)
{
    BL_PROFILE("HydroHybridHLL::HydroHybridHLL");

    hllc = HydroHLLC(def);
    hlle = HydroHLLE(def);

    const int n_cons = def["n_cons"];
    const int n_tracer = def["n_tracer"];

    // n_cons() is ConsIdx::NUM + n_tracers, i.e. it ALREADY counts the tracers;
    // 'n_cons + n_tracer' would double-count them. Latent while the blend loop
    // stopped at ConsIdx::NUM, live once it blends the full n_flux vector.
    n_flux = +HydroDef::ConsIdx::NUM + n_tracer;

    F_hlle.resize(n_flux);
    F_hllc.resize(n_flux);
}

void HydroHybridHLL::solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk)
{
    BL_PROFILE("HydroHybridHLL::solve");

    constexpr Real eps = 1e-14;

    if (*shk < eps) {
        hllc.solve(L, R, F, shk);
    } else if ((1.0 - *shk) < eps) {
        hlle.solve(L, R, F, shk);
    } else {
        hllc.solve(L, R, F_hllc, shk);
        hlle.solve(L, R, F_hlle, shk);
        // blend the FULL flux vector (conserved hydro + tracers); stopping at
        // ConsIdx::NUM leaves the tracer slots of F holding stale values from
        // the previous face, corrupting passive-scalar transport in every
        // shock-transition cell
        for (int i = 0; i < n_flux; ++i) {
            F[i] = (1.0 - *shk) * F_hllc[i] + *shk * F_hlle[i];
        }
    }
}
