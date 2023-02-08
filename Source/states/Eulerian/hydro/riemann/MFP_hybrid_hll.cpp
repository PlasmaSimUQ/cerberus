#include "MFP_hybrid_hll.H"
#include "MFP_utility.H"
#include "MFP.H"
#include "MFP_state.H"

//================================================================================

std::string HydroHybridHLL::tag = "HLLE/HLLC";
bool HydroHybridHLL::registered = GetHydroRiemannSolverFactory().Register(HydroHybridHLL::tag, HydroRiemannSolverBuilder<HydroHybridHLL>);

HydroHybridHLL::HydroHybridHLL() {}
HydroHybridHLL::HydroHybridHLL(const sol::table &def)
{
    hllc = HydroHLLC(def);
    hlle = HydroHLLE(def);

    const int n_cons = def["n_cons"];
    const int n_tracer = def["n_tracer"];

    F_hlle.resize(n_cons+n_tracer);
    F_hllc.resize(n_cons+n_tracer);

    n_flux = n_cons + n_tracer;
}

void HydroHybridHLL::solve(Vector<Real> &L,
                                Vector<Real> &R,
                                Vector<Real> &F,
                                Real* shk)
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
        for (int i=0; i<+HydroDef::ConsIdx::NUM; ++i) {
            F[i] = (1.0 - *shk)*F_hllc[i] + *shk*F_hlle[i];
        }
    }
}
