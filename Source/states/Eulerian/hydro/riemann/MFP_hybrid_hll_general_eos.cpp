#include "MFP_hybrid_hll_general_eos.H"

#include "MFP.H"
#include "MFP_hydro_gas.H"
#include "MFP_state.H"
#include "MFP_utility.H"

//================================================================================

std::string HydroHybridHLLGeneralEOS::tag = "HLLE/HLLC_general_eos";
bool HydroHybridHLLGeneralEOS::registered = GetHydroRiemannSolverFactory().Register(
  HydroHybridHLLGeneralEOS::tag, HydroRiemannSolverBuilder<HydroHybridHLLGeneralEOS>);

HydroHybridHLLGeneralEOS::HydroHybridHLLGeneralEOS() {}
HydroHybridHLLGeneralEOS::HydroHybridHLLGeneralEOS(const sol::table& def)
{
    BL_PROFILE("HydroHybridHLLGeneralEOS::HydroHybridHLLGeneralEOS");

    hllc = HydroHLLCGeneralEOS(def);
    hlle = HydroHLLEGeneralEOS(def);

    const int n_cons = def["n_cons"];
    const int n_tracer = def["n_tracer"];

    // n_cons() is ConsIdx::NUM + n_tracers, i.e. it ALREADY counts the tracers;
    // 'n_cons + n_tracer' would double-count them, and the blend loop below
    // writes n_flux entries into the caller's F -> heap overrun when tracers>0
    n_flux = +HydroDef::ConsIdx::NUM + n_tracer;

    F_hlle.resize(n_flux);
    F_hllc.resize(n_flux);
}

void HydroHybridHLLGeneralEOS::solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk)
{
    BL_PROFILE("HydroHybridHLLGeneralEOS::solve");

    constexpr Real eps = 1e-14;

    if (*shk < eps) {
        hllc.solve(L, R, F, shk);
    } else if ((1.0 - *shk) < eps) {
        hlle.solve(L, R, F, shk);
    } else {
        // transition band: evaluate the table ONCE per side and hand the
        // (e, a) pairs to both sub-solvers, instead of letting each re-invert
        // the table for the same L and R (which would be 4 inversions/face)
        AMREX_ASSERT(gas != nullptr);
        Real eL, aL;
        gas->get_face_eval_from_prim(L, eL, aL);
        Real eR, aR;
        gas->get_face_eval_from_prim(R, eR, aR);

        hllc.solve(L, R, F_hllc, shk, eL, aL, eR, aR);
        hlle.solve(L, R, F_hlle, shk, eL, aL, eR, aR);
        // blend the FULL flux vector (conserved hydro + tracers) so passive
        // scalars are advected consistently with the hydro state they ride on
        for (int i = 0; i < n_flux; ++i) {
            F[i] = (1.0 - *shk) * F_hllc[i] + *shk * F_hlle[i];
        }
    }
}
