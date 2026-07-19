#include "MFP_hlle_general_eos.H"

#include "MFP_hydro_gas.H"
#include "MFP_utility.H"

//================================================================================
std::string HydroHLLEGeneralEOS::tag = "HLLE_general_eos";
bool HydroHLLEGeneralEOS::registered = GetHydroRiemannSolverFactory().Register(
  HydroHLLEGeneralEOS::tag, HydroRiemannSolverBuilder<HydroHLLEGeneralEOS>);

HydroHLLEGeneralEOS::HydroHLLEGeneralEOS() {}
HydroHLLEGeneralEOS::HydroHLLEGeneralEOS(const sol::table& def)
{
    BL_PROFILE("HydroHLLEGeneralEOS::HydroHLLEGeneralEOS");

    const int n_cons = def["n_cons"];
    const int n_tracer = def["n_tracer"];

    fvL.resize(n_cons);
    svL.resize(n_cons);
    fvR.resize(n_cons);
    svR.resize(n_cons);
    trL.resize(n_tracer);
    trR.resize(n_tracer);

    // n_cons() is ConsIdx::NUM + n_tracers, i.e. it ALREADY counts the tracers;
    // 'n_cons + n_tracer' would double-count them and overrun the caller's F
    n_flux = +HydroDef::ConsIdx::NUM + n_tracer;
}

void HydroHLLEGeneralEOS::solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk)
{
    BL_PROFILE("HydroHLLEGeneralEOS::solve");

    AMREX_ASSERT(gas != nullptr);  // wired by HydroState::set_flux

    // swaps 1+2: specific internal energy and sound speed from ONE combined
    // gas-model face evaluation per side instead of p/(gam-1) and sqrt(gam*p/rho)
    Real eL, aL;
    gas->get_face_eval_from_prim(L, eL, aL);
    Real eR, aR;
    gas->get_face_eval_from_prim(R, eR, aR);

    solve(L, R, F, shk, eL, aL, eR, aR);
}

void HydroHLLEGeneralEOS::solve(Vector<Real>& L,
                                Vector<Real>& R,
                                Vector<Real>& F,
                                Real* shk,
                                Real eL,
                                Real aL,
                                Real eR,
                                Real aR)
{
    BL_PROFILE("HydroHLLEGeneralEOS::solve_reuse");

    const int n_alpha = L.size() - +HydroDef::PrimIdx::NUM;

    // get the data out of the passed in arrays
    Real rhoL = L[+HydroDef::PrimIdx::Density];
    Real uL = L[+HydroDef::PrimIdx::Xvel];
    Real vL = L[+HydroDef::PrimIdx::Yvel];
    Real wL = L[+HydroDef::PrimIdx::Zvel];
    Real pL = L[+HydroDef::PrimIdx::Prs];
    Real nrgL = rhoL * eL + 0.5 * rhoL * (uL * uL + vL * vL + wL * wL);

    for (int i = +HydroDef::PrimIdx::NUM; i < L.size(); ++i) {
        trL[i - +HydroDef::PrimIdx::NUM] = L[i] * rhoL;
    }

    // get the data out of the passed in arrays
    Real rhoR = R[+HydroDef::PrimIdx::Density];
    Real uR = R[+HydroDef::PrimIdx::Xvel];
    Real vR = R[+HydroDef::PrimIdx::Yvel];
    Real wR = R[+HydroDef::PrimIdx::Zvel];
    Real pR = R[+HydroDef::PrimIdx::Prs];
    Real nrgR = rhoR * eR + 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);

    for (int i = +HydroDef::PrimIdx::NUM; i < R.size(); ++i) {
        trR[i - +HydroDef::PrimIdx::NUM] = R[i] * rhoR;
    }

    // speeds
    Real sL = std::min(uL - aL, uR - aR);
    Real sR = std::max(uL + aL, uR + aR);

    // SPECULATIVE — this guard has never fired. Treat it as unproven code.
    //
    // It was added while chasing the case_1_eos failure on the hypothesis that
    // the HLL average below was dividing by a collapsing (sR - sL). That
    // hypothesis was REFUTED by measurement: the counter read zero across every
    // run, and the real cause turned out to be the p_star guard in
    // MFP_hllc_general_eos.cpp (since removed). It is kept only because the
    // divide is genuinely unguarded otherwise and the test is a few flops.
    //
    // The reasoning for the hazard still stands on paper: the spread collapses
    // when both sound speeds vanish and the normal velocities coincide — the
    // population HLLC's guard 2 diverts here — so being star-less removes the
    // p/(rho*(S_L-u)) blowup but not this one. It has simply never been
    // observed. If you are auditing this file, note that nothing here is
    // validated by a passing test; deleting it should be uncontroversial if it
    // is still cold after more cases have run.
    //
    // Widening (rather than substituting a different flux) is the conservative
    // direction: a wider fan is a more diffusive HLL flux and preserves the
    // scheme. The 1e-12 factor is a guess, not a measured threshold.
    const Real v_scale = std::max(std::max(std::abs(uL), std::abs(uR)), std::max(aL, aR));
    const Real ds_min = 1.0e-12 * std::max(v_scale, 1.0e-12);
    if ((sR - sL) < ds_min) {
        const Real mid = 0.5 * (sL + sR);
        sL = mid - 0.5 * ds_min;
        sR = mid + 0.5 * ds_min;
#ifdef MFP_SOLVER_DIAG
        ++n_fan_widened;
#endif
    }

    if (sL >= 0.0) {
        // flux vector L
        F[+HydroDef::ConsIdx::Density] = rhoL * uL;
        F[+HydroDef::ConsIdx::Xmom] = rhoL * uL * uL + pL;
        F[+HydroDef::ConsIdx::Ymom] = rhoL * uL * vL;
        F[+HydroDef::ConsIdx::Zmom] = rhoL * uL * wL;
        F[+HydroDef::ConsIdx::Eden] = uL * (nrgL + pL);
        for (int i = 0; i < n_alpha; ++i) { F[+HydroDef::ConsIdx::NUM + i] = uL * trL[i]; }
    } else if ((sL <= 0.0) && (sR >= 0.0)) {
        // flux vector L
        fvL[+HydroDef::ConsIdx::Density] = rhoL * uL;
        fvL[+HydroDef::ConsIdx::Xmom] = rhoL * uL * uL + pL;
        fvL[+HydroDef::ConsIdx::Ymom] = rhoL * uL * vL;
        fvL[+HydroDef::ConsIdx::Zmom] = rhoL * uL * wL;
        fvL[+HydroDef::ConsIdx::Eden] = uL * (nrgL + pL);
        for (int i = 0; i < n_alpha; ++i) { fvL[+HydroDef::ConsIdx::NUM + i] = uL * trL[i]; }

        // flux vector R
        fvR[+HydroDef::ConsIdx::Density] = rhoR * uR;
        fvR[+HydroDef::ConsIdx::Xmom] = rhoR * uR * uR + pR;
        fvR[+HydroDef::ConsIdx::Ymom] = rhoR * uR * vR;
        fvR[+HydroDef::ConsIdx::Zmom] = rhoR * uR * wR;
        fvR[+HydroDef::ConsIdx::Eden] = uR * (nrgR + pR);
        for (int i = 0; i < n_alpha; ++i) { fvR[+HydroDef::ConsIdx::NUM + i] = uR * trR[i]; }

        // state vector L
        svL[+HydroDef::ConsIdx::Density] = rhoL;
        svL[+HydroDef::ConsIdx::Xmom] = rhoL * uL;
        svL[+HydroDef::ConsIdx::Ymom] = rhoL * vL;
        svL[+HydroDef::ConsIdx::Zmom] = rhoL * wL;
        svL[+HydroDef::ConsIdx::Eden] = nrgL;
        for (int i = 0; i < n_alpha; ++i) { svL[+HydroDef::ConsIdx::NUM + i] = trL[i]; }

        // state vector R
        svR[+HydroDef::ConsIdx::Density] = rhoR;
        svR[+HydroDef::ConsIdx::Xmom] = rhoR * uR;
        svR[+HydroDef::ConsIdx::Ymom] = rhoR * vR;
        svR[+HydroDef::ConsIdx::Zmom] = rhoR * wR;
        svR[+HydroDef::ConsIdx::Eden] = nrgR;
        for (int i = 0; i < n_alpha; ++i) { svR[+HydroDef::ConsIdx::NUM + i] = trR[i]; }

        for (int i = 0; i < +HydroDef::ConsIdx::NUM + n_alpha; ++i) {
            F[i] = (sR * fvL[i] - sL * fvR[i] + sL * sR * (svR[i] - svL[i])) / (sR - sL);
        }
    } else {
        F[+HydroDef::ConsIdx::Density] = rhoR * uR;
        F[+HydroDef::ConsIdx::Xmom] = rhoR * uR * uR + pR;
        F[+HydroDef::ConsIdx::Ymom] = rhoR * uR * vR;
        F[+HydroDef::ConsIdx::Zmom] = rhoR * uR * wR;
        F[+HydroDef::ConsIdx::Eden] = uR * (nrgR + pR);
        for (int i = 0; i < n_alpha; ++i) { F[+HydroDef::ConsIdx::NUM + i] = uR * trR[i]; }
    }
}
