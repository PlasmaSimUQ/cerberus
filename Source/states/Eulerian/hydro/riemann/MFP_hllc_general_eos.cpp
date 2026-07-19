#include "MFP_hllc_general_eos.H"

#include "MFP_hydro_gas.H"
#include "MFP_utility.H"

#include <algorithm>
#include <cmath>

//================================================================================
std::string HydroHLLCGeneralEOS::tag = "HLLC_general_eos";
bool HydroHLLCGeneralEOS::registered = GetHydroRiemannSolverFactory().Register(
  HydroHLLCGeneralEOS::tag, HydroRiemannSolverBuilder<HydroHLLCGeneralEOS>);

HydroHLLCGeneralEOS::HydroHLLCGeneralEOS() {}
HydroHLLCGeneralEOS::HydroHLLCGeneralEOS(const sol::table& def)
{
    BL_PROFILE("HydroHLLCGeneralEOS::HydroHLLCGeneralEOS");

    const int n_cons = def["n_cons"];
    const int n_tracer = def["n_tracer"];

    svLs.resize(n_cons);
    fvL.resize(n_cons);
    svL.resize(n_cons);
    svRs.resize(n_cons);
    fvR.resize(n_cons);
    svR.resize(n_cons);
    trL.resize(n_tracer);
    trR.resize(n_tracer);

    n_flux = +HydroDef::ConsIdx::NUM + n_tracer;

    // embedded star-less fallback for degenerate faces, and its relative floor
    hlle = HydroHLLEGeneralEOS(def);
    eps_rel = def["fallback_eps"].get_or(1.0e-6);
}

void HydroHLLCGeneralEOS::solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk)
{
    BL_PROFILE("HydroHLLCGeneralEOS::solve");

    AMREX_ASSERT(gas != nullptr);  // wired by HydroState::set_flux (STAGE5.md D-d)

    // swaps 1+2 (STAGE5.md): specific internal energy and sound speed from
    // ONE combined gas-model face evaluation per side instead of p/(gam-1) and
    // sqrt(gam*p/rho) — for a tabulated gas both come from the same table
    // inversion, which is what keeps the G8 wall-time budget
    Real eL, aL;
    gas->get_face_eval_from_prim(L, eL, aL);
    Real eR, aR;
    gas->get_face_eval_from_prim(R, eR, aR);

    solve(L, R, F, shk, eL, aL, eR, aR);
}

void HydroHLLCGeneralEOS::solve(Vector<Real>& L,
                                Vector<Real>& R,
                                Vector<Real>& F,
                                Real* shk,
                                Real eL,
                                Real aL,
                                Real eR,
                                Real aR)
{
    BL_PROFILE("HydroHLLCGeneralEOS::solve_reuse");

    const size_t n_alpha = L.size() - +HydroDef::PrimIdx::NUM;
    const size_t n_flux = +HydroDef::ConsIdx::NUM + n_alpha;

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

    // ---- Degeneracy fallback guards (doc/general_eos_hllc_fallback_plan.md).
    // Each is PREDICTIVE — tested before the operation it protects — because
    // FPE traps (USE_PRIM_FLOOR builds) fire at the bad divide, before any
    // post-hoc isfinite check could run. On a trip, substitute the star-less
    // general-EOS HLLE (reusing the face evaluations above) and return. ----

    // guard 1: EOS evaluation finiteness/sign (mode D — off-hull face)
    if (!std::isfinite(eL) || !std::isfinite(eR) || !std::isfinite(aL) || !std::isfinite(aR) ||
        (aL < 0.0) || (aR < 0.0)) {
#ifdef MFP_SOLVER_DIAG
        ++n_fallback.eval;
#endif
        hlle.solve(L, R, F, shk, eL, aL, eR, aR);
        return;
    }

    // guard 2: relative sound-speed floor (modes A, B — cold hull). v_ref is the
    // local characteristic speed; a face whose sound speed is negligible against
    // it makes pL/(rhoL*aL) in the starred energy (and rho_bar*a_bar in S*)
    // blow up. Also catches the a_bar == 0 dead-cell corner.
    const Real a_bar = 0.5 * (aL + aR);
    const Real v_ref = std::max(std::max(std::abs(uL), std::abs(uR)), std::max(aL, aR));
    if ((std::min(aL, aR) < eps_rel * v_ref) || (a_bar <= 0.0)) {
#ifdef MFP_SOLVER_DIAG
        ++n_fallback.sound;
#endif
        hlle.solve(L, R, F, shk, eL, aL, eR, aR);
        return;
    }

    // Calculate wave speeds S_L, S_star and S_R

    Real rho_bar = 0.5 * (rhoL + rhoR);

    Real p_star = 0.5 * (pL + pR) - 0.5 * (uR - uL) * rho_bar * a_bar;

    // NOTE: there is deliberately NO guard on p_star <= 0 here. A negative PVRS
    // pressure looks like incipient vacuum, but p_star is never a divisor in
    // this solver — it only selects the q-factor branch below, and p_star <= pL
    // (trivially true when negative) yields q = 1, i.e. S_L = uL - aL and
    // S_R = uR + aR, the ordinary acoustic/rarefaction case HLLC is designed
    // for. An earlier revision did guard it and fell back to HLLE, which broke
    // case_1_eos outright: it diverted strong-expansion faces into the one
    // regime where HLLE is *worse* than HLLC, since HLL's single averaged star
    // state cannot represent two separating rarefactions and can produce
    // negative internal energy there. See Section 11 of
    // doc/general_eos_hllc_fallback_plan.md.

    // a_bar > 0 was enforced by guard 2, so this division is safe
    Real S_star = 0.5 * (uL + uR) - 0.5 * (pR - pL) / (rho_bar * a_bar);

    // swap 3: where the two-shock q-factor needs a gamma, use the local
    // acoustic gamma (Gamma_1 = a^2*rho/p) at the FACE state — already paid
    // for by the combined face evaluation, so the q-factor costs no extra
    // table calls (STAGE5.md D-c recorded alternative, adopted after the G8
    // measurement; the PVRS-p* variant was physics-identical within every
    // gate but ~1 extra inversion per face). On a gamma-law gas
    // Gamma_1 == gamma, so this matches MFP_hllc.cpp's q exactly — standard
    // HLLC also takes its q-factor gamma from the face, not the mid-state.
    Real qL;
    if (p_star <= pL) {
        qL = 1.0;
    } else {
        const Real g_loc = aL * aL * rhoL / pL;
        qL = std::sqrt(1.0 + ((g_loc + 1.0) / (2 * g_loc)) * (p_star / pL - 1.0));
    }

    Real S_L = uL - aL * qL;

    Real qR;
    if (p_star <= pR) {
        qR = 1.0;
    } else {
        const Real g_loc = aR * aR * rhoR / pR;
        qR = std::sqrt(1.0 + ((g_loc + 1.0) / (2 * g_loc)) * (p_star / pR - 1.0));
    }

    Real S_R = uR + aR * qR;

    // guard 4: wave ordering + relative gap floor (modes C, F). coeff and the
    // tracer star divide by (S_L - S_star)/(S_R - S_star); a collapsed fan
    // (S_R <= S_L) or an acoustic wave merging with the contact (gap -> 0), or
    // an S_star thrown outside [S_L, S_R] by a marginal estimate, makes those
    // singular. Enforcing the ordering here also guarantees coeff > 0, ruling
    // out negative star density/energy.
    const Real dS = S_R - S_L;
    if ((dS <= 0.0) || ((S_star - S_L) < eps_rel * dS) || ((S_R - S_star) < eps_rel * dS)) {
#ifdef MFP_SOLVER_DIAG
        ++n_fallback.wave;
#endif
        hlle.solve(L, R, F, shk, eL, aL, eR, aR);
        return;
    }

    if (S_L >= 0.0) {
        // flux vector L
        F[+HydroDef::ConsIdx::Density] = rhoL * uL;
        F[+HydroDef::ConsIdx::Xmom] = rhoL * uL * uL + pL;
        F[+HydroDef::ConsIdx::Ymom] = rhoL * uL * vL;
        F[+HydroDef::ConsIdx::Zmom] = rhoL * uL * wL;
        F[+HydroDef::ConsIdx::Eden] = uL * (nrgL + pL);

        for (int i = 0; i < n_alpha; ++i) { F[+HydroDef::ConsIdx::NUM + i] = uL * trL[i]; }

        return;
    } else if ((S_L <= 0.0) && (0.0 <= S_star)) {
        // flux vector L
        fvL[+HydroDef::ConsIdx::Density] = rhoL * uL;
        fvL[+HydroDef::ConsIdx::Xmom] = rhoL * uL * uL + pL;
        fvL[+HydroDef::ConsIdx::Ymom] = rhoL * uL * vL;
        fvL[+HydroDef::ConsIdx::Zmom] = rhoL * uL * wL;
        fvL[+HydroDef::ConsIdx::Eden] = uL * (nrgL + pL);

        for (int i = 0; i < n_alpha; ++i) { fvL[+HydroDef::ConsIdx::NUM + i] = uL * trL[i]; }

        // state vector L
        svL[+HydroDef::ConsIdx::Density] = rhoL;
        svL[+HydroDef::ConsIdx::Xmom] = rhoL * uL;
        svL[+HydroDef::ConsIdx::Ymom] = rhoL * vL;
        svL[+HydroDef::ConsIdx::Zmom] = rhoL * wL;
        svL[+HydroDef::ConsIdx::Eden] = nrgL;

        for (int i = 0; i < n_alpha; ++i) { svL[+HydroDef::ConsIdx::NUM + i] = trL[i]; }

        Real coeff = rhoL * ((S_L - uL) / (S_L - S_star));

        svLs[+HydroDef::ConsIdx::Density] = coeff;
        svLs[+HydroDef::ConsIdx::Xmom] = coeff * S_star;
        svLs[+HydroDef::ConsIdx::Ymom] = coeff * vL;
        svLs[+HydroDef::ConsIdx::Zmom] = coeff * wL;
        svLs[+HydroDef::ConsIdx::Eden] =
          coeff * (nrgL / rhoL + (S_star - uL) * (S_star + pL / (rhoL * (S_L - uL))));

        for (int i = 0; i < n_alpha; ++i) {
            svLs[+HydroDef::ConsIdx::NUM + i] = trL[i] * ((S_L - uL) / (S_L - S_star));
        }

        for (int i = 0; i < +HydroDef::ConsIdx::NUM + n_alpha; ++i) {
            F[i] = fvL[i] + S_L * (svLs[i] - svL[i]);
        }

    } else if ((S_star <= 0.0) && (0.0 <= S_R)) {
        // flux vector R
        fvR[+HydroDef::ConsIdx::Density] = rhoR * uR;
        fvR[+HydroDef::ConsIdx::Xmom] = rhoR * uR * uR + pR;
        fvR[+HydroDef::ConsIdx::Ymom] = rhoR * uR * vR;
        fvR[+HydroDef::ConsIdx::Zmom] = rhoR * uR * wR;
        fvR[+HydroDef::ConsIdx::Eden] = uR * (nrgR + pR);

        for (int i = 0; i < n_alpha; ++i) { fvR[+HydroDef::ConsIdx::NUM + i] = uR * trR[i]; }

        // state vector R
        svR[+HydroDef::ConsIdx::Density] = rhoR;
        svR[+HydroDef::ConsIdx::Xmom] = rhoR * uR;
        svR[+HydroDef::ConsIdx::Ymom] = rhoR * vR;
        svR[+HydroDef::ConsIdx::Zmom] = rhoR * wR;
        svR[+HydroDef::ConsIdx::Eden] = nrgR;

        for (int i = 0; i < n_alpha; ++i) { svR[+HydroDef::ConsIdx::NUM + i] = trR[i]; }

        Real coeff = rhoR * ((S_R - uR) / (S_R - S_star));

        svRs[+HydroDef::ConsIdx::Density] = coeff;
        svRs[+HydroDef::ConsIdx::Xmom] = coeff * S_star;
        svRs[+HydroDef::ConsIdx::Ymom] = coeff * vR;
        svRs[+HydroDef::ConsIdx::Zmom] = coeff * wR;
        svRs[+HydroDef::ConsIdx::Eden] =
          coeff * (nrgR / rhoR + (S_star - uR) * (S_star + pR / (rhoR * (S_R - uR))));

        for (int i = 0; i < n_alpha; ++i) {
            svRs[+HydroDef::ConsIdx::NUM + i] = trR[i] * ((S_R - uR) / (S_R - S_star));
        }

        for (int i = 0; i < +HydroDef::ConsIdx::NUM + n_alpha; ++i) {
            F[i] = fvR[i] + S_R * (svRs[i] - svR[i]);
        }

    } else {
        // flux vector R
        F[+HydroDef::ConsIdx::Density] = rhoR * uR;
        F[+HydroDef::ConsIdx::Xmom] = rhoR * uR * uR + pR;
        F[+HydroDef::ConsIdx::Ymom] = rhoR * uR * vR;
        F[+HydroDef::ConsIdx::Zmom] = rhoR * uR * wR;
        F[+HydroDef::ConsIdx::Eden] = uR * (nrgR + pR);

        for (int i = 0; i < n_alpha; ++i) { F[+HydroDef::ConsIdx::NUM + i] = trR[i] * uR; }
    }
}
