#include "MFP_ausmdv.H"
#include "MFP_utility.H"
#include "MFP_global.H"

using GD = GlobalData;
//================================================================================

std::string HydroAUSMDV::tag = "AUSMDV";
bool HydroAUSMDV::registered = GetRiemannSolverFactory().Register(HydroAUSMDV::tag, RiemannSolverBuilder<HydroAUSMDV>);


HydroAUSMDV::HydroAUSMDV(){}
HydroAUSMDV::HydroAUSMDV(const int i)
{
    idx = i;
    istate = GD::get_state_ptr(i);
}

void HydroAUSMDV::solve(Vector<Real> &L,
                        Vector<Real> &R,
                        Vector<Real> &F,
                        Real* shk)
{
    BL_PROFILE("HydroAUSMDV::solve");

    // get the data out of the passed in arrays
    Real apL = L[+HydroState::FluxIdx::Alpha];
    Real gamL = L[+HydroState::FluxIdx::Gamma];
    Real rL = L[+HydroState::FluxIdx::Density];
    Real uL = L[+HydroState::FluxIdx::Xvel];
    Real vL = L[+HydroState::FluxIdx::Yvel];
    Real wL = L[+HydroState::FluxIdx::Zvel];
    Real pL = L[+HydroState::FluxIdx::Prs];
    Real trL = apL*rL;

    // get the data out of the passed in arrays
    Real apR = R[+HydroState::FluxIdx::Alpha];
    Real gamR = R[+HydroState::FluxIdx::Gamma];
    Real rR = R[+HydroState::FluxIdx::Density];
    Real uR = R[+HydroState::FluxIdx::Xvel];
    Real vR = R[+HydroState::FluxIdx::Yvel];
    Real wR = R[+HydroState::FluxIdx::Zvel];
    Real pR = R[+HydroState::FluxIdx::Prs];
    Real trR = apR*rR;

    // constants
    const Real K_SWITCH = 10.0;
    const Real C_EFIX = 0.125;

    // derived quantities

    Real nrgL = pL/(gamL - 1) + 0.5*rL*(uL*uL + vL*vL + wL*wL);
    Real pLrL = pL / rL;
    Real aL = sqrt(gamL*pL/rL);
    Real HL = (nrgL + pL)/rL;

    Real nrgR = pR/(gamR-1) + 0.5*rR*(uR*uR + vR*vR + wR*wR);
    Real aR = sqrt(gamR*pR/rR);
    Real pRrR = pR / rR;
    Real HR = (nrgR + pR)/rR;


    //
    // This is the main part of the flux calculator.
    // Weighting parameters (eqn 32) for velocity splitting.
    Real alphaL = 2.0 * pLrL / (pLrL + pRrR);
    Real alphaR = 2.0 * pRrR / (pLrL + pRrR);
    // Common sound speed (eqn 33) and Mach numbers.
    Real am = std::max(aL, aR);
    Real ML = uL / am;
    Real MR = uR / am;
    // Left state:
    // pressure splitting (eqn 34)
    // and velocity splitting (eqn 30)
    Real duL = 0.5 * (uL + std::abs(uL));
    Real pLplus, uLplus;
    if (std::abs(ML) <= 1.0) {
        pLplus = pL * (ML + 1.0) * (ML + 1.0) * (2.0 - ML) * 0.25;
        uLplus = alphaL * ((uL + am) * (uL + am) / (4.0 * am) - duL) + duL;
    } else {
        pLplus = pL * duL / uL;
        uLplus = duL;
    }
    // Right state:
    // pressure splitting (eqn 34)
    // and velocity splitting (eqn 31)
    Real duR = 0.5 * (uR - std::abs(uR));
    Real pRminus, uRminus;
    if (std::abs(MR) <= 1.0) {
        pRminus = pR * (MR - 1.0) * (MR - 1.0) * (2.0 + MR) * 0.25;
        uRminus = alphaR * (-(uR - am) * (uR - am) / (4.0 * am) - duR) + duR;
    } else {
        pRminus = pR * duR / uR;
        uRminus = duR;
    }
    // Mass Flux (eqn 29)
    // The mass flux is relative to the moving interface.
    Real ru_half = uLplus * rL + uRminus * rR;
    Real tr_half = uLplus * trL + uRminus * trR;
    // Pressure flux (eqn 34)
    Real p_half = pLplus + pRminus;
    // Momentum flux: normal direction
    // Compute blending parameter s (eqn 37),
    // the momentum flux for AUSMV (eqn 21) and AUSMD (eqn 21)
    // and blend (eqn 36).
    Real dp = pL - pR;
    dp = K_SWITCH * std::abs(dp) / std::min(pL, pR);
    Real  s = 0.5 * std::min(1.0, dp);
    Real ru2_AUSMV = uLplus * rL * uL + uRminus * rR * uR;
    Real ru2_AUSMD = 0.5 * (ru_half * (uL + uR) - std::abs(ru_half) * (uR - uL));
    Real ru2_half = (0.5 + s) * ru2_AUSMV + (0.5 - s) * ru2_AUSMD;
    //
    // Assemble components of the flux vector.
    F[+HydroState::ConsIdx::Density] = ru_half;
    F[+HydroState::ConsIdx::Tracer] = tr_half;
    F[+HydroState::ConsIdx::Xmom] = ru2_half+p_half;
    if (ru_half >= 0.0) {
        // Wind is blowing from the left.
        F[+HydroState::ConsIdx::Ymom] = ru_half*vL;
        F[+HydroState::ConsIdx::Zmom] = ru_half*wL;
        F[+HydroState::ConsIdx::Eden] = ru_half*HL;
    } else {
        // Wind is blowing from the right.
        F[+HydroState::ConsIdx::Ymom] = ru_half*vR;
        F[+HydroState::ConsIdx::Zmom] = ru_half*wR;
        F[+HydroState::ConsIdx::Eden] = ru_half*HR;
    }
    //
    // Apply entropy fix (section 3.5 in Wada and Liou's paper)
    bool caseA = ((uL - aL) < 0.0) && ((uR - aR) > 0.0);
    bool caseB = ((uL + aL) < 0.0) && ((uR + aR) > 0.0);
    //
    Real d_ua = 0.0;
    if (caseA && !caseB) {
        d_ua = C_EFIX * ((uR - aR) - (uL - aL));
    }
    if (caseB && !caseA) {
        d_ua = C_EFIX * ((uR + aR) - (uL + aL));
    }
    //
    if (d_ua != 0.0) {
        F[+HydroState::ConsIdx::Density] -= d_ua*(rR - rL);
        F[+HydroState::ConsIdx::Tracer]  -= d_ua*(trR - trL);
        F[+HydroState::ConsIdx::Xmom]    -= d_ua*(rR*uR - rL*uL);
        F[+HydroState::ConsIdx::Ymom]    -= d_ua*(rR*vR - rL*vL);
        F[+HydroState::ConsIdx::Zmom]    -= d_ua*(rR*wR - rL*wL);
        F[+HydroState::ConsIdx::Eden]    -= d_ua*(rR*HR - rL*HL);
    } // end of entropy fix (d_ua != 0)
}

bool HydroAUSMDV::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isHydro) {
        return false;
    }
    return true;
}
