#include "MFP_ausmdv.H"
#include "MFP_utility.H"
#include "MFP.H"

//================================================================================

std::string HydroAUSMDV::tag = "AUSMDV";
bool HydroAUSMDV::registered = GetHydroRiemannSolverFactory().Register(HydroAUSMDV::tag, HydroRiemannSolverBuilder<HydroAUSMDV>);


HydroAUSMDV::HydroAUSMDV(){}
HydroAUSMDV::HydroAUSMDV(const sol::table &def)
{
    const int n_tracer = def["n_tracer"];

    trL.resize(n_tracer);
    trR.resize(n_tracer);

    n_flux = +HydroDef::ConsIdx::NUM + n_tracer;
}

void HydroAUSMDV::solve(Vector<Real> &L,
                                Vector<Real> &R,
                                Vector<Real> &F,
                                Real* shk)
{
    BL_PROFILE("HydroAUSMDV::solve");

    const int n_alpha = L.size() - +HydroDef::PrimIdx::NUM;

    // get the data out of the passed in arrays
    Real gamL = L[+HydroDef::PrimIdx::Gamma];
    Real rL = L[+HydroDef::PrimIdx::Density];
    Real uL = L[+HydroDef::PrimIdx::Xvel];
    Real vL = L[+HydroDef::PrimIdx::Yvel];
    Real wL = L[+HydroDef::PrimIdx::Zvel];
    Real pL = L[+HydroDef::PrimIdx::Prs];

    for (int i=+HydroDef::PrimIdx::NUM; i<L.size(); ++ i) {
        trL[i-+HydroDef::PrimIdx::NUM]= L[i]*rL;
    }

    // get the data out of the passed in arrays
    Real gamR = R[+HydroDef::PrimIdx::Gamma];
    Real rR = R[+HydroDef::PrimIdx::Density];
    Real uR = R[+HydroDef::PrimIdx::Xvel];
    Real vR = R[+HydroDef::PrimIdx::Yvel];
    Real wR = R[+HydroDef::PrimIdx::Zvel];
    Real pR = R[+HydroDef::PrimIdx::Prs];


    for (int i=+HydroDef::PrimIdx::NUM; i<R.size(); ++ i) {
        trR[i-+HydroDef::PrimIdx::NUM]= R[i]*rR;
    }

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

    for (int i=0; i<n_alpha; ++i) {
          F[+HydroDef::ConsIdx::NUM+i]  = uLplus * trL[i] + uRminus * trR[i];
    }

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
    F[+HydroDef::ConsIdx::Density] = ru_half;
    F[+HydroDef::ConsIdx::Xmom] = ru2_half+p_half;
    if (ru_half >= 0.0) {
        // Wind is blowing from the left.
        F[+HydroDef::ConsIdx::Ymom] = ru_half*vL;
        F[+HydroDef::ConsIdx::Zmom] = ru_half*wL;
        F[+HydroDef::ConsIdx::Eden] = ru_half*HL;
    } else {
        // Wind is blowing from the right.
        F[+HydroDef::ConsIdx::Ymom] = ru_half*vR;
        F[+HydroDef::ConsIdx::Zmom] = ru_half*wR;
        F[+HydroDef::ConsIdx::Eden] = ru_half*HR;
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
        F[+HydroDef::ConsIdx::Density] -= d_ua*(rR - rL);
        F[+HydroDef::ConsIdx::Xmom]    -= d_ua*(rR*uR - rL*uL);
        F[+HydroDef::ConsIdx::Ymom]    -= d_ua*(rR*vR - rL*vL);
        F[+HydroDef::ConsIdx::Zmom]    -= d_ua*(rR*wR - rL*wL);
        F[+HydroDef::ConsIdx::Eden]    -= d_ua*(rR*HR - rL*HL);

        for (int i=0; i<n_alpha; ++i) {
            F[+HydroDef::ConsIdx::NUM+i]  -= d_ua*(trR[i] - trL[i]);
        }
    } // end of entropy fix (d_ua != 0)
}
