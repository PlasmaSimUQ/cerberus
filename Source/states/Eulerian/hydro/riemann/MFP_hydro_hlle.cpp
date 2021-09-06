#include "MFP_hydro_hlle.H"

#include "MFP_utility.H"
#include "MFP.H"
#include "MFP_state.H"

//================================================================================
std::string HydroHLLE::tag = "HLLE";
bool HydroHLLE::registered = GetHydroRiemannSolverFactory().Register(HydroHLLE::tag, HydroRiemannSolverBuilder<HydroHLLE>);

HydroHLLE::HydroHLLE(){}
HydroHLLE::HydroHLLE(const sol::table &def)
{
    const int n_cons = def["n_cons"];
    const int n_tracer = def["n_tracer"];


    fvL.resize(n_cons);
    svL.resize(n_cons);
    fvR.resize(n_cons);
    svR.resize(n_cons);
    trL.resize(n_tracer);
    trR.resize(n_tracer);

}

void HydroHLLE::solve(Vector<Real> &L,
                                Vector<Real> &R,
                                Vector<Real> &F,
                                Real* shk)
{
    BL_PROFILE("HydroHLLE::solve");

    const int n_alpha = L.size() - +HydroDef::PrimIdx::NUM;

    // get the data out of the passed in arrays
    Real rhoL = L[+HydroDef::PrimIdx::Density];
    Real uL = L[+HydroDef::PrimIdx::Xvel];
    Real vL = L[+HydroDef::PrimIdx::Yvel];
    Real wL = L[+HydroDef::PrimIdx::Zvel];
    Real pL = L[+HydroDef::PrimIdx::Prs];
    Real gamL = L[+HydroDef::PrimIdx::Gamma];
    Real nrgL = pL/(gamL - 1) + 0.5*rhoL*(uL*uL + vL*vL + wL*wL);
    Real aL = sqrt(gamL*pL/rhoL);

    for (int i=+HydroDef::PrimIdx::NUM; i<L.size(); ++ i) {
        trL[i-+HydroDef::PrimIdx::NUM]= L[i]*rhoL;
    }

    // get the data out of the passed in arrays
    Real rhoR = R[+HydroDef::PrimIdx::Density];
    Real uR = R[+HydroDef::PrimIdx::Xvel];
    Real vR = R[+HydroDef::PrimIdx::Yvel];
    Real wR = R[+HydroDef::PrimIdx::Zvel];
    Real pR = R[+HydroDef::PrimIdx::Prs];
    Real gamR = R[+HydroDef::PrimIdx::Gamma];
    Real nrgR = pR/(gamR-1) + 0.5*rhoR*(uR*uR + vR*vR + wR*wR);
    Real aR = sqrt(gamR*pR/rhoR);

    for (int i=+HydroDef::PrimIdx::NUM; i<R.size(); ++ i) {
        trR[i-+HydroDef::PrimIdx::NUM]= R[i]*rhoR;
    }

    // speeds
    Real sL = std::min(uL - aL, uR - aR);
    Real sR = std::max(uL + aL, uR + aR);

    if (sL >= 0.0) {
        // flux vector L
        F[+HydroDef::ConsIdx::Density]  = rhoL*uL;
        F[+HydroDef::ConsIdx::Xmom]   = rhoL*uL*uL + pL;
        F[+HydroDef::ConsIdx::Ymom]   = rhoL*uL*vL;
        F[+HydroDef::ConsIdx::Zmom]   = rhoL*uL*wL;
        F[+HydroDef::ConsIdx::Eden]   = uL*(nrgL + pL);
        for (int i=0; i<n_alpha; ++i) {
            F[+HydroDef::ConsIdx::NUM+i]  = uL * trL[i];
        }
    } else if ((sL <= 0.0) && (sR >= 0.0)) {

        // flux vector L
        fvL[+HydroDef::ConsIdx::Density]  = rhoL*uL;
        fvL[+HydroDef::ConsIdx::Xmom]   = rhoL*uL*uL + pL;
        fvL[+HydroDef::ConsIdx::Ymom]   = rhoL*uL*vL;
        fvL[+HydroDef::ConsIdx::Zmom]   = rhoL*uL*wL;
        fvL[+HydroDef::ConsIdx::Eden]   = uL*(nrgL + pL);
        for (int i=0; i<n_alpha; ++i) {
            fvL[+HydroDef::ConsIdx::NUM+i]  = uL * trL[i];
        }

        // flux vector R
        fvR[+HydroDef::ConsIdx::Density]  = rhoR*uR;
        fvR[+HydroDef::ConsIdx::Xmom]   = rhoR*uR*uR + pR;
        fvR[+HydroDef::ConsIdx::Ymom]   = rhoR*uR*vR;
        fvR[+HydroDef::ConsIdx::Zmom]   = rhoR*uR*wR;
        fvR[+HydroDef::ConsIdx::Eden]   = uR*(nrgR + pR);
        for (int i=0; i<n_alpha; ++i) {
            fvR[+HydroDef::ConsIdx::NUM+i]  = uR * trR[i];
        }

        // state vector L
        svL[+HydroDef::ConsIdx::Density]  = rhoL;
        svL[+HydroDef::ConsIdx::Xmom]   = rhoL*uL;
        svL[+HydroDef::ConsIdx::Ymom]   = rhoL*vL;
        svL[+HydroDef::ConsIdx::Zmom]   = rhoL*wL;
        svL[+HydroDef::ConsIdx::Eden]   = nrgL;
        for (int i=0; i<n_alpha; ++i) {
            svL[+HydroDef::ConsIdx::NUM+i]  = trL[i];
        }

        // state vector R
        svR[+HydroDef::ConsIdx::Density]  = rhoR;
        svR[+HydroDef::ConsIdx::Xmom]   = rhoR*uR;
        svR[+HydroDef::ConsIdx::Ymom]   = rhoR*vR;
        svR[+HydroDef::ConsIdx::Zmom]   = rhoR*wR;
        svR[+HydroDef::ConsIdx::Eden]   = nrgR;
        for (int i=0; i<n_alpha; ++i) {
            svR[+HydroDef::ConsIdx::NUM+i]  = trR[i];
        }

        for (int i=0; i<+HydroDef::ConsIdx::NUM + n_alpha; ++i) {
            F[i] = (sR*fvL[i] - sL*fvR[i] + sL*sR*(svR[i] - svL[i]))/(sR - sL);
        }
    } else {
        F[+HydroDef::ConsIdx::Density]  = rhoR*uR;
        F[+HydroDef::ConsIdx::Xmom]   = rhoR*uR*uR + pR;
        F[+HydroDef::ConsIdx::Ymom]   = rhoR*uR*vR;
        F[+HydroDef::ConsIdx::Zmom]   = rhoR*uR*wR;
        F[+HydroDef::ConsIdx::Eden]   = uR*(nrgR + pR);
        for (int i=0; i<n_alpha; ++i) {
            F[+HydroDef::ConsIdx::NUM+i]  = uR*trR[i];
        }
    }
}
