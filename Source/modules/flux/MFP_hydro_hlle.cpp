#include "MFP_hydro_hlle.H"

#include "MFP_utility.H"
#include "MFP_global.H"

using GD = GlobalData;

//================================================================================
std::string HydroHLLE::tag = "HLLE";
bool HydroHLLE::registered = GetRiemannSolverFactory().Register(HydroHLLE::tag, RiemannSolverBuilder<HydroHLLE>);

HydroHLLE::HydroHLLE(){}
HydroHLLE::HydroHLLE(const int i)
{
    idx = i;

}

void HydroHLLE::solve(Vector<Real> &L,
                      Vector<Real> &R,
                      Vector<Real> &F,
                      Real* shk) const
{

    State &istate = GD::get_state(idx);

    // get the data out of the passed in arrays
    Real rhoL = L[+HydroState::PrimIdx::Density];
    Real uL = L[+HydroState::PrimIdx::Xvel];
    Real vL = L[+HydroState::PrimIdx::Yvel];
    Real wL = L[+HydroState::PrimIdx::Zvel];
    Real pL = L[+HydroState::PrimIdx::Prs];
    Real apL = L[+HydroState::PrimIdx::Alpha];

    Real gamL = istate.get_gamma(apL);
    Real nrgL = pL/(gamL - 1) + 0.5*rhoL*(uL*uL + vL*vL + wL*wL);
    Real tL = apL*rhoL;
    Real aL = sqrt(gamL*pL/rhoL);

    // get the data out of the passed in arrays
    Real rhoR = R[+HydroState::PrimIdx::Density];
    Real uR = R[+HydroState::PrimIdx::Xvel];
    Real vR = R[+HydroState::PrimIdx::Yvel];
    Real wR = R[+HydroState::PrimIdx::Zvel];
    Real pR = R[+HydroState::PrimIdx::Prs];
    Real apR = R[+HydroState::PrimIdx::Alpha];

    Real gamR = istate.get_gamma(apR);
    Real nrgR = pR/(gamR-1) + 0.5*rhoR*(uR*uR + vR*vR + wR*wR);
    Real tR = apR*rhoR;
    Real aR = sqrt(gamR*pR/rhoR);

    // flux vector L
    Vector<Real> fvL(+HydroState::ConsIdx::NUM);
    fvL[+HydroState::ConsIdx::Density]  = rhoL*uL;
    fvL[+HydroState::ConsIdx::Xmom]   = rhoL*uL*uL + pL;
    fvL[+HydroState::ConsIdx::Ymom]   = rhoL*uL*vL;
    fvL[+HydroState::ConsIdx::Zmom]   = rhoL*uL*wL;
    fvL[+HydroState::ConsIdx::Eden]   = uL*(nrgL + pL);
    fvL[+HydroState::ConsIdx::Tracer] = tL*uL;

    // flux vector R
    Vector<Real> fvR(+HydroState::ConsIdx::NUM);
    fvR[+HydroState::ConsIdx::Density]  = rhoR*uR;
    fvR[+HydroState::ConsIdx::Xmom]   = rhoR*uR*uR + pR;
    fvR[+HydroState::ConsIdx::Ymom]   = rhoR*uR*vR;
    fvR[+HydroState::ConsIdx::Zmom]   = rhoR*uR*wR;
    fvR[+HydroState::ConsIdx::Eden]   = uR*(nrgR + pR);
    fvR[+HydroState::ConsIdx::Tracer] = tR*uR;

    // state vector L
    Vector<Real> svL(+HydroState::ConsIdx::NUM);
    svL[+HydroState::ConsIdx::Density]  = rhoL;
    svL[+HydroState::ConsIdx::Xmom]   = rhoL*uL;
    svL[+HydroState::ConsIdx::Ymom]   = rhoL*vL;
    svL[+HydroState::ConsIdx::Zmom]   = rhoL*wL;
    svL[+HydroState::ConsIdx::Eden]   = nrgL;
    svL[+HydroState::ConsIdx::Tracer] = tL;

    // state vector R
    Vector<Real> svR(+HydroState::ConsIdx::NUM);
    svR[+HydroState::ConsIdx::Density]  = rhoR;
    svR[+HydroState::ConsIdx::Xmom]   = rhoR*uR;
    svR[+HydroState::ConsIdx::Ymom]   = rhoR*vR;
    svR[+HydroState::ConsIdx::Zmom]   = rhoR*wR;
    svR[+HydroState::ConsIdx::Eden]   = nrgR;
    svR[+HydroState::ConsIdx::Tracer] = tR;

    // speeds
    Real sL = std::min(uL - aL, uR - aR);
    Real sR = std::max(uL + aL, uR + aR);

    if (sL >= 0.0) {
        F = fvL;
    } else if ((sL <= 0.0) && (sR >= 0.0)) {
        for (int i=0; i<+HydroState::ConsIdx::NUM; ++i) {
            F[i] = (sR*fvL[i] - sL*fvR[i] + sL*sR*(svR[i] - svL[i]))/(sR - sL);
        }
    } else {
        F = fvR;
    }
}

bool HydroHLLE::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isHydro) {
        return false;
    }
    return true;
}
