#include "MFP_hydro_hlle2p.H"
#include "MFP_utility.H"
#include "MFP_global.H"

using GD = GlobalData;

//================================================================================

std::string HLLE2P::tag = "HLLE";
bool HLLE2P::registered = GetRiemannSolverFactory().Register(HLLE2P::tag, RiemannSolverBuilder<HLLE2P>);

HLLE2P::HLLE2P(){}
HLLE2P::HLLE2P(const int i)
{
    idx = i;
    istate = GD::get_state_ptr(i);
}

void HLLE2P::solve(Vector<Real> &L,
                   Vector<Real> &R,
                   Vector<Real> &F,
                   Real* shk)
{
    BL_PROFILE("HLLE2P::solve");

    // get the data out of the passed in arrays
    Real apL = L[+Hydro2PState::FluxIdx::Alpha];
    //    Real qL = GD::get_charge(apL, idx);
    Real gamL = L[+Hydro2PState::FluxIdx::Gamma];
    Real rhoL = L[+Hydro2PState::FluxIdx::Density];
    Real uL = L[+Hydro2PState::FluxIdx::Xvel];
    Real vL = L[+Hydro2PState::FluxIdx::Yvel];
    Real wL = L[+Hydro2PState::FluxIdx::Zvel];
    Real pL = L[+Hydro2PState::FluxIdx::Prs];
    Real paL = L[+Hydro2PState::FluxIdx::PrsP];
    Real peL = 0.5*(3*pL - paL);
    Real dpL = paL - peL;
    Real nrgL = pL/(gamL - 1) + 0.5*rhoL*(uL*uL + vL*vL + wL*wL);
    Real tL = apL*rhoL;
    Real aL = sqrt(gamL*pL/rhoL);

    Real BxL = L[+Hydro2PState::FluxIdx::Bx];
    Real ByL = L[+Hydro2PState::FluxIdx::By];
    Real BzL = L[+Hydro2PState::FluxIdx::Bz];

    Real BL = BxL*BxL + ByL*ByL + BzL*BzL;
    if (BL <= 0.0) {
        BL = 0.0;
    } else {
        BL = 1/std::sqrt(BL);
    }

    BxL *= BL;
    ByL *= BL;
    BzL *= BL;


    // get the data out of the passed in arrays
    Real apR = R[+Hydro2PState::FluxIdx::Alpha];
    //    Real qR = GD::get_charge(apR, idx);
    Real gamR = R[+Hydro2PState::FluxIdx::Gamma];
    Real rhoR = R[+Hydro2PState::FluxIdx::Density];
    Real uR = R[+Hydro2PState::FluxIdx::Xvel];
    Real vR = R[+Hydro2PState::FluxIdx::Yvel];
    Real wR = R[+Hydro2PState::FluxIdx::Zvel];
    Real pR = R[+Hydro2PState::FluxIdx::Prs];
    Real paR = R[+Hydro2PState::FluxIdx::PrsP];
    Real peR = 0.5*(3*pR - paR);
    Real dpR = paR - peR;
    Real nrgR = pR/(gamR-1) + 0.5*rhoR*(uR*uR + vR*vR + wR*wR);
    Real tR = apR*rhoR;
    Real aR = sqrt(gamR*pR/rhoR);

    Real BxR = R[+Hydro2PState::FluxIdx::Bx];
    Real ByR = R[+Hydro2PState::FluxIdx::By];
    Real BzR = R[+Hydro2PState::FluxIdx::Bz];

    Real BR = BxL*BxL + ByL*ByL + BzL*BzL;
    if (BR <= 0.0) {
        BR = 0.0;
    } else {
        BR = 1/std::sqrt(BR);
    }

    BxR *= BR;
    ByR *= BR;
    BzR *= BR;

    // speeds
    Real sL = std::min(uL - aL, uR - aR);
    Real sR = std::max(uL + aL, uR + aR);

    if (sL >= 0.0) {
        F[+Hydro2PState::ConsIdx::Density]  = rhoL*uL;
        F[+Hydro2PState::ConsIdx::Xmom]   = rhoL*uL*uL + BxL*BxL*dpL + peL;
        F[+Hydro2PState::ConsIdx::Ymom]   = rhoL*uL*vL + BxL*ByL*dpL;
        F[+Hydro2PState::ConsIdx::Zmom]   = rhoL*uL*wL + BxL*BzL*dpL;
        F[+Hydro2PState::ConsIdx::Eden]   = uL*(nrgL + peL) + BxL*dpL*(uL*BxL + vL*ByL + wL*BzL);
        F[+Hydro2PState::ConsIdx::PrsP]   = uL*paL;
        F[+Hydro2PState::ConsIdx::Tracer] = tL*uL;
    } else if ((sL <= 0.0) && (sR >= 0.0)) {

        Array<Real, +Hydro2PState::ConsIdx::NUM> fvL, fvR, svL, svR;

        // flux vector L
        fvL[+Hydro2PState::ConsIdx::Density]  = rhoL*uL;
        fvL[+Hydro2PState::ConsIdx::Xmom]   = rhoL*uL*uL + BxL*BxL*dpL + peL;
        fvL[+Hydro2PState::ConsIdx::Ymom]   = rhoL*uL*vL + BxL*ByL*dpL;
        fvL[+Hydro2PState::ConsIdx::Zmom]   = rhoL*uL*wL + BxL*BzL*dpL;
        fvL[+Hydro2PState::ConsIdx::Eden]   = uL*(nrgL + peL) + BxL*dpL*(uL*BxL + vL*ByL + wL*BzL);
        fvL[+Hydro2PState::ConsIdx::PrsP]   = uL*paL;
        fvL[+Hydro2PState::ConsIdx::Tracer] = tL*uL;

        // flux vector R
        fvR[+Hydro2PState::ConsIdx::Density]  = rhoR*uR;
        fvR[+Hydro2PState::ConsIdx::Xmom]   = rhoR*uR*uR + BxR*BxR*dpR + peR;
        fvR[+Hydro2PState::ConsIdx::Ymom]   = rhoR*uR*vR + BxR*ByR*dpR;
        fvR[+Hydro2PState::ConsIdx::Zmom]   = rhoR*uR*wR + BxR*BzR*dpR;
        fvR[+Hydro2PState::ConsIdx::Eden]   = uR*(nrgR + peR) + BxR*dpR*(uR*BxR + vR*ByR + wR*BzR);
        fvR[+Hydro2PState::ConsIdx::PrsP]   = uR*paR;
        fvR[+Hydro2PState::ConsIdx::Tracer] = tR*uR;

        // state vector L
        svL[+Hydro2PState::ConsIdx::Density]  = rhoL;
        svL[+Hydro2PState::ConsIdx::Xmom]   = rhoL*uL;
        svL[+Hydro2PState::ConsIdx::Ymom]   = rhoL*vL;
        svL[+Hydro2PState::ConsIdx::Zmom]   = rhoL*wL;
        svL[+Hydro2PState::ConsIdx::Eden]   = nrgL;
        svL[+Hydro2PState::ConsIdx::PrsP]   = paL;
        svL[+Hydro2PState::ConsIdx::Tracer] = tL;

        // state vector R
        svR[+Hydro2PState::ConsIdx::Density]  = rhoR;
        svR[+Hydro2PState::ConsIdx::Xmom]   = rhoR*uR;
        svR[+Hydro2PState::ConsIdx::Ymom]   = rhoR*vR;
        svR[+Hydro2PState::ConsIdx::Zmom]   = rhoR*wR;
        svR[+Hydro2PState::ConsIdx::Eden]   = nrgR;
        svR[+Hydro2PState::ConsIdx::PrsP]   = paR;
        svR[+Hydro2PState::ConsIdx::Tracer] = tR;


        for (int i=0; i<+Hydro2PState::ConsIdx::NUM; ++i) {
            F[i] = (sR*fvL[i] - sL*fvR[i] + sL*sR*(svR[i] - svL[i]))/(sR - sL);
        }
    } else {
        F[+Hydro2PState::ConsIdx::Density]  = rhoR*uR;
        F[+Hydro2PState::ConsIdx::Xmom]   = rhoR*uR*uR + BxR*BxR*dpR + peR;
        F[+Hydro2PState::ConsIdx::Ymom]   = rhoR*uR*vR + BxR*ByR*dpR;
        F[+Hydro2PState::ConsIdx::Zmom]   = rhoR*uR*wR + BxR*BzR*dpR;
        F[+Hydro2PState::ConsIdx::Eden]   = uR*(nrgR + peR) + BxR*dpR*(uR*BxR + vR*ByR + wR*BzR);
        F[+Hydro2PState::ConsIdx::PrsP]   = uR*paR;
        F[+Hydro2PState::ConsIdx::Tracer] = tR*uR;
    }
}

bool HLLE2P::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isHydro2P) {
        return false;
    }
    return true;
}
