#include "MFP_mhd_hlle.H"
#include "MFP_utility.H"
#include "MFP_global.H"

using GD = GlobalData;

//================================================================================
std::string MhdHLLE::tag = "HLLE";
bool MhdHLLE::registered = GetRiemannSolverFactory().Register(MhdHLLE::tag, RiemannSolverBuilder<MhdHLLE>);

MhdHLLE::MhdHLLE(){}
MhdHLLE::MhdHLLE(const int i)
{
    idx = i;
}

void MhdHLLE::solve(Vector<Real> &L,
                    Vector<Real> &R,
                    Vector<Real> &F,
                    Real* shk) const
{
    BL_PROFILE("MhdHLLE::solve");
    State &istate = GD::get_state(idx);

    // get the data out of the passed in arrays
    Real rhoL = L[+MhdState::PrimIdx::Density];
    Real uL =  L[+MhdState::PrimIdx::Xvel];
    Real vL =  L[+MhdState::PrimIdx::Yvel];
    Real wL =  L[+MhdState::PrimIdx::Zvel];
    Real pL =   L[+MhdState::PrimIdx::Prs];
    Real apL =  L[+MhdState::PrimIdx::Alpha];
    Real BxL = L[+MhdState::PrimIdx::Bx];
    Real ByL = L[+MhdState::PrimIdx::By];
    Real BzL = L[+MhdState::PrimIdx::Bz];
    Real psiL = L[+MhdState::PrimIdx::psi];

    Real rhoR = R[+MhdState::PrimIdx::Density];
    Real uR =  R[+MhdState::PrimIdx::Xvel];
    Real vR =  R[+MhdState::PrimIdx::Yvel];
    Real wR =  R[+MhdState::PrimIdx::Zvel];
    Real pR =   R[+MhdState::PrimIdx::Prs];
    Real apR =  R[+MhdState::PrimIdx::Alpha];
    Real BxR = R[+MhdState::PrimIdx::Bx];
    Real ByR = R[+MhdState::PrimIdx::By];
    Real BzR = R[+MhdState::PrimIdx::Bz];
    Real psiR = R[+MhdState::PrimIdx::psi];

    // compute some stuff
    Real mL =   istate.get_mass(apL);
    Real gamL = istate.get_gamma(apL);

    Real trL = apL*rhoL;
    Real tL = pL/(rhoL/mL);
    Real RgasL = pL / (rhoL * tL);
    Real cvL = RgasL/(gamL - 1);
    Real rLsqrt = sqrt(rhoL);

    Real mR =   istate.get_mass(apR);
    Real gamR = istate.get_gamma(apR);

    Real trR = apR*rhoR;
    Real tR = pR/(rhoR/mR);
    Real RgasR = pR / (rhoR * tR);
    Real cvR = RgasR/(gamR - 1);
    Real rRsqrt = sqrt(rhoR);

    // combined properties

    Real alpha = rLsqrt / (rLsqrt + rRsqrt);

    Real cv = alpha * cvL + (1.0 - alpha) * cvR;
    Real Rgas = alpha * RgasL + (1.0 - alpha) * RgasR;

    Real cp = cv + Rgas;
    Real gam = cp / cv;

    /////////////////////////

    Real rL = rhoL;
    Real rR = rhoR;

    // average state
    Real rho = 0.5*(rL + rR);
    Real p   = 0.5*(pL + pR);
    Real u   = 0.5*(uL + uR);
    Real bx  = 0.5*(BxL + BxR);
    Real by  = 0.5*(ByL + ByR);
    Real bz  = 0.5*(BzL + BzR);

    // eigenvalues
    Real a2 = gam*p/rho;
    Real Bx2 = bx*bx;
    Real Bt2 = by*by + bz*bz;
    Real BB = Bx2 + Bt2;
    Real ca2 = Bx2/rho;
    Real alf = a2+BB/rho;
    Real als = std::sqrt(alf*alf-4*a2*ca2);
    Real cf2 = 0.5*(alf+als);
    Real cf = std::sqrt(cf2);

    Real wp = u + cf;
    Real wm = u - cf;

    // Compute the Jump in Conserved Variables between L and R
    Real BxL2 = BxL*BxL;
    Real BtL2 = ByL*ByL + BzL*BzL;
    Real BBL = BxL2 + BtL2;
    Real ptL = pL + 0.5*BBL;
    Real uL2 = uL*uL;
    Real uuL = uL2 + vL*vL + wL*wL;
    Real aL2 = gam*pL/rL;
    Real caL2 = BxL2/rL;
    Real alfL = aL2+BBL/rL;
    Real alsL = std::sqrt(alfL*alfL-4*aL2*caL2);
    Real cfL2 = 0.5*(alfL+alsL);
    Real cfL = std::sqrt(cfL2);
    Real wmL = uL-cfL;

    Real BxR2 = BxR*BxR;
    Real BtR2 = ByR*ByR + BzR*BzR;
    Real BBR = BxR2 + BtR2;
    Real ptR = pR + 0.5*BBR;
    Real uR2 = uR*uR;
    Real uuR = uR2 + vR*vR + wR*wR;
    Real aR2 = gam*pR/rR;
    Real caR2 = BxR2/rR;
    Real alfR = aR2+BBR/rR;
    Real alsR = std::sqrt(alfR*alfR-4*aR2*caR2);
    Real cfR2 = 0.5*(alfR+alsR);
    Real cfR = std::sqrt(cfR2);
    Real wpR = uR+cfR;

    Vector<Real> dU(+MhdState::ConsIdx::NUM);

    dU[+MhdState::ConsIdx::Density] = rR - rL;
    dU[+MhdState::ConsIdx::Xmom] = rR*uR - rL*uL;
    dU[+MhdState::ConsIdx::Ymom] = rR*vR - rL*vL;
    dU[+MhdState::ConsIdx::Zmom] = rR*wR - rL*wL;
    dU[+MhdState::ConsIdx::Eden] = (pR - pL)/(gam-1) + 0.5*(rR*uuR+BBR) - 0.5*(rL*uuL+BBL);
    dU[+MhdState::ConsIdx::Tracer] = trR - trL;
    dU[+MhdState::ConsIdx::Bx] = BxR - BxL;
    dU[+MhdState::ConsIdx::By] = ByR - ByL;
    dU[+MhdState::ConsIdx::Bz] = BzR - BzL;
    dU[+MhdState::ConsIdx::psi] = psiR - psiL;


    Real bl = std::min(wmL, wm);
    Real br = std::max(wpR, wp);
    Real blm = std::min(bl, 0.0);
    Real brp = std::max(br, 0.0);

    Real ch = istate.div_speed;
    blm = std::min(blm, -ch);
    brp = std::max(brp, ch);

    Real iden = 1.0/(brp - blm);
    Real fac1 = brp*blm;



    Vector<Real> fLU(+MhdState::ConsIdx::NUM, 0.0);

    fLU[+MhdState::ConsIdx::Density] = rL*uL;
    fLU[+MhdState::ConsIdx::Xmom] = rL*uL2 - BxL2 + ptL;
    fLU[+MhdState::ConsIdx::Ymom] = rL*uL*vL - BxL*ByL;
    fLU[+MhdState::ConsIdx::Zmom] = rL*uL*wL - BxL*BzL;
    fLU[+MhdState::ConsIdx::Eden] = (pL/(gam-1)+0.5*(rL*uuL+BBL)+ptL)*uL - (uL*BxL+vL*ByL+wL*BzL)*BxL;
    fLU[+MhdState::ConsIdx::Tracer] = trL*uL;
    fLU[+MhdState::ConsIdx::Bx]= psiL;
    fLU[+MhdState::ConsIdx::By]= uL*ByL - vL*BxL;
    fLU[+MhdState::ConsIdx::Bz] = uL*BzL - wL*BxL;
    fLU[+MhdState::ConsIdx::psi] = ch*ch*BxL;

    Vector<Real> fRU(+MhdState::ConsIdx::NUM);

    fRU[+MhdState::ConsIdx::Density] = rR*uR;
    fRU[+MhdState::ConsIdx::Xmom] = rR*uR2 - BxR2 + ptR;
    fRU[+MhdState::ConsIdx::Ymom] = rR*uR*vR - BxR*ByR;
    fRU[+MhdState::ConsIdx::Zmom] = rR*uR*wR - BxR*BzR;
    fRU[+MhdState::ConsIdx::Eden] = (pR/(gam-1)+0.5*(rR*uuR+BBR)+ptR)*uR - (uR*BxR+vR*ByR+wR*BzR)*BxR;
    fRU[+MhdState::ConsIdx::Tracer] = trR*uR;
    fRU[+MhdState::ConsIdx::Bx]= psiR;
    fRU[+MhdState::ConsIdx::By] = uR*ByR - vR*BxR;
    fRU[+MhdState::ConsIdx::Bz] = uR*BzR - wR*BxR;
    fRU[+MhdState::ConsIdx::psi] = ch*ch*BxR;

    for (int i=0; i<+MhdState::ConsIdx::NUM; ++i) {
        F[i] = (brp*fLU[i] - blm*fRU[i] + fac1*dU[i])*iden;
    }
}

bool MhdHLLE::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isMHD) {
        return false;
    }
    return true;
}
