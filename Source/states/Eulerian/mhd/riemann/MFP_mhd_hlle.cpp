#include "MFP_mhd_hlle.H"

#include "MFP.H"
#include "MFP_state.H"
#include "MFP_utility.H"

//================================================================================
std::string MHDHLLE::tag = "HLLE";
bool MHDHLLE::registered =
  GetMHDRiemannSolverFactory().Register(MHDHLLE::tag, MHDRiemannSolverBuilder<MHDHLLE>);

MHDHLLE::MHDHLLE() {}
MHDHLLE::MHDHLLE(const sol::table& def)
{
    BL_PROFILE("MHDHLLE::MHDHLLE");

    gamma = def["gamma"];
    div_transport = def["div_transport"];
}

void MHDHLLE::solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk)
{
    BL_PROFILE("MHDHLLE::solve");

    // get the data out of the passed in arrays
    Real rL = L[+MHDDef::PrimIdx::Density];
    Real uL = L[+MHDDef::PrimIdx::Xvel];
    Real vL = L[+MHDDef::PrimIdx::Yvel];
    Real wL = L[+MHDDef::PrimIdx::Zvel];
    Real pL = L[+MHDDef::PrimIdx::Prs];
    Real BxL = L[+MHDDef::PrimIdx::Bx];
    Real ByL = L[+MHDDef::PrimIdx::By];
    Real BzL = L[+MHDDef::PrimIdx::Bz];
    Real psiL = L[+MHDDef::PrimIdx::psi];

    Real rR = R[+MHDDef::PrimIdx::Density];
    Real uR = R[+MHDDef::PrimIdx::Xvel];
    Real vR = R[+MHDDef::PrimIdx::Yvel];
    Real wR = R[+MHDDef::PrimIdx::Zvel];
    Real pR = R[+MHDDef::PrimIdx::Prs];
    Real BxR = R[+MHDDef::PrimIdx::Bx];
    Real ByR = R[+MHDDef::PrimIdx::By];
    Real BzR = R[+MHDDef::PrimIdx::Bz];
    Real psiR = R[+MHDDef::PrimIdx::psi];

    /////////////////////////

    // average state
    Real rho = 0.5 * (rL + rR);
    Real p = 0.5 * (pL + pR);
    Real u = 0.5 * (uL + uR);
    Real bx = 0.5 * (BxL + BxR);
    Real by = 0.5 * (ByL + ByR);
    Real bz = 0.5 * (BzL + BzR);

    // eigenvalues
    Real a2 = gamma * p / rho;
    Real Bx2 = bx * bx;
    Real Bt2 = by * by + bz * bz;
    Real BB = Bx2 + Bt2;
    Real ca2 = Bx2 / rho;
    Real alf = a2 + BB / rho;
    Real als = std::sqrt(alf * alf - 4 * a2 * ca2);
    Real cf2 = 0.5 * (alf + als);
    Real cf = std::sqrt(cf2);

    Real wp = u + cf;
    Real wm = u - cf;

    // Compute the Jump in Conserved Variables between L and R
    Real BxL2 = BxL * BxL;
    Real BtL2 = ByL * ByL + BzL * BzL;
    Real BBL = BxL2 + BtL2;
    Real ptL = pL + 0.5 * BBL;
    Real uL2 = uL * uL;
    Real uuL = uL2 + vL * vL + wL * wL;
    Real aL2 = gamma * pL / rL;
    Real caL2 = BxL2 / rL;
    Real alfL = aL2 + BBL / rL;
    Real alsL = std::sqrt(alfL * alfL - 4 * aL2 * caL2);
    Real cfL2 = 0.5 * (alfL + alsL);
    Real cfL = std::sqrt(cfL2);
    Real wmL = uL - cfL;

    Real BxR2 = BxR * BxR;
    Real BtR2 = ByR * ByR + BzR * BzR;
    Real BBR = BxR2 + BtR2;
    Real ptR = pR + 0.5 * BBR;
    Real uR2 = uR * uR;
    Real uuR = uR2 + vR * vR + wR * wR;
    Real aR2 = gamma * pR / rR;
    Real caR2 = BxR2 / rR;
    Real alfR = aR2 + BBR / rR;
    Real alsR = std::sqrt(alfR * alfR - 4 * aR2 * caR2);
    Real cfR2 = 0.5 * (alfR + alsR);
    Real cfR = std::sqrt(cfR2);
    Real wpR = uR + cfR;

    Array<Real, +MHDDef::ConsIdx::NUM> dU, fLU, fRU;

    dU[+MHDDef::ConsIdx::Density] = rR - rL;
    dU[+MHDDef::ConsIdx::Xmom] = rR * uR - rL * uL;
    dU[+MHDDef::ConsIdx::Ymom] = rR * vR - rL * vL;
    dU[+MHDDef::ConsIdx::Zmom] = rR * wR - rL * wL;
    dU[+MHDDef::ConsIdx::Eden] =
      (pR - pL) / (gamma - 1) + 0.5 * (rR * uuR + BBR) - 0.5 * (rL * uuL + BBL);
    dU[+MHDDef::ConsIdx::Bx] = BxR - BxL;
    dU[+MHDDef::ConsIdx::By] = ByR - ByL;
    dU[+MHDDef::ConsIdx::Bz] = BzR - BzL;
    dU[+MHDDef::ConsIdx::psi] = psiR - psiL;

    Real bl = std::min(wmL, wm);
    Real br = std::max(wpR, wp);
    Real blm = std::min(bl, 0.0);
    Real brp = std::max(br, 0.0);

    Real div_speed = div_transport * max_speed;

    blm = std::min(blm, -div_speed);
    brp = std::max(brp, div_speed);

    Real iden = 1.0 / (brp - blm);
    Real fac1 = brp * blm;

    div_speed *= div_speed;

    fLU[+MHDDef::ConsIdx::Density] = rL * uL;
    fLU[+MHDDef::ConsIdx::Xmom] = rL * uL2 - BxL2 + ptL;
    fLU[+MHDDef::ConsIdx::Ymom] = rL * uL * vL - BxL * ByL;
    fLU[+MHDDef::ConsIdx::Zmom] = rL * uL * wL - BxL * BzL;
    fLU[+MHDDef::ConsIdx::Eden] = (pL / (gamma - 1) + 0.5 * (rL * uuL + BBL) + ptL) * uL -
                                  (uL * BxL + vL * ByL + wL * BzL) * BxL;
    fLU[+MHDDef::ConsIdx::Bx] = psiL;
    fLU[+MHDDef::ConsIdx::By] = uL * ByL - vL * BxL;
    fLU[+MHDDef::ConsIdx::Bz] = uL * BzL - wL * BxL;
    fLU[+MHDDef::ConsIdx::psi] = div_speed * BxL;

    fRU[+MHDDef::ConsIdx::Density] = rR * uR;
    fRU[+MHDDef::ConsIdx::Xmom] = rR * uR2 - BxR2 + ptR;
    fRU[+MHDDef::ConsIdx::Ymom] = rR * uR * vR - BxR * ByR;
    fRU[+MHDDef::ConsIdx::Zmom] = rR * uR * wR - BxR * BzR;
    fRU[+MHDDef::ConsIdx::Eden] = (pR / (gamma - 1) + 0.5 * (rR * uuR + BBR) + ptR) * uR -
                                  (uR * BxR + vR * ByR + wR * BzR) * BxR;
    fRU[+MHDDef::ConsIdx::Bx] = psiR;
    fRU[+MHDDef::ConsIdx::By] = uR * ByR - vR * BxR;
    fRU[+MHDDef::ConsIdx::Bz] = uR * BzR - wR * BxR;
    fRU[+MHDDef::ConsIdx::psi] = div_speed * BxR;

    for (int i = 0; i < +MHDDef::ConsIdx::NUM; ++i) {
        F[i] = (brp * fLU[i] - blm * fRU[i] + fac1 * dU[i]) * iden;
    }
}
