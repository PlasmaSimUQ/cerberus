#include "MFP_field_hlle.H"
#include "MFP_utility.H"
#include "MFP_global.H"

using GD = GlobalData;

//================================================================================

std::string FieldHLLE::tag = "HLLE";
bool FieldHLLE::registered = GetRiemannSolverFactory().Register(FieldHLLE::tag, RiemannSolverBuilder<FieldHLLE>);

FieldHLLE::FieldHLLE(){}
FieldHLLE::FieldHLLE(const int i)
{
    idx = i;
    istate = GD::get_state_ptr(i);
}

void FieldHLLE::solve(Vector<Real> &L,
                      Vector<Real> &R,
                      Vector<Real> &F,
                      Real* shk)
{
    BL_PROFILE("FieldHLLE::solve");

    // vvv this could all be static, maybe move to init?? vvv
    Real c0 = GD::lightspeed;
    Real c2 = c0*c0;
    Real ch = istate->div_speed;
    Real ch2 = ch*ch;

    Real cc = ch2/c2;
    // ^^^                                                ^^^

    Array<Real, +FieldState::ConsIdx::NUM> FL, FR;

    FL[+FieldState::ConsIdx::Bx]   =   L[+FieldState::PrimIdx::psi];
    FL[+FieldState::ConsIdx::By]   = - L[+FieldState::PrimIdx::Dz]/L[+FieldState::PrimIdx::ep];
    FL[+FieldState::ConsIdx::Bz]   =   L[+FieldState::PrimIdx::Dy]/L[+FieldState::PrimIdx::ep];
    FL[+FieldState::ConsIdx::Dx]   =   L[+FieldState::PrimIdx::phi];
    FL[+FieldState::ConsIdx::Dy]   =   L[+FieldState::PrimIdx::Bz]/L[+FieldState::PrimIdx::mu];
    FL[+FieldState::ConsIdx::Dz]   = - L[+FieldState::PrimIdx::By]/L[+FieldState::PrimIdx::mu];
    FL[+FieldState::ConsIdx::psi] =   L[+FieldState::PrimIdx::Bx]*cc;
    FL[+FieldState::ConsIdx::phi] =   L[+FieldState::PrimIdx::Dx]*cc;

    FR[+FieldState::ConsIdx::Bx]   =   R[+FieldState::PrimIdx:: psi];
    FR[+FieldState::ConsIdx::By]   = - R[+FieldState::PrimIdx:: Dz]/R[+FieldState::PrimIdx::ep];
    FR[+FieldState::ConsIdx::Bz]   =   R[+FieldState::PrimIdx:: Dy]/R[+FieldState::PrimIdx::ep];
    FR[+FieldState::ConsIdx::Dx]   =   R[+FieldState::PrimIdx:: phi];
    FR[+FieldState::ConsIdx::Dy]   =   R[+FieldState::PrimIdx:: Bz]/R[+FieldState::PrimIdx::mu];
    FR[+FieldState::ConsIdx::Dz]   = - R[+FieldState::PrimIdx:: By]/R[+FieldState::PrimIdx::mu];
    FR[+FieldState::ConsIdx::psi] =   R[+FieldState::PrimIdx::Bx]*cc;
    FR[+FieldState::ConsIdx::phi] =   R[+FieldState::PrimIdx::Dx]*cc;

    RealArray speed = istate->get_speed_from_prim(L); // only care about speed in x

    for (int n=0; n<+FieldState::ConsIdx::NUM; ++n) {
        F[n] = 0.5*(FL[n]*c0 + FR[n]*c0 + (L[n] - R[n])*speed[0]);
    }

}

bool FieldHLLE::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isField) {
        return false;
    }
    return true;
}
