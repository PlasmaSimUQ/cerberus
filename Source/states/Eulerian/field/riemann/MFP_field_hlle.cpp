#include "MFP_field_hlle.H"

#include "MFP.H"
#include "MFP_field.H"
#include "MFP_utility.H"

//================================================================================

std::string FieldHLLE::tag = "HLLE";
bool FieldHLLE::registered =
  GetFieldRiemannSolverFactory().Register(FieldHLLE::tag, FieldRiemannSolverBuilder<FieldHLLE>);

FieldHLLE::FieldHLLE() {}
FieldHLLE::FieldHLLE(const int i)
{
    BL_PROFILE("FieldHLLE::FieldHLLE");

    idx = i;

    FieldState& istate = static_cast<FieldState&>(MFP::get_state(i));

    c0 = MFP::lightspeed;
    Real c2 = c0 * c0;
    Real ch = istate.div_speed;
    Real ch2 = ch * ch;
    cc = ch2 / c2;
}

void FieldHLLE::solve(Vector<Real>& L, Vector<Real>& R, Vector<Real>& F, Real* shk)
{
    BL_PROFILE("FieldHLLE::solve");

    Array<Real, +FieldDef::ConsIdx::NUM> FL, FR;

    FL[+FieldDef::ConsIdx::Bx] = L[+FieldDef::ConsIdx::psi];
    FL[+FieldDef::ConsIdx::By] = -L[+FieldDef::ConsIdx::Dz] / L[+FieldDef::ConsIdx::ep];
    FL[+FieldDef::ConsIdx::Bz] = L[+FieldDef::ConsIdx::Dy] / L[+FieldDef::ConsIdx::ep];
    FL[+FieldDef::ConsIdx::Dx] = L[+FieldDef::ConsIdx::phi];
    FL[+FieldDef::ConsIdx::Dy] = L[+FieldDef::ConsIdx::Bz] / L[+FieldDef::ConsIdx::mu];
    FL[+FieldDef::ConsIdx::Dz] = -L[+FieldDef::ConsIdx::By] / L[+FieldDef::ConsIdx::mu];
    FL[+FieldDef::ConsIdx::psi] = L[+FieldDef::ConsIdx::Bx] * cc;
    FL[+FieldDef::ConsIdx::phi] = L[+FieldDef::ConsIdx::Dx] * cc;

    FR[+FieldDef::ConsIdx::Bx] = R[+FieldDef::ConsIdx::psi];
    FR[+FieldDef::ConsIdx::By] = -R[+FieldDef::ConsIdx::Dz] / R[+FieldDef::ConsIdx::ep];
    FR[+FieldDef::ConsIdx::Bz] = R[+FieldDef::ConsIdx::Dy] / R[+FieldDef::ConsIdx::ep];
    FR[+FieldDef::ConsIdx::Dx] = R[+FieldDef::ConsIdx::phi];
    FR[+FieldDef::ConsIdx::Dy] = R[+FieldDef::ConsIdx::Bz] / R[+FieldDef::ConsIdx::mu];
    FR[+FieldDef::ConsIdx::Dz] = -R[+FieldDef::ConsIdx::By] / R[+FieldDef::ConsIdx::mu];
    FR[+FieldDef::ConsIdx::psi] = R[+FieldDef::ConsIdx::Bx] * cc;
    FR[+FieldDef::ConsIdx::phi] = R[+FieldDef::ConsIdx::Dx] * cc;

    Real speed = c0;  // UPDATE THIS TO ACCOMMODATE DIV SPEED

    for (int n = 0; n < +FieldDef::ConsIdx::NUM; ++n) {
        F[n] = 0.5 * (FL[n] * c0 + FR[n] * c0 + (L[n] - R[n]) * speed);
    }
}

bool FieldHLLE::valid_state(const int idx)
{
    BL_PROFILE("FieldHLLE::valid_state");

    if (MFP::get_state(idx).get_type() != State::StateType::Field) { return false; }
    return true;
}
