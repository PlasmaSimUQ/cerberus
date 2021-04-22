#include "MFP_field_RH.H"
#include "MFP_utility.H"
#include "MFP_global.H"

using GD = GlobalData;

//================================================================================

// J. Moreno, E. Oliva, P. Velarde, J.C.P. 2020, In Press

std::string FieldRH::tag = "RankineHugoniot";
bool FieldRH::registered = GetRiemannSolverFactory().Register(FieldRH::tag, RiemannSolverBuilder<FieldRH>);

FieldRH::FieldRH(){}
FieldRH::FieldRH(const int i)
{
    idx = i;
}

void FieldRH::solve(Vector<Real> &L,
                    Vector<Real> &R,
                    Vector<Real> &F,
                    Real* shk) const
{
    BL_PROFILE("FieldRH::solve");
    std::fill(F.begin(), F.end(), 0);

    State &istate = GD::get_state(idx);

    Real c0 = GD::lightspeed;
    Real ch = istate.div_speed;
    Real ch2 = ch*ch;

    Real Dr, Dl; // y- & z- D fields
    Real Br, Bl; // y- & z- B fields
    Real Pl, Pr; // div correction factors
    Real fl, fr; // x- fields
    Real epr, epl;
    Real mur, mul;
    Real cr, cl;

    // D-wave

    Dl = L[+FieldState::PrimIdx::Dy];
    Dr = R[+FieldState::PrimIdx::Dy];

    epr = R[+FieldState::PrimIdx::ep];
    epl = L[+FieldState::PrimIdx::ep];

    Bl = L[+FieldState::PrimIdx::Bz];
    Br = R[+FieldState::PrimIdx::Bz];

    mur = R[+FieldState::PrimIdx::mu];
    mul = L[+FieldState::PrimIdx::mu];

    cr = 1/std::sqrt(mur*epr);
    cl = 1/std::sqrt(mul*epl);

    F[+FieldState::PrimIdx::Dy] = c0*((Bl*cl + Br*cr) + (Dl/epl - Dr/epr))/(cl*mul + cr*mur);
    F[+FieldState::PrimIdx::Bz] = c0*((Dl*cl + Dr*cr) + (Bl/mul - Br/mur))/(cl*epl + cr*epr);


    // div clean
    if (ch > 0) {
        Pl = L[+FieldState::PrimIdx::phi];
        Pr = R[+FieldState::PrimIdx::phi];

        fl = L[+FieldState::PrimIdx::Dx];
        fr = R[+FieldState::PrimIdx::Dx];

        F[+FieldState::PrimIdx::Dx] =   0.5*c0*(Pr + Pl)  - 0.5*ch*(fr - fl);
        F[+FieldState::PrimIdx::phi] = 0.5*ch2/c0*(fr + fl) - 0.5*ch*(Pr - Pl);
    }

    // B-wave

    Dl = L[+FieldState::PrimIdx::Dz];
    Dr = R[+FieldState::PrimIdx::Dz];

    epr = R[+FieldState::PrimIdx::ep];
    epl = L[+FieldState::PrimIdx::ep];

    Bl = L[+FieldState::PrimIdx::By];
    Br = R[+FieldState::PrimIdx::By];

    F[+FieldState::PrimIdx::Dz] = c0*(-(Bl*cl + Br*cr) + (Dl/epl - Dr/epr))/(cl*mul + cr*mur);
    F[+FieldState::PrimIdx::By] = c0*(-(Dl*cl + Dr*cr) + (Bl/mul - Br/mur))/(cl*epl + cr*epr);

    // div clean
    if (ch > 0) {
        Pl = L[+FieldState::PrimIdx::psi];
        Pr = R[+FieldState::PrimIdx::psi];

        fl = L[+FieldState::PrimIdx::Bx];
        fr = R[+FieldState::PrimIdx::Bx];

        F[+FieldState::PrimIdx::Bx] =   0.5*c0*(Pr + Pl)  - 0.5*ch*(fr - fl);
        F[+FieldState::PrimIdx::psi] = 0.5*ch2/c0*(fr + fl) - 0.5*ch*(Pr - Pl);
    }


    return;
}

bool FieldRH::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isField) {
        return false;
    }
    return true;
}
