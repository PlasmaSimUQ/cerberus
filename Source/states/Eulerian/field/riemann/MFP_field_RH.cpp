#include "MFP_field_RH.H"
#include "MFP_utility.H"
#include "MFP.H"
#include "MFP_field.H"

//================================================================================

// J. Moreno, E. Oliva, P. Velarde, J.C.P. 2020, In Press

std::string FieldRH::tag = "RankineHugoniot";
bool FieldRH::registered = GetFieldRiemannSolverFactory().Register(FieldRH::tag, FieldRiemannSolverBuilder<FieldRH>);

FieldRH::FieldRH(){}
FieldRH::FieldRH(const int i)
{
    idx = i;

    c0 = MFP::lightspeed;
}

void FieldRH::solve(Array<Real,+FieldDef::ConsIdx::NUM> &L,
                    Array<Real,+FieldDef::ConsIdx::NUM> &R,
                    Array<Real,+FieldDef::ConsIdx::NUM> &F) const
{
    BL_PROFILE("FieldRH::solve");
    std::fill(F.begin(), F.end(), 0);

    Real Dr, Dl; // y- & z- D fields
    Real Br, Bl; // y- & z- B fields
    Real Pl, Pr; // div correction factors
    Real fl, fr; // x- fields
    Real epr, epl;
    Real mur, mul;
    Real cr, cl;

    // D-wave

    Dl = L[+FieldDef::ConsIdx::Dy];
    Dr = R[+FieldDef::ConsIdx::Dy];

    epr = R[+FieldDef::ConsIdx::ep];
    epl = L[+FieldDef::ConsIdx::ep];

    Bl = L[+FieldDef::ConsIdx::Bz];
    Br = R[+FieldDef::ConsIdx::Bz];

    mur = R[+FieldDef::ConsIdx::mu];
    mul = L[+FieldDef::ConsIdx::mu];

    cr = 1/std::sqrt(mur*epr);
    cl = 1/std::sqrt(mul*epl);

    F[+FieldDef::ConsIdx::Dy] = c0*((Bl*cl + Br*cr) + (Dl/epl - Dr/epr))/(cl*mul + cr*mur);
    F[+FieldDef::ConsIdx::Bz] = c0*((Dl*cl + Dr*cr) + (Bl/mul - Br/mur))/(cl*epl + cr*epr);

    // B-wave

    Dl = L[+FieldDef::ConsIdx::Dz];
    Dr = R[+FieldDef::ConsIdx::Dz];

    epr = R[+FieldDef::ConsIdx::ep];
    epl = L[+FieldDef::ConsIdx::ep];

    Bl = L[+FieldDef::ConsIdx::By];
    Br = R[+FieldDef::ConsIdx::By];

    F[+FieldDef::ConsIdx::Dz] = c0*(-(Bl*cl + Br*cr) + (Dl/epl - Dr/epr))/(cl*mul + cr*mur);
    F[+FieldDef::ConsIdx::By] = c0*(-(Dl*cl + Dr*cr) + (Bl/mul - Br/mur))/(cl*epl + cr*epr);


    return;
}

bool FieldRH::valid_state(const int idx)
{

    if (MFP::get_state(idx).get_type() != State::StateType::Field) {
        return false;
    }
    return true;
}
