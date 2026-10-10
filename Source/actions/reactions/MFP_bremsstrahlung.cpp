#include "MFP_bremsstrahlung.H"

#include "MFP.H"

std::string Bremsstrahlung::tag = "bremsstrahlung";
bool Bremsstrahlung::registered =
  GetActionFactory().Register(Bremsstrahlung::tag, ActionBuilder<Bremsstrahlung>);

Bremsstrahlung::Bremsstrahlung() {}
Bremsstrahlung::~Bremsstrahlung() {}

Bremsstrahlung::Bremsstrahlung(const int idx, const sol::table& def) : Collisions(idx, def)
{
    BL_PROFILE("Bremsstrahlung::Bremsstrahlung");

    const Real gaunt = def.get_or("gaunt", 1.0);
    max_sink_fraction = def.get_or("max_sink_fraction", 0.5);
    e_idx = electron_species();

    // P_code = P * t_ref / prs_ref with n = n_code * n_ref, T_eV = T_code * kT_eV
    const Real kT_eV = 1.38064852e-23 * MFP::T_ref / 1.602176634e-19;
    C_br = gaunt * 1.69e-38 * MFP::n_ref * MFP::n_ref * std::sqrt(kT_eV) * MFP::t_ref / MFP::prs_ref;
}

void Bremsstrahlung::calc_update()
{
    Real n_e = 0.0, Z2n = 0.0, T_e;
    int sink;

    if (e_idx >= 0) {
        n_e = species_info[e_idx].n;
        T_e = species_info[e_idx].T;
        sink = e_idx;
    } else {
        T_e = species_info[0].T;
        sink = 0;
    }

    for (const SpeciesInfo& s : species_info) {
        if (s.q > 0.0) {
            Z2n += s.q2 * s.n;
            if (e_idx < 0) n_e += s.q * s.n;  // quasi-neutral single fluid
        }
    }

    if ((n_e <= 0.0) || (T_e <= 0.0) || (Z2n <= 0.0)) return;

    Real P = C_br * n_e * std::sqrt(T_e) * Z2n;

    const Real E_th = 1.5 * n_e * T_e;  // electron thermal energy density (ideal-gas estimate)
    if (dt_stage > 0.0) P = std::min(P, max_sink_fraction * E_th / dt_stage);

    species_info[sink].delta[+HydroDef::ConsIdx::Eden] -= P;

    tau_min = std::min(tau_min, E_th / P);
}
