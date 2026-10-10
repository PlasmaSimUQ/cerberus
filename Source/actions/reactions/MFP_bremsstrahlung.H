#ifndef BREMSSTRAHLUNG_H
#define BREMSSTRAHLUNG_H

#include "MFP_collisions.H"

// Optically-thin, non-relativistic bremsstrahlung loss.
//
// Reference: NRL Plasma Formulary (J. D. Huba, Naval Research Laboratory,
// 2019 revision), "Radiation" section, total bremsstrahlung power from a
// Maxwellian plasma (Gaunt factor 1):
//   P_br = 1.69e-32 n_e[cm^-3] sqrt(T_e[eV]) sum_i Z_i^2 n_i[cm^-3]   W cm^-3
//        = 1.69e-38 n_e        sqrt(T_e[eV]) sum_i Z_i^2 n_i          W m^-3
// The optional `gaunt` multiplier stands in for the frequency-averaged
// free-free Gaunt factor (W. J. Karzas and R. Latter, Astrophys. J. Suppl. 6,
// 167 (1961)); ~1.2 is the usual thermal value. No reabsorption, no
// relativistic correction, no line/recombination radiation.
//
// Applied as an energy sink on the electron species (q < 0) if one is listed,
// otherwise on the single-fluid state with n_e from quasi-neutrality.
class Bremsstrahlung : public Collisions
{
  public:
    Bremsstrahlung();
    Bremsstrahlung(const int idx, const sol::table& def);
    ~Bremsstrahlung();

    virtual std::string get_tag() const override { return tag; }
    virtual ActionType get_type() const override { return ActionType::Collisions; }

    virtual void calc_update() override;

    static bool registered;
    static std::string tag;

    Real C_br;               // code-unit prefactor, precomputed
    Real max_sink_fraction;  // per-stage cap on E_th removed (positivity)
    int e_idx;               // electron species index or -1
};

#endif  // BREMSSTRAHLUNG_H
