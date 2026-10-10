#ifndef FUSION_H
#define FUSION_H

#include "MFP_collisions.H"

// =============================================================================
//  !!  WARNING: THE REACTION NETWORK AND FUEL/PRODUCT ACCOUNTING IN THIS  !!
//  !!  ACTION SUPPORT D + T -> He4 + n ONLY.  No D-D, D-He3, p-B11, no    !!
//  !!  secondary reactions, no burn-up of the product.  `fuel` must name  !!
//  !!  exactly two species and `product` exactly one; anything else is a  !!
//  !!  misuse that the code cannot detect.                                 !!
// =============================================================================
//
// D + T -> He4 (3.5 MeV) + n (14.1 MeV)
//
// Reactivity: NRL Plasma Formulary (J. D. Huba, Naval Research Laboratory,
// 2019 revision), "Fusion" section, Gamow-form fit valid for T <~ 25 keV:
//   <sv>_DT = 3.68e-12 T^-2/3 exp(-19.94 T^-1/3)  cm^3 s^-1  (T in keV)
//           = 3.68e-18 ...                       m^3 s^-1
// The standard higher-accuracy replacement (0.2-100 keV) is the
// parameterisation of H.-S. Bosch and G. M. Hale, Nucl. Fusion 32, 611 (1992);
// it is a drop-in change to the single reactivity expression in calc_update().
//
// Fuel accounting is the Collisions tracer write-back: D and T mass moves to
// the He4 sub-species; the fluid loses m_D + m_T - m_He4 per reaction (the
// neutron plus the Q/c^2 rest-mass defect, with nuclear masses).
//
// Alpha-particle energy is deposited locally (no transport, no range), split
// between the ion and electron populations either by a fixed `ion_fraction`
// or by the Fraley fit (`alpha_partition = 'fraley'`):
//   G. S. Fraley, E. J. Linnebur, R. J. Mason and R. L. Morse,
//   "Thermonuclear burn characteristics of compressed deuterium-tritium
//   microspheres", Phys. Fluids 17, 474 (1974):
//   f_i = 1 / (1 + fraley_constant / T_e[keV]),  fraley_constant = 32 (some
//   later authors quote 33). The neutron energy is lost from the system. The
// thermal energy of the consumed nuclei (~3kT per reaction) is neglected.
class Fusion : public Collisions
{
  public:
    Fusion();
    Fusion(const int idx, const sol::table& def);
    ~Fusion();

    virtual std::string get_tag() const override { return tag; }
    virtual ActionType get_type() const override { return ActionType::Collisions; }

    virtual void calc_update() override;

    static bool registered;
    static std::string tag;

    Real C_sv, kT_keV, E_alpha;  // code units, precomputed
    Real alpha_fraction, ion_fraction, max_burn_fraction;
    bool fraley;           // alpha_partition = 'fraley' (else fixed ion_fraction)
    Real fraley_constant;  // keV, 32 in the common fit
    int d_idx, t_idx, p_idx, e_idx;   // species indices (e_idx = -1 if none)
};

#endif  // FUSION_H
