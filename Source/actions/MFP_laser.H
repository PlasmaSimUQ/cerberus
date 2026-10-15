#ifndef LASER_H
#define LASER_H

#include "MFP_action.H"
#include "MFP_hydro.H"
#include "MFP_optional_func.H"

// Laser energy deposition by inverse bremsstrahlung along straight,
// unrefracted rays (geometric optics without the bending term).
//
// A beam is a lattice of parallel rays from `origin` along `direction`, with
// a `radius` cross-section and a 'tophat' or 'gaussian' profile; each ray
// carries P_k = I(t) * profile * dA_k. Along a ray dP/ds = -kappa P with
//   kappa = (nu_ei / c) (n_e / n_c) (1 - n_e / n_c)^-1/2
//   nu_ei = 2.91e-6 sum_i Z_i^2 n_i[cm^-3] lnLambda T_e[eV]^-3/2   s^-1
//   n_c   = eps0 m_e omega^2 / e^2
// (NRL Plasma Formulary, J. D. Huba, NRL 2019 revision: "Collisions and
// transport" for nu_ei and the electron-ion Coulomb logarithm, "Laser-plasma"
// for the critical density; the absorption coefficient is the standard
// collisional-damping form, e.g. W. L. Kruer, "The Physics of Laser Plasma
// Interactions", Addison-Wesley 1988, ch. 5). The power absorbed in each
// cell is deposited as electron heating; whatever reaches the first cell
// with n_e >= critical_fraction * n_c is deposited there (no reflection).
//
// Rays are marched cell by cell with the Amanatides-Woo traversal
// (J. Amanatides and A. Woo, "A fast voxel traversal algorithm for ray
// tracing", Eurographics 1987); kappa is piecewise constant per cell. The
// regular transverse lattice (spacing dx_min / rays_per_cell) keeps the
// deposition smooth; cf. the inverse-projection launch of HYDRA.
//
// Not modelled: refraction, resonance absorption, cross-beam energy
// transfer, the Langdon effect. Level 0 only (aborts if amr.max_level > 0);
// Cartesian only. The level is gathered to rank 0 for the march.
class Laser : public Action
{
  public:
    Laser();
    Laser(const int idx, const sol::table& def);
    ~Laser();

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;
    virtual void calc_time_derivative(MFP* mfp,
                                      Vector<UpdateData>& update,
                                      const Real time,
                                      const Real dt) override;

    virtual std::string get_tag() const override { return tag; }
    virtual ActionType get_type() const override { return ActionType::Laser; }

    static bool registered;
    static std::string tag;

  private:
    void build_lattice(const Real* dx);
    void trace(Array<Real, 3> o,
               Real P,
               const Array4<const Real>& f4,
               const Array4<Real>& d4,
               const Real* dx,
               const Real* plo,
               const Real* phi,
               const Box& dom,
               Real& absorbed,
               Real& dumped,
               Real& transmitted) const;

    HydroState* state;
    bool electron_state;  // every component of `state` has q < 0 (electron fluid)
    Real Z_eff;           // ion charge used in nu_ei and lnLambda

    Array<Real, 3> origin, dir;  // code units; unit vector
    Real radius, extent;         // beam radius and lattice half-extent, code units
    bool gaussian;
    int rays_per_cell;
    Vector<std::pair<Array<Real, 3>, Real>> lattice;  // (offset, profile * dA)

    Optional3D1VFunction intensity;  // W m^-2 on axis, function of t
    Real I_scale;                    // -> code power per code area

    Real n_c, crit_frac, C_kappa, lnLambda, kT_eV, n_ref_cm3;
    bool lnLambda_nrl;
};

#endif  // LASER_H
