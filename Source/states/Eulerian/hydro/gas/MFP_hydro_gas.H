#ifndef MFP_HYDRO_GAS_H
#define MFP_HYDRO_GAS_H

#include "MFP.H"
#include "MFP_factory.H"
#include "MFP_hydro_defs.H"

#include <cmath>

//-------------------------------------------------------------------------

class HydroGas
{
  public:
    HydroGas();
    virtual ~HydroGas();

    size_t n_species() const { return mass.size(); }
    size_t n_tracers() const { return n_species() - 1; }

    void get_alpha_fractions_from_prim(const Vector<Real>& Q,
                                       Vector<Real>& alpha,
                                       const int tracer_idx = +HydroDef::PrimIdx::NUM) const;
    void get_alpha_fractions_from_cons(const Vector<Real>& U,
                                       Vector<Real>& alpha,
                                       const int density_idx = +HydroDef::ConsIdx::Density,
                                       const int tracer_idx = +HydroDef::ConsIdx::NUM) const;

    Real get_mass_from_cons(const Vector<Real>& U,
                            const int density_idx = +HydroDef::ConsIdx::Density,
                            const int tracer_idx = +HydroDef::ConsIdx::NUM) const;
    Real get_charge_from_cons(const Vector<Real>& U,
                              const int density_idx = +HydroDef::ConsIdx::Density,
                              const int tracer_idx = +HydroDef::ConsIdx::NUM) const;

    virtual Real get_gamma_from_cons(const Vector<Real>& U,
                                     const int density_idx = +HydroDef::ConsIdx::Density,
                                     const int tracer_idx = +HydroDef::ConsIdx::NUM) const
    {
        return 0;
    }
    virtual Real get_cp_from_cons(const Vector<Real>& U,
                                  const int density_idx = +HydroDef::ConsIdx::Density,
                                  const int tracer_idx = +HydroDef::ConsIdx::NUM) const
    {
        return 0;
    }

    Real get_mass_from_prim(const Vector<Real>& Q, const int idx = +HydroDef::PrimIdx::NUM) const;
    Real get_charge_from_prim(const Vector<Real>& Q,
                              const int tracer_idx = +HydroDef::PrimIdx::NUM) const;
    virtual Real get_gamma_from_prim(const Vector<Real>& Q,
                                     const int idx = +HydroDef::PrimIdx::NUM) const
    {
        return 0;
    }
    virtual Real get_cp_from_prim(const Vector<Real>& Q,
                                  const int tracer_idx = +HydroDef::PrimIdx::NUM) const
    {
        return 0;
    }

    virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const { return true; }
    virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const {}

    // one-shot warm start for the NEXT cons2prim call — consumed and reset
    // on entry by models that use it, so unwired call sites stay cold. No-op
    // for closed-form gas models.
    virtual void set_cons2prim_T_seed(const Real T) const {}

    virtual void define_rho_p_T(Vector<Real>& Q) const {}
    Real get_density_from_number_density(const Real nd, Vector<Real>& Q) const;

    bool prim_valid(const Vector<Real>& Q) const;
    bool cons_valid(const Vector<Real>& U) const;
    Real get_energy_from_cons(const Vector<Real>& U) const;
    Real get_internal_energy_density_from_cons(const Vector<Real>& U) const;
    virtual Real get_temperature_from_cons(const Vector<Real>& U) const { return 0; }
    Real get_temperature_from_prim(const Vector<Real>& Q) const;
    virtual RealArray get_speed_from_cons(const Vector<Real>& U) const
    {
        return {AMREX_D_DECL(0, 0, 0)};
    }
    virtual RealArray get_speed_from_prim(const Vector<Real>& Q) const
    {
        return {AMREX_D_DECL(0, 0, 0)};
    }

    // Face evaluations for the general-EOS Riemann solver (plan W10,
    // Stage 5 — STAGE5.md D-a/D-b). `Q` is a reconstructed face primitive
    // vector in the HydroDef layout. The defaults below read the
    // reconstructed `Gamma` SLOT (not get_gamma_from_prim, which would
    // recompute gamma from the reconstructed tracers) so they reproduce the
    // algebra of MFP_hllc.cpp exactly — HLLC_general_eos on a gamma-law gas
    // then agrees with HLLC to round-off, which is the gate that isolates
    // solver bugs from EOS bugs.
    virtual Real get_internal_energy_from_prim(const Vector<Real>& Q) const
    {
        return Q[+HydroDef::PrimIdx::Prs] /
               ((Q[+HydroDef::PrimIdx::Gamma] - 1.0) * Q[+HydroDef::PrimIdx::Density]);
    }
    virtual Real get_sound_speed_from_prim_rp(const Vector<Real>& Q) const
    {
        return std::sqrt(Q[+HydroDef::PrimIdx::Gamma] * Q[+HydroDef::PrimIdx::Prs] /
                         Q[+HydroDef::PrimIdx::Density]);
    }
    // combined face evaluation: the solver needs BOTH e and a per face, and
    // for a tabulated gas they come from the same (rho, p) table inversion —
    // the override answers both for one inversion (the G8 wall-time budget
    // depends on this). The default is exactly the two calls above, so the
    // gamma-law path is unchanged.
    virtual void get_face_eval_from_prim(const Vector<Real>& Q, Real& e, Real& a) const
    {
        e = get_internal_energy_from_prim(Q);
        a = get_sound_speed_from_prim_rp(Q);
    }

    virtual std::string get_tag() const { return ""; }
    virtual void write_info(nlohmann::json& js) const;

    // true for gas models whose closure is not a gamma-law (tabulated,
    // tabulated_mixture): HydroState::set_flux defaults these to the
    // HLLC_general_eos solver, and the MHD compatibility guard keys off the
    // same predicate (plan W24 generalising W10/W14 — the tag-literal
    // comparisons would silently miss new general-EOS gas types).
    virtual bool needs_general_eos_solver() const { return false; }

    // clamp face density/pressure ahead of a Riemann solve; returns the
    // number of components that were clamped (all hydro gases share the
    // HydroDef primitive layout so this lives on the base). Virtual so a
    // gas model with a *physical* floor can override the absolute
    // effective_zero: the tabulated gas floors to its table-hull edge
    // (plan W8.1, Stage 4), which keeps the floored sound speed finite.
    virtual int apply_prim_floor(Vector<Real>& Q) const
    {
#ifdef MFP_PRIM_FLOOR
        int n = 0;
        if (Q[+HydroDef::PrimIdx::Density] < effective_zero) {
            Q[+HydroDef::PrimIdx::Density] = effective_zero;
            n += 1;
        }
        if (Q[+HydroDef::PrimIdx::Prs] < effective_zero) {
            Q[+HydroDef::PrimIdx::Prs] = effective_zero;
            n += 1;
        }
        return n;
#else
        return 0;
#endif
    }

    // diagnostic: number of table-hull clamps since the last call (always 0
    // for analytic gas models; the tabulated gas counts its inversion
    // clamps/non-convergences here — plan W8.3, Stage 4)
    virtual long get_and_reset_hull_clamps() const { return 0; }

    // composition repairs (multi-component gas models only): evaluations
    // whose raw tracer fractions needed clamping/renormalising into the
    // simplex before use -- silent for single-component models
    virtual long get_and_reset_alpha_fixes() const { return 0; }

    // B (T4b plan): floor-time momentum repairs by cons2prim, and the
    // per-call indicator that the last cons2prim modified U in place (the
    // caller then persists U to the conserved fab). Base gas models never
    // repair -- both are inert defaults.
    virtual long get_and_reset_mom_fixes() const { return 0; }
    virtual bool cons_was_repaired() const { return false; }

    // TEMPORARY DIAGNOSTIC (W31 investigation): per-component pin accounting
    // for the mixture assembly. Empty string = nothing to report, which is
    // every model except MixtureEOS under the Amagat rule.
    virtual std::string get_and_reset_pin_report() const { return {}; }

    int idx;
    bool mass_const, charge_const;
    Vector<Real> mass, charge;
    Vector<std::string> comp_names;

    // single source of truth for the positivity floor, mirroring the state's
    // Lua-configurable EulerianState::effective_zero (set by HydroGasBuilder)
    Real effective_zero = 1e-14;
};

template <typename D> std::unique_ptr<HydroGas> HydroGasBuilder(const sol::table& def)
{
    const int global_idx = def["global_idx"];
    const sol::table& gas = def["gas"].get_or(sol::table());

    if (gas.valid()) {
        // load the diffusion coefficients

        if (gas["type"] == D::tag) {
            std::unique_ptr<D> obj(new D(global_idx, gas));
            // the ctor only sees the 'gas' sub-table; the floor threshold is a
            // state-level key so it is picked up here for every gas type
            obj->effective_zero = def["effective_zero"].get_or(1e-14);
            return obj;
        } else {
            return nullptr;
        }
    } else {
        return nullptr;
    }
}

ClassFactory<HydroGas>& GetHydroGasFactory();

#endif  // MFP_HYDRO_GAS_H
