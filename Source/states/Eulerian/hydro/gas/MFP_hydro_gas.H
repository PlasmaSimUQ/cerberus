#ifndef MFP_HYDRO_GAS_H
#define MFP_HYDRO_GAS_H

#include "MFP.H"
#include "MFP_factory.H"
#include "MFP_hydro_defs.H"

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

    virtual std::string get_tag() const { return ""; }
    virtual void write_info(nlohmann::json& js) const;

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
