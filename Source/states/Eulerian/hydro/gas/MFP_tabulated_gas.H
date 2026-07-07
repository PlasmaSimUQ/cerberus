#ifndef MFP_TABULATED_GAS_H
#define MFP_TABULATED_GAS_H

// Tabulated-EOS gas model (plan W6, Stage 3 — doc/eos_implementation_plan.md;
// user documentation in doc/eos_table_reader.md §10).
//
// A HydroGas backend whose thermodynamic closure comes from a conditioned
// `.eostab` table (EosTable) instead of the gamma-law. Configured per state:
//
//   gas = {
//       type   = 'tabulated',
//       table  = 'data/D_fpeos.eostab',  -- path relative to the run dir
//       mass   = 1.0,                    -- scalar or per-tracer array,
//       charge = 0.0,                    --   TPG conventions (code units)
//       names  = {'D'},                  -- optional component names
//       ttol       = 1e-10,              -- inversion residual tolerance
//       max_newton = 100,                -- inversion iteration cap
//   }
//
// Composition semantics (v1): ONE thermodynamic material per state. The
// mass/charge arrays follow the ThermallyPerfectGas convention so the
// non-virtual base helpers (get_mass_from_cons etc., used by plasma5 /
// Lorentz / collisions) work unchanged with tracers present — but tracers
// are PASSIVE for thermodynamics: the single table closes the state
// regardless of the tracer fractions. Multi-material closure is the
// Stage-6 mixture driver.
//
// The `Gamma` primitive slot is filled with the energy-consistent effective
// gamma  γₑ = 1 + p/(ρ e)  so the existing Riemann solvers run unchanged
// (`effective_gamma` mode, plan D2): the face ENERGY they reconstruct from
// p/(γₑ−1) is exact; only their internal wave-speed estimates are
// approximate. The CFL time step always uses the true table sound speed
// (get_speed_from_*). The `SpHeat` slot carries cp from the general-EOS
// identity  cp = cv + (T/ρ²)(∂p/∂T)² / (∂p/∂ρ).
//
// NOT compatible with MHD states (which hard-code a constant gamma); a
// config-time guard aborts that combination (plan W14, lifting path in
// plan §6).

#include "MFP_eos_table.H"
#include "MFP_hydro_gas.H"

class TabulatedEOS : public HydroGas
{
  public:
    TabulatedEOS();
    TabulatedEOS(const int global_idx, const sol::table& def);
    ~TabulatedEOS();

    virtual Real get_gamma_from_cons(const Vector<Real>& U,
                                     const int density_idx = +HydroDef::ConsIdx::Density,
                                     const int tracer_idx = +HydroDef::ConsIdx::NUM) const override;
    virtual Real get_cp_from_cons(const Vector<Real>& U,
                                  const int density_idx = +HydroDef::ConsIdx::Density,
                                  const int tracer_idx = +HydroDef::ConsIdx::NUM) const override;

    virtual Real get_gamma_from_prim(const Vector<Real>& Q,
                                     const int idx = +HydroDef::PrimIdx::NUM) const override;
    virtual Real get_cp_from_prim(const Vector<Real>& Q,
                                  const int tracer_idx = +HydroDef::PrimIdx::NUM) const override;

    virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const override;
    virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const override;
    virtual void define_rho_p_T(Vector<Real>& Q) const override;

    virtual Real get_temperature_from_cons(const Vector<Real>& U) const override;
    virtual RealArray get_speed_from_cons(const Vector<Real>& U) const override;
    virtual RealArray get_speed_from_prim(const Vector<Real>& Q) const override;

    virtual std::string get_tag() const override { return tag; }
    virtual void write_info(nlohmann::json& js) const override;

    static bool registered;
    static std::string tag;

  private:
    // (rho, e_int) -> full evaluation via the re inversion; the workhorse
    // behind cons2prim and the from_cons getters
    void eval_from_rho_e(Real rho, Real e_int, EosEval& ev) const;
    // (rho, p) -> full evaluation via the rp inversion (from_prim getters)
    void eval_from_rho_p(Real rho, Real p, EosEval& ev) const;

    EosTable table;    // owner: loaded + nondimensionalised in the ctor
    EosTableView tv;   // cached non-owning view of `table`
};

#endif  // MFP_TABULATED_GAS_H
