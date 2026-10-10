#ifndef COLLISIONS_H
#define COLLISIONS_H

#include "MFP_action.H"
#include "MFP_hydro.H"

// Per-cell multi-species reactor base. Lua `species` names resolve to the
// alpha sub-components (gas.names) of the hydro `states`; calc_update() fills
// species_info[*].delta (d/dt of rho, mom, E per sub-species, code units) and
// the write-back applies them to the state AND its tracers. Ordering contract
// relied on by the tracer write-back: sub-species i == tracer i for
// i < n_tracers (gas.mass/charge/names order), the last entry is the
// implicit 1 - sum(alpha) component (HydroGas::get_alpha_fractions_from_prim).
// Sources are evaluated several times per step (Strang/Heun), so
// calc_update() must be a pure function of the gathered state: no
// side-effect accumulators.
class Collisions : public Action
{
  public:
    Collisions();
    Collisions(const int idx, const sol::table& def);
    ~Collisions();

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    virtual void calc_time_derivative(MFP* mfp,
                                      Vector<UpdateData>& update,
                                      const Real time,
                                      const Real dt) override;

    struct SpeciesInfo {
        std::string name;
        int idx;
        int state_idx;
        int alpha_idx;
        Real rho, m, q, T, q2, n;
        Array<Real, 3> vel;
        Array<Real, +HydroDef::ConsIdx::NUM> delta;
    };

    void init_species_info(const sol::table& def);

    virtual void calc_update() {}

    // index into species_info of the (first) q < 0 species, or -1
    int electron_species() const;

    // lagged, opt-in limit: dt_fraction * min(tau) over every evaluation
    // since the previous query (Lua `dt_fraction`, default 0 = off)
    virtual Real get_allowed_time_step(MFP* mfp) override;

    Vector<HydroState*> states;
    Vector<SpeciesInfo> species_info;

    Real dt_stage = 0.0;  // dt of the current evaluation (for positivity clamps)
    Real tau_min = std::numeric_limits<Real>::max();  // set by calc_update()
    Real dt_fraction = 0.0;
};

#endif  // LORENTZ_H
