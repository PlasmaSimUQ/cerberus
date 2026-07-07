#include "MFP_tabulated_gas.H"

#include "MFP_lua.H"

#include <cmath>
#include <limits>

std::string TabulatedEOS::tag = "tabulated";
bool TabulatedEOS::registered =
  GetHydroGasFactory().Register(TabulatedEOS::tag, HydroGasBuilder<TabulatedEOS>);

TabulatedEOS::TabulatedEOS() {}
TabulatedEOS::~TabulatedEOS() {}

TabulatedEOS::TabulatedEOS(const int global_idx, const sol::table& def)
{
    BL_PROFILE("TabulatedEOS::TabulatedEOS");

    idx = global_idx;

    const std::string name = MFP::state_names[idx];

    // mass/charge follow the ThermallyPerfectGas array conventions so the
    // non-virtual base helpers work unchanged; thermodynamically the table
    // is the single closure (see header)
    set_values(def["mass"], mass);
    set_values(def["charge"], charge);
    set_values(def.get_or("names", sol::object()), comp_names);

    if (any_equal(mass.begin(), mass.end(), 0.0))
        Abort("State: " + name + "; mass cannot be 0");
    if (mass.size() != charge.size())
        Abort("State: " + name + "; 'mass' and 'charge' must have the same number of components");

    mass_const = all_equal(mass.begin(), mass.end(), mass[0]);
    charge_const = all_equal(charge.begin(), charge.end(), charge[0]);

    if (comp_names.empty()) {
        if (n_species() == 1) {
            comp_names.push_back(name);
        } else {
            for (int i = 0; i < n_species(); ++i) { comp_names.push_back(name + "_" + num2str(i)); }
        }
    }

    const std::string table_path = def["table"].get_or<std::string>("");
    if (table_path.empty())
        Abort("State: " + name + "; gas type 'tabulated' requires a 'table' file path");

    table.load(table_path);

    table.ttol = def["ttol"].get_or(1.0e-10);
    table.max_newton = def["max_newton"].get_or(100);

    // Nondimensionalise once at load (plan D7). Ordering is guaranteed:
    // MFP::update_ref() runs before any state's init_from_lua()
    // (MFP_config.cpp), but assert rather than assume. The table stores CGS
    // (frozen .eostab spec) while the reference quantities are SI, so the
    // refs are converted to CGS before dividing:
    //   rho: kg/m^3 -> g/cc (1e-3), u: m/s -> cm/s (1e2), p: Pa -> barye (10)
    if (!(MFP::u_ref > 0.0) || !(MFP::prs_ref > 0.0) || !(MFP::rho_ref > 0.0))
        Abort("State: " + name + "; reference quantities not set before tabulated-gas construction");
    table.nondimensionalise(MFP::rho_ref * 1.0e-3,
                            MFP::T_ref,
                            MFP::prs_ref * 10.0,
                            MFP::u_ref * 1.0e2);

    tv = table.view();

    amrex::Print() << "TabulatedEOS[" << name << "]: hull (code units) rho=["
                   << table.rho_min() << "," << table.rho_max() << "] T=["
                   << table.T_min() << "," << table.T_max() << "]\n";
}

// ---------------------------------------------------------------------------
// internal drivers

void TabulatedEOS::eval_from_rho_e(Real rho, Real e_int, EosEval& ev) const
{
    EosInvertStats st;
    const Real T =
      EosTable::invert_T_from_e(tv, rho, e_int, -1.0, table.ttol, table.max_newton, st);
    EosTable::eval_rt(tv, rho, T, ev);
}

void TabulatedEOS::eval_from_rho_p(Real rho, Real p, EosEval& ev) const
{
    EosInvertStats st;
    const Real T =
      EosTable::invert_T_from_p(tv, rho, p, -1.0, table.ttol, table.max_newton, st);
    EosTable::eval_rt(tv, rho, T, ev);
}

// ---------------------------------------------------------------------------
// conversions

bool TabulatedEOS::cons2prim(Vector<Real>& U, Vector<Real>& Q) const
{
    BL_PROFILE("TabulatedEOS::cons2prim");

    Real rho = U[+HydroDef::ConsIdx::Density];
    Real mx = U[+HydroDef::ConsIdx::Xmom];
    Real my = U[+HydroDef::ConsIdx::Ymom];
    Real mz = U[+HydroDef::ConsIdx::Zmom];
    Real ed = U[+HydroDef::ConsIdx::Eden];

    Real rhoinv = 1 / rho;
    Real u = mx * rhoinv;
    Real v = my * rhoinv;
    Real w = mz * rhoinv;
    Real e_int = ed * rhoinv - 0.5 * (u * u + v * v + w * w);

    // re inversion (inverse-map seed + Newton polish, plan D5/D8); an
    // unattainable e clamps to the hull edge inside the driver (flagged)
    EosEval ev;
    eval_from_rho_e(rho, e_int, ev);

    // energy-consistent effective gamma: what the (unchanged) Riemann
    // solvers reconstruct face energy from. Use the physical e_int so
    // p/(gamma_e - 1) = rho*e_int is exact where the inversion converged.
    const Real ge = 1.0 + ev.p / std::max(rho * e_int, std::numeric_limits<Real>::min());

    // general-EOS specific heat identity (reduces to gamma-law cp on the
    // synthetic table)
    const Real cp = ev.cv + (ev.T / (rho * rho)) * ev.dpdT * ev.dpdT /
                              std::max(ev.dpdrho, std::numeric_limits<Real>::min());

    Q[+HydroDef::PrimIdx::Density] = rho;
    Q[+HydroDef::PrimIdx::Xvel] = u;
    Q[+HydroDef::PrimIdx::Yvel] = v;
    Q[+HydroDef::PrimIdx::Zvel] = w;
    Q[+HydroDef::PrimIdx::Prs] = ev.p;
    Q[+HydroDef::PrimIdx::Temp] = ev.T;
    Q[+HydroDef::PrimIdx::Gamma] = ge;
    Q[+HydroDef::PrimIdx::SpHeat] = cp;

    for (int i = 0; i < n_tracers(); ++i) {
        Q[+HydroDef::PrimIdx::NUM + i] = U[+HydroDef::ConsIdx::NUM + i] * rhoinv;
    }

#ifdef MFP_PRIM_FLOOR
    // same floor policy as ThermallyPerfectGas::cons2prim
    if (Q[+HydroDef::PrimIdx::Density] < effective_zero) {
        Q[+HydroDef::PrimIdx::Density] = effective_zero;
    }
    if (Q[+HydroDef::PrimIdx::Prs] < effective_zero) {
        Q[+HydroDef::PrimIdx::Prs] = effective_zero;
    }
    if (Q[+HydroDef::PrimIdx::Temp] < effective_zero) {
        Q[+HydroDef::PrimIdx::Temp] = effective_zero;
    }
#endif

    return prim_valid(Q);
}

void TabulatedEOS::prim2cons(Vector<Real>& Q, Vector<Real>& U) const
{
    BL_PROFILE("TabulatedEOS::prim2cons");

    Real rho = Q[+HydroDef::PrimIdx::Density];
    Real u = Q[+HydroDef::PrimIdx::Xvel];
    Real v = Q[+HydroDef::PrimIdx::Yvel];
    Real w = Q[+HydroDef::PrimIdx::Zvel];
    Real p = Q[+HydroDef::PrimIdx::Prs];

    Real mx = u * rho;
    Real my = v * rho;
    Real mz = w * rho;
    Real ke = 0.5 * rho * (u * u + v * v + w * w);

    // rp inversion -> T -> e (the Temp slot is not trusted here; the seed
    // map makes the guess irrelevant anyway)
    EosEval ev;
    eval_from_rho_p(rho, p, ev);

    U[+HydroDef::ConsIdx::Density] = rho;
    U[+HydroDef::ConsIdx::Xmom] = mx;
    U[+HydroDef::ConsIdx::Ymom] = my;
    U[+HydroDef::ConsIdx::Zmom] = mz;
    U[+HydroDef::ConsIdx::Eden] = rho * ev.e + ke;

    for (int i = 0; i < n_tracers(); ++i) {
        U[+HydroDef::ConsIdx::NUM + i] = Q[+HydroDef::PrimIdx::NUM + i] * rho;
    }
}

void TabulatedEOS::define_rho_p_T(Vector<Real>& Q) const
{
    BL_PROFILE("TabulatedEOS::define_rho_p_T");

    Real rho = Q[+HydroDef::PrimIdx::Density];
    Real p = Q[+HydroDef::PrimIdx::Prs];
    Real T = Q[+HydroDef::PrimIdx::Temp];

    // same given-ness convention as ThermallyPerfectGas: positive = given,
    // same priority order
    EosInvertStats st;
    if ((rho > 0.0) && (p > 0.0)) {
        T = EosTable::invert_T_from_p(tv, rho, p, T, table.ttol, table.max_newton, st);
    } else if ((p > 0.0) && (T > 0.0)) {
        rho = EosTable::invert_rho_from_p(tv, T, p, rho, table.ttol, table.max_newton, st);
    } else if ((rho > 0.0) && (T > 0.0)) {
        EosEval ev;
        EosTable::eval_rt(tv, rho, T, ev);
        p = ev.p;
    }

    Q[+HydroDef::PrimIdx::Density] = rho;
    Q[+HydroDef::PrimIdx::Prs] = p;
    Q[+HydroDef::PrimIdx::Temp] = T;
}

// ---------------------------------------------------------------------------
// getters

Real TabulatedEOS::get_temperature_from_cons(const Vector<Real>& U) const
{
    BL_PROFILE("TabulatedEOS::get_temperature_from_cons");

    const Real rho = U[+HydroDef::ConsIdx::Density];
    const Real rhoinv = 1 / rho;
    const Real u = U[+HydroDef::ConsIdx::Xmom] * rhoinv;
    const Real v = U[+HydroDef::ConsIdx::Ymom] * rhoinv;
    const Real w = U[+HydroDef::ConsIdx::Zmom] * rhoinv;
    const Real e_int = U[+HydroDef::ConsIdx::Eden] * rhoinv - 0.5 * (u * u + v * v + w * w);

    EosInvertStats st;
    return EosTable::invert_T_from_e(tv, rho, e_int, -1.0, table.ttol, table.max_newton, st);
}

Real TabulatedEOS::get_gamma_from_cons(const Vector<Real>& U,
                                       const int density_idx,
                                       const int tracer_idx) const
{
    BL_PROFILE("TabulatedEOS::get_gamma_from_cons");
    amrex::ignore_unused(tracer_idx);

    const Real rho = U[density_idx];
    const Real rhoinv = 1 / rho;
    const Real u = U[+HydroDef::ConsIdx::Xmom] * rhoinv;
    const Real v = U[+HydroDef::ConsIdx::Ymom] * rhoinv;
    const Real w = U[+HydroDef::ConsIdx::Zmom] * rhoinv;
    const Real e_int = U[+HydroDef::ConsIdx::Eden] * rhoinv - 0.5 * (u * u + v * v + w * w);

    EosEval ev;
    eval_from_rho_e(rho, e_int, ev);
    return 1.0 + ev.p / std::max(rho * e_int, std::numeric_limits<Real>::min());
}

Real TabulatedEOS::get_gamma_from_prim(const Vector<Real>& Q, const int idx) const
{
    BL_PROFILE("TabulatedEOS::get_gamma_from_prim");
    amrex::ignore_unused(idx);

    const Real rho = Q[+HydroDef::PrimIdx::Density];
    const Real p = Q[+HydroDef::PrimIdx::Prs];

    EosEval ev;
    eval_from_rho_p(rho, p, ev);
    return 1.0 + p / std::max(rho * ev.e, std::numeric_limits<Real>::min());
}

Real TabulatedEOS::get_cp_from_cons(const Vector<Real>& U,
                                    const int density_idx,
                                    const int tracer_idx) const
{
    BL_PROFILE("TabulatedEOS::get_cp_from_cons");
    amrex::ignore_unused(tracer_idx);

    const Real rho = U[density_idx];
    const Real rhoinv = 1 / rho;
    const Real u = U[+HydroDef::ConsIdx::Xmom] * rhoinv;
    const Real v = U[+HydroDef::ConsIdx::Ymom] * rhoinv;
    const Real w = U[+HydroDef::ConsIdx::Zmom] * rhoinv;
    const Real e_int = U[+HydroDef::ConsIdx::Eden] * rhoinv - 0.5 * (u * u + v * v + w * w);

    EosEval ev;
    eval_from_rho_e(rho, e_int, ev);
    return ev.cv + (ev.T / (rho * rho)) * ev.dpdT * ev.dpdT /
                     std::max(ev.dpdrho, std::numeric_limits<Real>::min());
}

Real TabulatedEOS::get_cp_from_prim(const Vector<Real>& Q, const int tracer_idx) const
{
    BL_PROFILE("TabulatedEOS::get_cp_from_prim");
    amrex::ignore_unused(tracer_idx);

    const Real rho = Q[+HydroDef::PrimIdx::Density];
    const Real p = Q[+HydroDef::PrimIdx::Prs];

    EosEval ev;
    eval_from_rho_p(rho, p, ev);
    return ev.cv + (ev.T / (rho * rho)) * ev.dpdT * ev.dpdT /
                     std::max(ev.dpdrho, std::numeric_limits<Real>::min());
}

RealArray TabulatedEOS::get_speed_from_cons(const Vector<Real>& U) const
{
    BL_PROFILE("TabulatedEOS::get_speed_from_cons");

    Real rho = U[+HydroDef::ConsIdx::Density];

#ifdef MFP_PRIM_FLOOR
    // raw conserved data that has not passed through the cons2prim floors
    rho = std::max(rho, effective_zero);
#endif

    const Real rhoinv = 1 / rho;
    const Real u = U[+HydroDef::ConsIdx::Xmom] * rhoinv;
    const Real v = U[+HydroDef::ConsIdx::Ymom] * rhoinv;
    const Real w = U[+HydroDef::ConsIdx::Zmom] * rhoinv;
    Real e_int = U[+HydroDef::ConsIdx::Eden] * rhoinv - 0.5 * (u * u + v * v + w * w);

#ifdef MFP_PRIM_FLOOR
    e_int = std::max(e_int, effective_zero);
#endif

    // true table sound speed -> the CFL time step is exact even in
    // effective_gamma flux mode (plan D2)
    EosEval ev;
    eval_from_rho_e(rho, e_int, ev);
    const Real a = ev.cs;

    RealArray s = {AMREX_D_DECL(a + std::abs(u), a + std::abs(v), a + std::abs(w))};

    return s;
}

RealArray TabulatedEOS::get_speed_from_prim(const Vector<Real>& Q) const
{
    BL_PROFILE("TabulatedEOS::get_speed_from_prim");

#ifdef MFP_PRIM_FLOOR
    const Real rho = std::max(Q[+HydroDef::PrimIdx::Density], effective_zero);
    const Real p = std::max(Q[+HydroDef::PrimIdx::Prs], effective_zero);
#else
    const Real rho = Q[+HydroDef::PrimIdx::Density];
    const Real p = Q[+HydroDef::PrimIdx::Prs];
#endif

    EosEval ev;
    eval_from_rho_p(rho, p, ev);
    const Real a = ev.cs;

    RealArray s = {AMREX_D_DECL(a + std::abs(Q[+HydroDef::PrimIdx::Xvel]),
                                a + std::abs(Q[+HydroDef::PrimIdx::Yvel]),
                                a + std::abs(Q[+HydroDef::PrimIdx::Zvel]))};

    return s;
}

void TabulatedEOS::write_info(nlohmann::json& js) const
{
    BL_PROFILE("TabulatedEOS::write_info");

    HydroGas::write_info(js);

    js["type"] = tag;
    for (const auto& kv : table.provenance()) { js["table_" + kv.first] = kv.second; }
}
