#include "MFP_hydro_gas.H"

#include "MFP.H"
#include "MFP_state.H"
#include "MFP_hydro.H"

// ====================================================================================

HydroGas::HydroGas()
{
}

HydroGas::~HydroGas()
{
}

void HydroGas::get_alpha_fractions_from_prim(const Vector<Real> &Q, Vector<Real> &alpha, const int tracer_idx) const
{
    alpha.resize(n_species(), 0.0);

    std::copy(Q.begin()+tracer_idx, Q.begin()+tracer_idx+n_tracers(), alpha.begin());

    Real sum_alpha = 0.0;
    for (size_t i=0; i<n_tracers(); ++i) {
        sum_alpha += alpha[i];
    }

    alpha[n_species()-1] = 1.0 - sum_alpha;
}

void HydroGas::get_alpha_fractions_from_cons(const Vector<Real> &U,
                                               Vector<Real> &alpha,
                                               const int density_idx,
                                               const int tracer_idx) const
{
    alpha.resize(n_species(), 0.0);

    std::copy(U.begin()+tracer_idx, U.begin()+tracer_idx+n_tracers(), alpha.begin());

    Real rho = U[density_idx];
    Real sum_alpha = 0.0;
    for (size_t i=0; i<n_tracers(); ++i) {
        alpha[i] /= rho;
        sum_alpha += alpha[i];
    }

    alpha[n_species()-1] = 1.0 - sum_alpha;
}

Real HydroGas::get_mass_from_prim(const Vector<Real> &Q, const int tracer_idx) const
{
    BL_PROFILE("HydroGas::get_mass_from_prim");

    if (mass_const) return mass[0];

    Real S_alphai_mi = 0.0;
    Real S_alphas = 0.0;

    Real alphai = 0.0;
    for (int i = 0; i < n_tracers(); ++i) {
        alphai = Q[tracer_idx + i];
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphas += alphai;
        S_alphai_mi += alphai / mass[i];
    }

    S_alphai_mi += (1.0-S_alphas) / mass[n_tracers()];

    return 1.0 / S_alphai_mi;
}

Real HydroGas::get_mass_from_cons(const Vector<Real> &U, const int density_idx, const int tracer_idx) const
{
    BL_PROFILE("HydroGas::get_mass_from_cons");

    if (mass_const) return mass[0];

    Real rho = U[density_idx];

    Real S_alphai_mi = 0.0;
    Real S_alphas = 0.0;

    Real alphai = 0.0;
    for (int i = 0; i < n_tracers(); ++i) {
        alphai = U[tracer_idx + i] / rho;
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphas += alphai;
        S_alphai_mi += alphai / mass[i];
    }

    S_alphai_mi += (1.0-S_alphas) / mass[n_tracers()];

    return 1.0 / S_alphai_mi;
}

Real HydroGas::get_charge_from_prim(const Vector<Real> &Q, const int tracer_idx) const
{
    BL_PROFILE("HydroGas::get_charge_from_prim");

    if (charge_const) return charge[0];

    Real S_alphaiqi_mi = 0;
    Real S_alphai_mi = 0;
    Real S_alphas = 0;

    Real alphai = 0.0;
    for (int i = 0; i < n_tracers(); ++i) {
        alphai = Q[tracer_idx + i];
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphaiqi_mi += alphai * charge[i] / mass[i];
        S_alphai_mi += alphai / mass[i];
        S_alphas += alphai;
    }

    S_alphaiqi_mi += (1.0-S_alphas) * charge[n_tracers()] / mass[n_tracers()];
    S_alphai_mi += (1.0-S_alphas) / mass[n_tracers()];

    return S_alphaiqi_mi / S_alphai_mi;
}

Real HydroGas::get_charge_from_cons(const Vector<Real> &U, const int density_idx, const int tracer_idx) const
{
    BL_PROFILE("HydroGas::get_charge_from_cons");

    if (charge_const) return charge[0];

    Real rho = U[density_idx];

    Real S_alphaiqi_mi = 0;
    Real S_alphai_mi = 0;
    Real S_alphas = 0;

    Real alphai = 0.0;
    for (int i = 0; i < n_tracers(); ++i) {
        alphai = U[tracer_idx + i] / rho;
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphaiqi_mi += alphai * charge[i] / mass[i];
        S_alphai_mi += alphai / mass[i];
        S_alphas += alphai;
    }

    S_alphaiqi_mi += (1.0-S_alphas) * charge[n_tracers()] / mass[n_tracers()];
    S_alphai_mi += (1.0-S_alphas) / mass[n_tracers()];

    return S_alphaiqi_mi / S_alphai_mi;
}

bool HydroGas::prim_valid(const Vector<Real> &Q) const
{
    if ((Q[+HydroDef::PrimIdx::Density] <= 0.0) ||  (Q[+HydroDef::PrimIdx::Prs] <= 0.0)
            ||  (Q[+HydroDef::PrimIdx::Temp] <= 0.0)
            ) {
        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

bool HydroGas::cons_valid(const Vector<Real> &U) const
{
    if ((U[+HydroDef::ConsIdx::Density] <= 0.0) ||  (U[+HydroDef::ConsIdx::Eden] <= 0.0)
            ) {
        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

Real HydroGas::get_energy_from_cons(const Vector<Real> &U) const
{
    return U[+HydroDef::ConsIdx::Eden];
}

Real HydroGas::get_temperature_from_prim(const Vector<Real> &Q) const
{
    return Q[+HydroDef::PrimIdx::Temp];
}

Real HydroGas::get_internal_energy_density_from_cons(const Vector<Real>& U) const
{

    BL_PROFILE("ThermallyPerfectGas::get_specific_internal_energy_from_cons");

    Real rho = U[+HydroDef::ConsIdx::Density];
    Real mx = U[+HydroDef::ConsIdx::Xmom];
    Real my = U[+HydroDef::ConsIdx::Ymom];
    Real mz = U[+HydroDef::ConsIdx::Zmom];
    Real ed = U[+HydroDef::ConsIdx::Eden];

    Real rhoinv = 1/rho;
    Real u = mx*rhoinv;
    Real v = my*rhoinv;
    Real w = mz*rhoinv;
    Real ke = 0.5*rho*(u*u + v*v + w*w);

    Real ied = ed - ke; // internal energy density

    return ied;
}

Real HydroGas::get_density_from_number_density(const Real nd, Vector<Real>& Q) const
{

    Real mass_inv = 0.0;
    if (mass_const)
        return nd * mass[0];

    Real alphai = 0.0;
    Real S_alphas = 0.0;
    for (int i=0; i < n_tracers(); ++i) {
        alphai = Q[+HydroDef::PrimIdx::NUM + i];
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphas += alphai;
        mass_inv += alphai / mass[i];
    }
    mass_inv += (1.0-S_alphas) / mass[n_tracers()]; // deduce alpha value for final component

    return nd / mass_inv;
}

void HydroGas::write_info(nlohmann::json &js) const
{
    js["mass"] = mass;
    js["charge"] = charge;
    js["comp_names"] = comp_names;
}

ClassFactory<HydroGas>& GetHydroGasFactory() {
    static ClassFactory<HydroGas> F;
    return F;
}
