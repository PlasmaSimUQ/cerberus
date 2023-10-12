#include "MFP_thermally_perfect_gas.H"
#include "MFP_lua.H"

std::string ThermallyPerfectGas::tag = "thermally_perfect";
bool ThermallyPerfectGas::registered = GetHydroGasFactory().Register(ThermallyPerfectGas::tag, HydroGasBuilder<ThermallyPerfectGas>);


ThermallyPerfectGas::ThermallyPerfectGas(){}
ThermallyPerfectGas::~ThermallyPerfectGas(){}

ThermallyPerfectGas::ThermallyPerfectGas(const int global_idx, const sol::table &def)
{
    idx = global_idx;

    std::string name = MFP::state_names[idx];

    //
    // get mass, charge, gamma
    //

    set_values(def["mass"], mass);
    set_values(def["charge"], charge);
    set_values(def["gamma"], gamma);
    set_values(def.get_or("names", sol::object()), comp_names);

    if (any_equal(mass.begin(), mass.end(), 0.0) or any_equal(gamma.begin(), gamma.end(), 0.0))
        Abort("State: "+name+"; mass and gamma cannot be 0");

    mass_const = all_equal(mass.begin(), mass.end(), mass[0]);
    charge_const = all_equal(charge.begin(), charge.end(), charge[0]);
    gamma_const = all_equal(gamma.begin(), gamma.end(), gamma[0]);

    if ((mass.size() != charge.size()) or (mass.size() != gamma.size()) or (charge.size() != gamma.size()))
        Abort("State: "+name+"; 'mass', 'charge' and 'gamma' must have the same number of components");

    // handle the names of the sub components if they haven't been provided
    if (comp_names.empty()) {
        if (n_species() == 1) {
            comp_names.push_back(name);
        } else {
            for (int i=0; i<n_species(); ++i) {
                comp_names.push_back(name+"_"+num2str(i));
            }
        }
    }

}

Real ThermallyPerfectGas::get_gamma_from_prim(const Vector<Real> &Q, const int idx) const
{
    BL_PROFILE("ThermallyPerfectGas::get_gamma_from_prim");

    if (gamma_const) return gamma[0];

    Real S_alphaicpi = 0.0;
    Real S_alphaicvi = 0.0;
    Real S_alphai = 0.0;

    Real cvi = 0.0;
    Real cpi = 0.0;
    Real alphai = 0.0;
    for (int i = 0; i < n_tracers(); ++i) {
        alphai = Q[idx + i];
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphai += alphai;

        cpi = gamma[i] / (mass[i] * (gamma[i] - 1.0));
        cvi = 1.0 / (mass[i] * (gamma[i] - 1.0));

        S_alphaicpi += alphai * cpi;
        S_alphaicvi += alphai * cvi;
    }

    cpi = gamma[n_tracers()] / (mass[n_tracers()] * (gamma[n_tracers()] - 1.0));
    cvi = 1.0 / (mass[n_tracers()] * (gamma[n_tracers()] - 1.0));

    S_alphaicpi += (1.0 - S_alphai) * cpi;
    S_alphaicvi += (1.0 - S_alphai) * cvi;

    return S_alphaicpi / S_alphaicvi;
}

Real ThermallyPerfectGas::get_gamma_from_cons(const Vector<Real> &U, const int density_idx, const int tracer_idx) const
{
    BL_PROFILE("ThermallyPerfectGas::get_gamma_from_cons");

    if (gamma_const) return gamma[0];

    Real rho = U[density_idx];

    Real S_alphaicpi = 0.0;
    Real S_alphaicvi = 0.0;
    Real S_alphai = 0.0;

    Real cvi = 0.0;
    Real cpi = 0.0;
    Real alphai = 0.0;
    for (int i = 0; i < n_tracers(); ++i) {
        alphai = U[tracer_idx + i] / rho;
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphai += alphai;

        cpi = gamma[i] / (mass[i] * (gamma[i] - 1.0));
        cvi = 1.0 / (mass[i] * (gamma[i] - 1.0));

        S_alphaicpi += alphai * cpi;
        S_alphaicvi += alphai * cvi;
    }

    cpi = gamma[n_tracers()] / (mass[n_tracers()] * (gamma[n_tracers()] - 1.0));
    cvi = 1.0 / (mass[n_tracers()] * (gamma[n_tracers()] - 1.0));

    S_alphaicpi += (1.0 - S_alphai) * cpi;
    S_alphaicvi += (1.0 - S_alphai) * cvi;

    return S_alphaicpi / S_alphaicvi;
}

Real ThermallyPerfectGas::get_cp_from_prim(const Vector<Real> &Q, const int tracer_idx) const
{
    BL_PROFILE("ThermallyPerfectGas::get_cp_from_prim");

    if (gamma_const && mass_const) return gamma[0]/(mass[0]*(gamma[0]-1));

    Real S_alphai = 0.0;
    Real S_alphaicpi = 0.0;

    Real cpi = 0.0;
    Real alphai = 0.0;
    for (int i = 0; i < n_tracers(); ++i) {
        alphai = Q[tracer_idx + i];
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphai += alphai;

        cpi = gamma[i] / (mass[i] * (gamma[i] - 1.0));
        S_alphaicpi += alphai * cpi;
    }

    cpi = gamma[n_tracers()] / (mass[n_tracers()] * (gamma[n_tracers()] - 1.0));
    S_alphaicpi += (1.0- S_alphai) * cpi;

    return S_alphaicpi;
}

Real ThermallyPerfectGas::get_cp_from_cons(const Vector<Real> &U, const int density_idx, const int tracer_idx) const
{
    BL_PROFILE("ThermallyPerfectGas::get_cp_from_cons");

    if (gamma_const && mass_const) return gamma[0]/(mass[0]*(gamma[0]-1));

    Real rho = U[density_idx];

    Real S_alphai = 0.0;
    Real S_alphaicpi = 0.0;

    Real cpi = 0.0;
    Real alphai = 0.0;
    for (int i = 0; i < n_tracers(); ++i) {
        alphai = U[tracer_idx + i] / rho;
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphai += alphai;

        cpi = gamma[i] / (mass[i] * (gamma[i] - 1.0));
        S_alphaicpi += alphai * cpi;
    }

    cpi = gamma[n_tracers()] / (mass[n_tracers()] * (gamma[n_tracers()] - 1.0));
    S_alphaicpi += (1.0- S_alphai) * cpi;

    return S_alphaicpi;

}


// in place conversion from conserved to primitive
bool ThermallyPerfectGas::cons2prim(Vector<Real>& U, Vector<Real>& Q) const
{
    BL_PROFILE("ThermallyPerfectGas::cons2prim");

    Real rho = U[+HydroDef::ConsIdx::Density];
    Real mx = U[+HydroDef::ConsIdx::Xmom];
    Real my = U[+HydroDef::ConsIdx::Ymom];
    Real mz = U[+HydroDef::ConsIdx::Zmom];
    Real ed = U[+HydroDef::ConsIdx::Eden];
  
    //Print() << "\n" << rho << "\t" << mx <<  "\t" << my << "\t" << mz << ed << "\t" ;  //TODO delete 

    Real rhoinv = 1/rho;
    Real u = mx*rhoinv;
    Real v = my*rhoinv;
    Real w = mz*rhoinv;
    Real ke = 0.5*rho*(u*u + v*v + w*w);
    Real m = get_mass_from_cons(U);
    Real g = get_gamma_from_cons(U);
    Real cp = get_cp_from_cons(U);
    Real p = (ed - ke)*(g - 1);
    Real T = p*rhoinv*m;

    Q[+HydroDef::PrimIdx::Density] = rho;
    Q[+HydroDef::PrimIdx::Xvel] = u;
    Q[+HydroDef::PrimIdx::Yvel] = v;
    Q[+HydroDef::PrimIdx::Zvel] = w;
    Q[+HydroDef::PrimIdx::Prs] = p;
    Q[+HydroDef::PrimIdx::Temp] = T;
    Q[+HydroDef::PrimIdx::Gamma] = g;
    Q[+HydroDef::PrimIdx::SpHeat] = cp;

    for (int i = 0; i < n_tracers(); ++i) {
        Q[+HydroDef::PrimIdx::NUM + i] = U[+HydroDef::ConsIdx::NUM + i] * rhoinv;
    }

    return prim_valid(Q);
}

void ThermallyPerfectGas::prim2cons(Vector<Real>& Q, Vector<Real>& U) const
{
    BL_PROFILE("ThermallyPerfectGas::prim2cons");

    Real rho = Q[+HydroDef::PrimIdx::Density];
    Real u = Q[+HydroDef::PrimIdx::Xvel];
    Real v = Q[+HydroDef::PrimIdx::Yvel];
    Real w = Q[+HydroDef::PrimIdx::Zvel];
    Real p = Q[+HydroDef::PrimIdx::Prs];

    Real mx = u*rho;
    Real my = v*rho;
    Real mz = w*rho;
    Real ke = 0.5*rho*(u*u + v*v + w*w);
    Real g = get_gamma_from_prim(Q);
    Real ed = p/(g - 1) + ke;

    U[+HydroDef::ConsIdx::Density] = rho;
    U[+HydroDef::ConsIdx::Xmom] = mx;
    U[+HydroDef::ConsIdx::Ymom] = my;
    U[+HydroDef::ConsIdx::Zmom] = mz;
    U[+HydroDef::ConsIdx::Eden] = ed;

    for (int i = 0; i < n_tracers(); ++i) {
        U[+HydroDef::ConsIdx::NUM + i] = Q[+HydroDef::PrimIdx::NUM + i] * rho;
    }

}

void ThermallyPerfectGas::define_rho_p_T(Vector<Real>& Q) const
{
    BL_PROFILE("ThermallyPerfectGas::prim2cons");

    Real rho = Q[+HydroDef::PrimIdx::Density];
    Real p = Q[+HydroDef::PrimIdx::Prs];
    Real T = Q[+HydroDef::PrimIdx::Temp];

    Real m = get_mass_from_prim(Q);

    if ((rho > 0.0) && (p > 0.0)) {
        T = p*m/rho;
    } else if ((p > 0.0) && (T > 0.0)) {
        rho = p*m/T;
    } else if ((rho > 0.0) && (T > 0.0)) {
        p = rho*T/m;
    }

    Q[+HydroDef::PrimIdx::Density] = rho;
    Q[+HydroDef::PrimIdx::Prs] = p;
    Q[+HydroDef::PrimIdx::Temp] = T;

}

Real ThermallyPerfectGas::get_temperature_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("ThermallyPerfectGas::get_temperature_from_cons");

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
    Real m = get_mass_from_cons(U);
    Real g = get_gamma_from_cons(U);
    Real p = (ed - ke)*(g - 1);
    Real T = p*rhoinv*m;

    return T;

}

RealArray ThermallyPerfectGas::get_speed_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("ThermallyPerfectGas::get_speed_from_cons");
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
    Real g = get_gamma_from_cons(U);
    Real p = (ed - ke)*(g - 1);

    Real a = std::sqrt(g*p*rhoinv);

    RealArray s = {AMREX_D_DECL(a + std::abs(u), a + std::abs(v), a + std::abs(w))};

    return s;

}

RealArray ThermallyPerfectGas::get_speed_from_prim(const Vector<Real> &Q) const
{
    BL_PROFILE("ThermallyPerfectGas::get_speed_from_prim");

    Real g = get_gamma_from_prim(Q);

    Real a = std::sqrt(g*Q[+HydroDef::PrimIdx::Prs]/Q[+HydroDef::PrimIdx::Density]);

    RealArray s = {AMREX_D_DECL(a + std::abs(Q[+HydroDef::PrimIdx::Xvel]),
                   a + std::abs(Q[+HydroDef::PrimIdx::Yvel]),
                   a + std::abs(Q[+HydroDef::PrimIdx::Zvel]))};


    return s;

}

void ThermallyPerfectGas::write_info(nlohmann::json &js) const
{

    HydroGas::write_info(js);

    js["type"] = tag;
    js["gamma"] = gamma;

}
