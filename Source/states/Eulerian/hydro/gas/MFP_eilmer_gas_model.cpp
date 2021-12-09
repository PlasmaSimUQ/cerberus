#ifdef EILMER_GAS
#include "MFP_eilmer_gas_model.H"
#include "MFP_lua.H"
#include "MFP_wrap_gas.H"

std::string EilmerGasModel::tag = "eilmer";
bool EilmerGasModel::registered = GetHydroGasFactory().Register(EilmerGasModel::tag, HydroGasBuilder<EilmerGasModel>);


EilmerGasModel::EilmerGasModel(){}
EilmerGasModel::~EilmerGasModel(){}

EilmerGasModel::EilmerGasModel(const int global_idx, const sol::table &def)
{
    idx = global_idx;

    std::string name = MFP::state_names[idx];

    // grab stuff from the config
    set_values(def.get_or("names", sol::object()), comp_names);
    set_values(def["charge"], charge);
    mass.resize(comp_names.size(), 0.0);

    if (mass.size() != charge.size())
        Abort("State: "+name+"; 'charge' and 'names' must have the same number of entries");

    // set up the Eilmer gas model

    EilmerGas::initialise();

    std::string gas_model_file = def.get_or<std::string>("gas_model","");
    if (gas_model_file.empty()) Abort("Action '"+name+"' requires 'gas_model' to be defined (a lua file)");

    gas_model_id = gas_model_new(gas_model_file.data());

    if (gas_model_id < 0) Abort("State '"+name+"' has failed when trying to create a new gas model");

    n_species = gas_model_n_species(gas_model_id);

    // get the molecular masses
    Vector<Real> mm(n_species);
    gas_model_mol_masses(gas_model_id, mm.data());

    // get the mapping between the MFP gas state and the Eilmer gas model

    species_info.resize(n_species);

    std::pair<bool, int > found;
    int n;
    for (int i=0; i<n_species; ++i) {
        std::string& sname = species_info[i].name;
        sname.resize(10);
        gas_model_species_name_and_length(gas_model_id, i, sname.data(), &n);
        sname.resize(n);

        found = findInVector(comp_names, sname);
        if (found.first) {
            species_info[i].alpha_idx = found.second;
            // update the mass to match that given by the gas model
            mass[found.second] = (mm[i]/6.02214076e23)/MFP::m_ref;
        } else {
            species_info[i].alpha_idx = -1;
        }
    }

    gas_state_id = gas_state_new(gas_model_id);

    if (gas_state_id < 0) Abort("Action '"+name+"' has failed when trying to create a new gas state");



    if (any_equal(mass.begin(), mass.end(), 0.0)) {
        Abort("State "+name+" has species masses of zero: "+vec2str(mass));
    }

    mass_const = all_equal(mass.begin(), mass.end(), mass[0]);
    charge_const = all_equal(charge.begin(), charge.end(), charge[0]);

}

void EilmerGasModel::set_eilmer_gas_state_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("EilmerGasModel::set_eilmer_gas_state_from_cons");

    Real rho = U[+HydroDef::ConsIdx::Density];

    // get the mass fractions
    Vector<Real> alpha;
    get_alpha_fractions_from_cons(U, alpha);

    Vector<Real> massf(n_species, 0.0);
    for (int i=0; i<n_species; ++i) {
        const int alpha_idx = species_info[i].alpha_idx;
        if (alpha_idx > -1) {
            massf[i] = alpha[alpha_idx];
        }
    }


    Real specific_internal_nrg = get_internal_energy_density_from_cons(U)/rho;

    gas_state_set_scalar_field(gas_model_id, "rho", rho*MFP::rho_ref);
    gas_state_set_scalar_field(gas_model_id, "u", specific_internal_nrg*MFP::prs_ref/MFP::rho_ref);
    gas_state_set_array_field(gas_state_id, "massf", massf.data(), n_species);

    // need an initial guess at the temperature
//    gas_state_set_scalar_field(gas_model_id, "T", guess_T*MFP::T_ref);

    // update model and solve for update
    gas_model_gas_state_update_thermo_from_rhou(gas_model_id, gas_state_id);
}

void EilmerGasModel::set_eilmer_gas_state_from_prim(const Vector<Real> &Q) const
{
    BL_PROFILE("EilmerGasModel::set_eilmer_gas_state_from_prim");

    Real rho = Q[+HydroDef::PrimIdx::Density];
    Real p = Q[+HydroDef::PrimIdx::Prs];

    // get the mass fractions
    Vector<Real> alpha;
    get_alpha_fractions_from_prim(Q, alpha);

    Vector<Real> massf(n_species, 0.0);
    for (int i=0; i<n_species; ++i) {
        const int alpha_idx = species_info[i].alpha_idx;
        if (alpha_idx > -1) {
            massf[i] = alpha[alpha_idx];
        }
    }

    gas_state_set_scalar_field(gas_model_id, "rho", rho*MFP::rho_ref);
    gas_state_set_scalar_field(gas_model_id, "p", p*MFP::prs_ref);
    gas_state_set_array_field(gas_state_id, "massf", massf.data(), n_species);

    // update model and solve for update
    gas_model_gas_state_update_thermo_from_rhop(gas_model_id, gas_state_id);
}

Real EilmerGasModel::get_gamma_from_prim(const Vector<Real> &Q, const int idx) const
{
    BL_PROFILE("EilmerGasModel::get_gamma_from_prim");

    Real gamma;

    set_eilmer_gas_state_from_prim(Q);
    gas_model_gas_state_gamma(gas_model_id, gas_state_id, &gamma);

    return gamma;

}

Real EilmerGasModel::get_gamma_from_cons(const Vector<Real> &U, const int density_idx, const int tracer_idx) const
{
    BL_PROFILE("EilmerGasModel::get_gamma_from_cons");

    Real gamma;

    set_eilmer_gas_state_from_cons(U);
    gas_model_gas_state_gamma(gas_model_id, gas_state_id, &gamma);

    return gamma;
}

Real EilmerGasModel::get_cp_from_prim(const Vector<Real> &Q, const int tracer_idx) const
{
    BL_PROFILE("EilmerGasModel::get_cp_from_prim");

    Real cp;

    set_eilmer_gas_state_from_prim(Q);
    gas_model_gas_state_Cp(gas_model_id, gas_state_id, &cp);

    return cp;
}

Real EilmerGasModel::get_cp_from_cons(const Vector<Real> &U, const int density_idx, const int tracer_idx) const
{
    BL_PROFILE("EilmerGasModel::get_cp_from_cons");

    Real cp;

    set_eilmer_gas_state_from_cons(U);
    gas_model_gas_state_Cp(gas_model_id, gas_state_id, &cp);

    return cp;

}


// conversion from conserved to primitive
bool EilmerGasModel::cons2prim(Vector<Real>& U, Vector<Real>& Q) const
{
    BL_PROFILE("EilmerGasModel::cons2prim");

    set_eilmer_gas_state_from_cons(U);

    Real rho = U[+HydroDef::ConsIdx::Density];
    Real mx = U[+HydroDef::ConsIdx::Xmom];
    Real my = U[+HydroDef::ConsIdx::Ymom];
    Real mz = U[+HydroDef::ConsIdx::Zmom];

    Real rhoinv = 1/rho;
    Real u = mx*rhoinv;
    Real v = my*rhoinv;
    Real w = mz*rhoinv;

    Real g, cp, p, T;

    gas_model_gas_state_gamma(gas_model_id, gas_state_id, &g);
    gas_model_gas_state_Cp(gas_model_id, gas_state_id, &cp);
    gas_state_get_scalar_field(gas_state_id, "p", &p);
    gas_state_get_scalar_field(gas_state_id, "T", &T);

    Q[+HydroDef::PrimIdx::Density] = rho;
    Q[+HydroDef::PrimIdx::Xvel] = u;
    Q[+HydroDef::PrimIdx::Yvel] = v;
    Q[+HydroDef::PrimIdx::Zvel] = w;
    Q[+HydroDef::PrimIdx::Prs] = p/MFP::prs_ref;
    Q[+HydroDef::PrimIdx::Temp] = T/MFP::T_ref;
    Q[+HydroDef::PrimIdx::Gamma] = g;
    Q[+HydroDef::PrimIdx::SpHeat] = cp/(MFP::u_ref*MFP::u_ref/MFP::T_ref);

    for (int i = 0; i < n_tracers(); ++i) {
        Q[+HydroDef::PrimIdx::NUM + i] = U[+HydroDef::ConsIdx::NUM + i] * rhoinv;
    }

    return prim_valid(Q);
}

void EilmerGasModel::prim2cons(Vector<Real>& Q, Vector<Real>& U) const
{
    BL_PROFILE("EilmerGasModel::prim2cons");

    set_eilmer_gas_state_from_prim(Q);

    Real rho = Q[+HydroDef::PrimIdx::Density];
    Real u = Q[+HydroDef::PrimIdx::Xvel];
    Real v = Q[+HydroDef::PrimIdx::Yvel];
    Real w = Q[+HydroDef::PrimIdx::Zvel];

    Real mx = u*rho;
    Real my = v*rho;
    Real mz = w*rho;
    Real ke = 0.5*rho*(u*u + v*v + w*w);

    Real specific_internal_nrg;
    gas_state_get_scalar_field(gas_state_id, "u", &specific_internal_nrg);
    specific_internal_nrg /= MFP::prs_ref/MFP::rho_ref;

    Real ed = rho*specific_internal_nrg + ke;

    U[+HydroDef::ConsIdx::Density] = rho;
    U[+HydroDef::ConsIdx::Xmom] = mx;
    U[+HydroDef::ConsIdx::Ymom] = my;
    U[+HydroDef::ConsIdx::Zmom] = mz;
    U[+HydroDef::ConsIdx::Eden] = ed;

    for (int i = 0; i < n_tracers(); ++i) {
        U[+HydroDef::ConsIdx::NUM + i] = Q[+HydroDef::PrimIdx::NUM + i] * rho;
    }

}

void EilmerGasModel::define_rho_p_T(Vector<Real>& Q) const
{
    BL_PROFILE("EilmerGasModel::prim2cons");

    // get the mass fractions
    Vector<Real> alpha;
    get_alpha_fractions_from_prim(Q, alpha);

    Vector<Real> massf(n_species, 0.0);
    for (int i=0; i<n_species; ++i) {
        const int alpha_idx = species_info[i].alpha_idx;
        if (alpha_idx > -1) {
            massf[i] = alpha[alpha_idx];
        }
    }

    gas_state_set_array_field(gas_state_id, "massf", massf.data(), n_species);

    Real rho = Q[+HydroDef::PrimIdx::Density];
    Real p = Q[+HydroDef::PrimIdx::Prs];
    Real T = Q[+HydroDef::PrimIdx::Temp];

    if ((rho > 0.0) && (p > 0.0)) {

        gas_state_set_scalar_field(gas_model_id, "rho", rho*MFP::rho_ref);
        gas_state_set_scalar_field(gas_model_id, "p", p*MFP::prs_ref);


        // update model and solve for update
        gas_model_gas_state_update_thermo_from_rhop(gas_model_id, gas_state_id);

        gas_state_get_scalar_field(gas_state_id, "T", &T);

        Q[+HydroDef::PrimIdx::Temp] = T/MFP::T_ref;

    } else if ((p > 0.0) && (T > 0.0)) {
        gas_state_set_scalar_field(gas_model_id, "T", T*MFP::T_ref);
        gas_state_set_scalar_field(gas_model_id, "p", p*MFP::prs_ref);


        // update model and solve for update
        gas_model_gas_state_update_thermo_from_pT(gas_model_id, gas_state_id);

        gas_state_get_scalar_field(gas_state_id, "rho", &rho);

        Q[+HydroDef::PrimIdx::Density] = rho/MFP::rho_ref;
    } else if ((rho > 0.0) && (T > 0.0)) {
        gas_state_set_scalar_field(gas_model_id, "rho", rho*MFP::rho_ref);
        gas_state_set_scalar_field(gas_model_id, "T", T*MFP::T_ref);


        // update model and solve for update
        gas_model_gas_state_update_thermo_from_rhoT(gas_model_id, gas_state_id);

        gas_state_get_scalar_field(gas_state_id, "p", &p);

        Q[+HydroDef::PrimIdx::Prs] = p/MFP::prs_ref;
    }

}

Real EilmerGasModel::get_temperature_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("EilmerGasModel::get_temperature_from_cons");

    set_eilmer_gas_state_from_cons(U);

    Real T;
    gas_state_get_scalar_field(gas_state_id, "T", &T);
    return T/MFP::T_ref;

}

RealArray EilmerGasModel::get_speed_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("EilmerGasModel::get_speed_from_cons");
    Real rho = U[+HydroDef::ConsIdx::Density];
    Real mx = U[+HydroDef::ConsIdx::Xmom];
    Real my = U[+HydroDef::ConsIdx::Ymom];
    Real mz = U[+HydroDef::ConsIdx::Zmom];

    Real rhoinv = 1/rho;
    Real u = mx*rhoinv;
    Real v = my*rhoinv;
    Real w = mz*rhoinv;

    set_eilmer_gas_state_from_cons(U);
    gas_model_gas_state_update_sound_speed(gas_model_id, gas_state_id);

    Real a;
    gas_state_get_scalar_field(gas_state_id, "a", &a);

    a /= MFP::u_ref;

    RealArray s = {AMREX_D_DECL(a + std::abs(u), a + std::abs(v), a + std::abs(w))};

    return s;

}

RealArray EilmerGasModel::get_speed_from_prim(const Vector<Real> &Q) const
{
    BL_PROFILE("EilmerGasModel::get_speed_from_prim");

    set_eilmer_gas_state_from_prim(Q);
    gas_model_gas_state_update_sound_speed(gas_model_id, gas_state_id);

    Real a;
    gas_state_get_scalar_field(gas_state_id, "a", &a);

    a /= MFP::u_ref;

    RealArray s = {AMREX_D_DECL(a + std::abs(Q[+HydroDef::PrimIdx::Xvel]),
                   a + std::abs(Q[+HydroDef::PrimIdx::Yvel]),
                   a + std::abs(Q[+HydroDef::PrimIdx::Zvel]))};


    return s;

}

void EilmerGasModel::write_info(nlohmann::json &js) const
{

    HydroGas::write_info(js);

    js["type"] = tag;

}
#endif
