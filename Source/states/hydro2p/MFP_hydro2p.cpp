#include "MFP_hydro2p.H"

#include "MFP_utility.H"
#include "MFP_global.H"
#include "MFP_lua.H"
#include "MFP_Riemann_solvers.H"
#include "MFP_sources.H"

using GD = GlobalData;

Vector<std::string> Hydro2PState::cons_names = {
    "rho",
    "x_mom",
    "y_mom",
    "z_mom",
    "nrg",
    "p_prs",
    "tracer"
};

Vector<std::string> Hydro2PState::prim_names = {
    "rho",
    "x_vel",
    "y_vel",
    "z_vel",
    "p",
    "pp",
    "alpha"
};

Vector<set_bc> Hydro2PState::bc_set = {
    &set_scalar_bc,
    &set_x_vel_bc,
    &set_y_vel_bc,
    &set_z_vel_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc
};

Vector<int> Hydro2PState::flux_vector_idx = {+Hydro2PState::FluxIdx::Xvel, +Hydro2PState::FluxIdx::Bx};
Vector<int> Hydro2PState::cons_vector_idx = {+Hydro2PState::ConsIdx::Xmom};
Vector<int> Hydro2PState::prim_vector_idx = {+Hydro2PState::PrimIdx::Xvel};

std::string Hydro2PState::tag = "hydro_2p";
bool Hydro2PState::registered = GetStateFactory().Register(Hydro2PState::tag, StateBuilder<Hydro2PState>);

Hydro2PState::Hydro2PState() : State(){}
Hydro2PState::Hydro2PState(const sol::table &def)
{
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
    sum_cons.resize(+ConsIdx::NUM);
}
Hydro2PState::~Hydro2PState(){}

void Hydro2PState::init_from_lua()
{
    BL_PROFILE("Hydro2PState::init_from_lua");

    State::init_from_lua();

    sol::state& lua = GD::lua;

    const sol::table state_def = lua["states"][name];

    GD::num_fluid += 1;

    //
    // get mass, charge, and density
    //
    expand(state_def["mass"], mass);
    expand(state_def["charge"], charge);
    expand(state_def["gamma"], gamma);

    mass_const = mass[0] == mass[1];
    charge_const = charge[0] == charge[1];
    gamma_const = gamma[0] == gamma[1];

    //
    // viscous terms coefficients
    //
    set_viscosity();

    //
    // two pressure relaxation
    //
    pressure_relaxation_rate = state_def["pressure_relaxation"];


    //
    // hydro boundary conditions
    //

    const Vector<std::string> dir_name = {"x", "y", "z"};
    const Vector<std::string> side_name = {"lo", "hi"};
    const Vector<std::string>& hydro_var = get_prim_names();
    const int N = hydro_var.size();

    BoundaryState &bs = boundary_conditions;
    bs.phys_fill_bc.resize(+Hydro2PState::PrimIdx::NUM);

    sol::optional<sol::table> bc_def = state_def["bc"];
    if (!bc_def)
        Warning("Warning: No boundary conditions defined for state '"+name+"', using defaults");


    for (int ax = 0; ax < AMREX_SPACEDIM; ++ax) {
        for (int lh=0; lh<2; ++lh) {

            std::string side_bc = state_def["bc"][dir_name[ax]][side_name[lh]]["fill_hydro_bc"].get_or<std::string>("outflow");
            int i_side_bc = bc_names.at(side_bc);

            // get any custom values/functions
            for (int j=0; j<N; ++j) {

                if (lh==0) {
                    bs.phys_fill_bc[j].setLo(ax,i_side_bc);
                } else {
                    bs.phys_fill_bc[j].setHi(ax,i_side_bc);
                }


                const sol::object v = state_def["bc"][dir_name[ax]][side_name[lh]][hydro_var[j]].get_or(sol::object());
                Optional3D1VFunction f = get_udf(v);
                bs.set(ax, hydro_var[j] ,lh,f);

                // special case for inflow condition
                if (i_side_bc == PhysBCType::inflow && !f.is_valid()) {
                    Abort("Setting 'fill_hydro_bc = inflow' requires all primitive variables to be defined, '" + hydro_var[j] + "' is not defined");
                }
            }

#ifdef AMREX_USE_EB
            bool is_symmetry = (i_side_bc == PhysBCType::symmetry) || (i_side_bc == PhysBCType::slipwall) || (i_side_bc == PhysBCType::noslipwall);
            if (lh==0) {
                bs.eb_bc.setLo(ax,is_symmetry ? BCType::reflect_even : BCType::foextrap);
            } else {
                bs.eb_bc.setHi(ax,is_symmetry ? BCType::reflect_even : BCType::foextrap);
            }
#endif
        }
    }

    // check validity of inflow bc
    boundary_conditions.post_init();


    //
    // shock detector threshold
    //
    set_shock_threshold();

    //
    // particles
    //

    particle_init = state_def["particles"].get_or<std::string>("");


    //
    // positivity
    //
    enforce_positivity = state_def["enforce_positivity"].get_or(0);

}

void Hydro2PState::post_init_from_lua()
{

    State::post_init_from_lua();

    if (associated_state.at(AssociatedType::Field).size() == 1) {
        linked_em = associated_state[AssociatedType::Field][0];
    } else {
        Abort("State '"+name+"' has more than one associated field state");
    }
}

AssociatedType Hydro2PState::get_association_type() const
{
    if (charge[0] < 0.0 && charge[1] < 0.0) {
        return AssociatedType::Electron;
    } else if (charge[0] > 0.0 && charge[1] > 0.0) {
        return AssociatedType::Ion;
    } else if (charge[0] == 0.0 && charge[1] == 0.0) {
        return AssociatedType::Neutral;
    } else {
        Abort("Charge of state '"+name+"' is not uniformly positive, negative or neutral");
        return AssociatedType::isNull; //  keep the syntax checker happy by returning something here
    }
}


Real Hydro2PState::init_from_number_density(std::map<std::string, Real> data)
{
    Real nd = functions["nd"](data);
    Real alpha = functions["alpha"](data);
    Real m = get_mass(alpha);

    return nd*m;
}

void Hydro2PState::set_udf()
{
    using namespace std::placeholders;

    sol::state& lua = GD::lua;

    sol::table state_def = lua["states"][name];

    // check if we have 'value' defined
    const sol::table value = state_def["value"].get_or(sol::table());


    if (!value.valid())
        Abort("State "+name+" does not have 'value' defined for initial conditions");

    // get a list of any initialisation functions that need to be called during the run
    sol::table dynamic = state_def["dynamic"].get_or(sol::table());

    const auto ignore = init_with_value();

    bool check, success;

    const Vector<std::string> &prim_names = get_prim_names();
    for (int i = 0; i<prim_names.size(); ++i) {

        std::string comp = prim_names[i];

        // is there a variable with this name?
        success = false;
        check = value[comp].valid();

        // doesn't exist, is there an alternative?
        if (!check) {

            // use number density instead of density
            if (comp.compare("rho") == 0) {

                check = value["nd"].valid();

                if (!check)
                    Abort("State "+name+" does not have 'rho' or 'nd' defined for initial conditions");

                Optional3D1VFunction nd;
                success = get_udf(value["nd"], nd, 0.0);

                functions["nd"] = nd;

                Optional3D1VFunction rho;

                rho.set_func(std::bind(&Hydro2PState::init_from_number_density, this, _1));

                functions[comp] = rho;
            }

        }

        if (!success) {

            Optional3D1VFunction v;

            success = get_udf(value[comp], v, 0.0);

            if (!success) {
                for (const auto &j : ignore) {
                    if (i == j.first) {
                        v.set_value(j.second);
                        success = true;
                        break;
                    }
                }
            }

            functions[comp] = v;
        }

        if (dynamic.valid()) {
            for (const auto &d : dynamic) {
                if (d.second.as<std::string>().compare(comp) == 0) {
                    dynamic_functions[i] = &functions[comp];
                }
            }
        }

    }

    return;
}

const Vector<std::string>& Hydro2PState::get_cons_names() const
{
    return cons_names;
}

const Vector<std::string>& Hydro2PState::get_prim_names() const
{
    return prim_names;
}

const Vector<set_bc>& Hydro2PState::get_bc_set() const
{
    return bc_set;
}


Real Hydro2PState::get_mass(Real alpha) const
{
    BL_PROFILE("Hydro2PState::get_mass");

    if (mass_const) return mass[0];

    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    return (m0*m1)/(m0*alpha + m1*(1-alpha));
}

Real Hydro2PState::get_mass(const Vector<Real> &U) const
{
    BL_PROFILE("Hydro2PState::get_mass");

    if (mass_const) return mass[0];

    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_mass(alpha);
}

Real Hydro2PState::get_charge(Real alpha) const
{
    BL_PROFILE("Hydro2PState::get_charge");

    if (charge_const) return charge[0];

    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    Real q0 = charge[0];
    Real q1 = charge[1];

    return (alpha*m0*q1 + (1-alpha)*m1*q0)/(m0*alpha + m1*(1-alpha));
}

Real Hydro2PState::get_charge(const Vector<Real> &U) const
{
    BL_PROFILE("Hydro2PState::get_charge");

    if (charge_const) return charge[0];

    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_charge(alpha);
}

Real Hydro2PState::get_gamma(Real alpha) const
{
    BL_PROFILE("Hydro2PState::get_gamma");

    if (gamma_const) return gamma[0];

    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    Real g0 = gamma[0];
    Real g1 = gamma[1];

    Real cp0 = g0/(m0*(g0-1));
    Real cp1 = g1/(m1*(g1-1));

    Real cv0 = 1/(m0*(g0-1));
    Real cv1 = 1/(m1*(g1-1));

    return ((1-alpha)*cp0 + alpha*cp1)/((1-alpha)*cv0 + alpha*cv1);
}

Real Hydro2PState::get_gamma(const Vector<Real> &U) const
{
    BL_PROFILE("Hydro2PState::get_gamma");

    if (gamma_const) return gamma[0];

    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_gamma(alpha);
}

Real Hydro2PState::get_cp(Real alpha) const
{
    BL_PROFILE("Hydro2PState::get_cp");

    if (gamma_const && mass_const) return gamma[0]/(mass[0]*(gamma[0]-1));

    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    Real g0 = gamma[0];
    Real g1 = gamma[1];

    Real cp0 = g0/(m0*(g0-1));
    Real cp1 = g1/(m1*(g1-1));

    return (1-alpha)*cp0 + alpha*cp1;
}

Real Hydro2PState::get_cp(const Vector<Real> &U) const
{
    BL_PROFILE("Hydro2PState::get_cp");

    if (gamma_const && mass_const) return gamma[0]/(mass[0]*(gamma[0]-1));

    Real alpha = get_alpha_from_cons(U);
    return get_cp(alpha);
}


// in place conversion from conserved to primitive
bool Hydro2PState::cons2prim(Vector<Real>& U) const
{
    BL_PROFILE("Hydro2PState::cons2prim");

    Real rhoinv = 1/U[+ConsIdx::Density];
    U[+PrimIdx::Xvel] = U[+ConsIdx::Xmom]*rhoinv;
    U[+PrimIdx::Yvel] = U[+ConsIdx::Ymom]*rhoinv;
    U[+PrimIdx::Zvel] = U[+ConsIdx::Zmom]*rhoinv;
    Real kineng = 0.5*U[+ConsIdx::Density]*(
                      U[+PrimIdx::Xvel]*U[+PrimIdx::Xvel]
                  + U[+PrimIdx::Yvel]*U[+PrimIdx::Yvel]
                  + U[+PrimIdx::Zvel]*U[+PrimIdx::Zvel]);

    U[+PrimIdx::Alpha] *= rhoinv;

    Real g = get_gamma(U[+PrimIdx::Alpha]);

    U[+PrimIdx::Prs] = (U[+ConsIdx::Eden] - kineng)*(g - 1);

    return Hydro2PState::prim_valid(U);
}

// in-place conversion from primitive to conserved variables
void Hydro2PState::prim2cons(Vector<Real>& U) const
{
    BL_PROFILE("Hydro2PState::prim2cons");

    Vector<Real> Q = U;

    U[+ConsIdx::Xmom] *= Q[+PrimIdx::Density];
    U[+ConsIdx::Ymom] *= Q[+PrimIdx::Density];
    U[+ConsIdx::Zmom] *= Q[+PrimIdx::Density];

    Real kineng = 0.5*Q[+PrimIdx::Density]*(
                      Q[+PrimIdx::Xvel]*Q[+PrimIdx::Xvel]
                  + Q[+PrimIdx::Yvel]*Q[+PrimIdx::Yvel]
                  + Q[+PrimIdx::Zvel]*Q[+PrimIdx::Zvel]);

    Real g = get_gamma(Q[+PrimIdx::Alpha]);


    U[+ConsIdx::Eden] = Q[+PrimIdx::Prs]/(g-1) + kineng;

    U[+ConsIdx::Tracer] *= U[+PrimIdx::Density];

}


bool Hydro2PState::prim_valid(Vector<Real> &Q) const
{
    BL_PROFILE("Hydro2PState::prim_valid");

    if (   (Q[+PrimIdx::Density] <= 0.0)
           || (Q[+PrimIdx::Prs] <= 0.0)
           || (Q[+PrimIdx::PrsP] <= 0.0)
           ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

bool Hydro2PState::cons_valid(Vector<Real> &U) const
{
    BL_PROFILE("Hydro2PState::cons_valid");

    if ((U[+ConsIdx::Density] <= 0.0) ||  (U[+ConsIdx::Eden] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

Real Hydro2PState::get_energy_from_cons(const Vector<Real> &U) const
{
    return U[+ConsIdx::Eden];
}

Real Hydro2PState::Hydro2PState::get_temperature_from_cons(const Vector<Real> &U) const
{
    Vector<Real> Q = U;

    cons2prim(Q);

    return get_temperature_from_prim(Q);

}

Real Hydro2PState::Hydro2PState::get_temperature_from_prim(const Vector<Real> &Q) const
{
    Real alpha = Q[+PrimIdx::Alpha];
    Real rho = Q[+PrimIdx::Density];
    Real p = Q[+PrimIdx::Prs];

    Real m = get_mass(alpha);

    return p/(rho/m);

}

dual Hydro2PState::get_temperature_from_prim(const Vector<dual> &Q) const
{
    Real alpha = Q[+PrimIdx::Alpha].val;
    dual rho = Q[+PrimIdx::Density];
    dual p = Q[+PrimIdx::Prs];

    Real m = get_mass(alpha);

    return p/(rho/m);

}

RealArray Hydro2PState::Hydro2PState::get_speed_from_cons(const Vector<Real> &U) const
{
    Vector<Real> Q = U;

    Hydro2PState::cons2prim(Q);

    return get_speed_from_prim(Q);

}

RealArray Hydro2PState::get_speed_from_prim(const Vector<Real> &Q) const
{
    BL_PROFILE("Hydro2PState::get_speed_from_prim");
    Real g = get_gamma(Q[+PrimIdx::Alpha]);

    Real a = std::sqrt(g*Q[+PrimIdx::Prs]/Q[+PrimIdx::Density]);

    RealArray s;
    for (int d = 0; d<AMREX_SPACEDIM; ++d) {
        s[d] = a + std::abs(Q[+PrimIdx::Xvel+d]);
    }

    return s;

}

RealArray Hydro2PState::get_current_from_cons(const Vector<Real> &U) const
{
    BL_PROFILE("Hydro2PState::get_current_from_cons");
    Real q = get_charge(U);
    Real m = get_mass(U);
    Real r = q/m;

    RealArray c = {AMREX_D_DECL(
                   r*U[+ConsIdx::Xmom],
                   r*U[+ConsIdx::Ymom],
                   r*U[+ConsIdx::Zmom]
                   )};

    return c;
}

Real Hydro2PState::local_shock_detector(const Vector<Real> &L,
                                        const Vector<Real> &R) const
{
    BL_PROFILE("Hydro2PState::local_shock_detector");

    Real pL = L[+PrimIdx::Prs];
    Real pR = R[+PrimIdx::Prs];
    Real varphi = std::abs(pR - pL)/(pR + pL);

    Real shk = 0.5 + 0.5*std::tanh(5*(varphi - 0.75*shock_threshold)/shock_threshold);

    return shk;
}

void Hydro2PState::get_state_values(const Box& box,
                                    const FArrayBox& src,
                                    std::map<std::string,FArrayBox>& out,
                                    Vector<std::string>& updated
                                    EB_OPTIONAL(,const FArrayBox& vfrac)
                                    ) const
{
    BL_PROFILE("Hydro2PState::get_state_values");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
    Array4<const Real> const& vf4 = vfrac.array();
#endif

    updated.resize(0);

    // check conserved variables
    std::map<std::string,int> cons_tags;
    for (int i=0; i<n_cons(); ++i) {
        const std::string s = cons_names[i];
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        cons_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // check primitive variables
    std::map<std::string,int> prim_tags;
    for (int i=0; i<n_prim(); ++i) {
        const std::string s = prim_names[i];
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        prim_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // additional variables

    Vector<std::string> other;

    const std::string charge_name = "charge-"+name;
    bool load_charge = out.find(charge_name) != out.end();
    if (load_charge) other.push_back(charge_name);

    const std::string mass_name = "mass-"+name;
    bool load_mass = out.find(mass_name) != out.end();
    if (load_mass) other.push_back(mass_name);

    const std::string gamma_name = "gamma-"+name;
    bool load_gamma = out.find(gamma_name) != out.end();
    if (load_gamma) other.push_back(gamma_name);

#ifdef AMREX_USE_EB
    const std::string vfrac_name = "vfrac-"+name;
    bool load_vfrac = out.find(vfrac_name) != out.end();
    if (load_vfrac) other.push_back(vfrac_name);
#endif

    updated.insert(updated.end(), other.begin(), other.end());

    std::map<std::string,Array4<Real>> out4;
    for (const std::string& s : updated) {
        out[s].resize(box, 1);
        out4[s] = out[s].array();
    }

    // temporary storage for retrieving the state data
    Vector<Real> S;

    Array4<const Real> const& src4 = src.array();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (vf4(i,j,k) == 0.0) {
                    for (const std::string& s : updated) {
                        out4[s](i,j,k) = 0.0;
                    }
                    continue;
                }
#endif

                S = get_state_vector(src, i, j, k);

                if (load_charge) out4[charge_name](i,j,k) = get_charge(S);
                if (load_mass)   out4[mass_name](i,j,k)   = get_mass(S);
                if (load_gamma)  out4[gamma_name](i,j,k)  = get_gamma(S);
#ifdef AMREX_USE_EB
                if (load_vfrac)  out4[vfrac_name](i,j,k)  = vf4(i,j,k);
#endif

                if (!cons_tags.empty()) {
                    for (const auto& var : cons_tags) {
                        out4[var.first](i,j,k) = S[var.second];
                    }
                }

                if (!prim_tags.empty()) {
                    cons2prim(S);

                    for (const auto& var : prim_tags) {
                        out4[var.first](i,j,k) = S[var.second];
                    }
                }
            }
        }
    }

    return;
}

void Hydro2PState::calc_velocity(const Box& box,
                                 FArrayBox& src,
                                 FArrayBox& vel
                                 EB_OPTIONAL(,const EBCellFlagFab& flag)
                                 ) const
{
    BL_PROFILE("Hydro2PState::calc_velocity");

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<Real> const& src4 = src.array();
    Array4<Real> const& vel4 = vel.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    int nc = vel.nComp();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    for (int n=0; n<nc; ++n) {
                        vel4(i,j,k,n) = 0.0;
                    }
                }
#endif

                Real invrho = 1.0/src4(i,j,k,+ConsIdx::Density);

                for (int n=0; n<nc; ++n) {
                    vel4(i,j,k,n) = src4(i,j,k,+ConsIdx::Xmom+n)*invrho;
                }
            }
        }
    }

    return;
}


// given all of the available face values load the ones expected by the flux calc into a vector
Vector<Real> Hydro2PState::load_state_for_flux(Vector<Array4<const Real>> &face,
                                               int i, int j, int k) const
{
    BL_PROFILE("Hydro2PState::load_state_for_flux");

    const int nf = +FluxIdx::NUM;
    Vector<Real> S(nf);

    // first get the primitives of this state
    Array4<const Real> const &f4 = face[global_idx];

    S[+FluxIdx::Density] = f4(i,j,k,+PrimIdx::Density);
    S[+FluxIdx::Xvel] = f4(i,j,k,+PrimIdx::Xvel);
    S[+FluxIdx::Yvel] = f4(i,j,k,+PrimIdx::Yvel);
    S[+FluxIdx::Zvel] = f4(i,j,k,+PrimIdx::Zvel);
    S[+FluxIdx::Prs] = f4(i,j,k,+PrimIdx::Prs);
    S[+FluxIdx::PrsP] = f4(i,j,k,+PrimIdx::PrsP);
    S[+FluxIdx::Alpha] = f4(i,j,k,+PrimIdx::Alpha);
    S[+FluxIdx::Gamma] = get_gamma(S[+FluxIdx::Alpha]);

    // then get the magnetic field from the attached field state
    S[+FluxIdx::Bx] = face[linked_em](i,j,k,+FieldState::PrimIdx::Bx);
    S[+FluxIdx::By] = face[linked_em](i,j,k,+FieldState::PrimIdx::By);
    S[+FluxIdx::Bz] = face[linked_em](i,j,k,+FieldState::PrimIdx::Bz);


    return S;
}

void Hydro2PState::write_info(nlohmann::json& js) const
{

    State::write_info(js);

    js["mass"] = mass;
    js["charge"] = charge;
    js["gamma"] = gamma;

    js["pressure_relaxation_rate"] = pressure_relaxation_rate;

}
