#include "MFP_viscous.H"

#include <cmath>

#include "MFP_global.H"
#include "MFP_source.H"

using GD = GlobalData;

Viscous::Viscous(){}
Viscous::~Viscous(){}

int Viscous::get_type(){return -1;}
int Viscous::get_num(){return 0;}

void Viscous::get_neutral_coeffs(const Vector<Real> &Q, Real &T, Real &mu, 
                                 Real &kappa){return;}


Real Viscous::get_max_speed(const Vector<Vector<Real>> &U){return 0.0;}

void Viscous::update_linked_states()
{
    State& istate = GD::get_state(idx);
    istate.set_num_grow(2);
}

PhysicsFactory<Viscous>& GetViscousFactory() {
    static PhysicsFactory<Viscous> F;
    return F;
}

std::string Viscous::str() const
{
    std::stringstream msg;

    const auto coeffs = get_refs();

    msg << get_tag() << "(";

    for (const auto& cf : coeffs) {
        msg << cf.first << "=" << cf.second << ", ";
    }

    msg.seekp(-2,msg.cur);
    msg << ")";

    return msg.str();
}

// ====================================================================================

Sutherland::Sutherland(){}
Sutherland::~Sutherland(){}

std::string Sutherland::tag = "Sutherland";
bool Sutherland::registered = GetViscousFactory().Register(Sutherland::tag, ViscousBuilder<Sutherland>);

Sutherland::Sutherland(const int global_idx, const sol::table& def)
{

    // check valid state
    if (!valid_state(global_idx))
        Abort("Sutherland viscosity is not valid for "+GD::get_state(global_idx).name);

    Real mu_ref = def["mu0"];
    Real T_ref = def["T0"];
    Real S_ref = def["S"];
    Prandtl = def["Pr"];
    cfl = def.get_or("cfl",1.0);

    idx = global_idx;

    Real r_r = GD::rho_ref;
    Real u_r = GD::u_ref;
    Real x_r = GD::x_ref;
    Real T_r = GD::T_ref;

    mu_0 = mu_ref/(r_r*x_r*u_r);
    T0 = T_ref/T_r;
    S = S_ref/T_r;

    if ((mu_0 <= 0) || (T0 <= 0) || (S <= 0)) {
        amrex::Abort("Sutherland input coefficients non-physical");
    }
}

int Sutherland::get_type(){return Neutral;}
int Sutherland::get_num(){return NUM_NEUTRAL_DIFF_COEFFS;}

void Sutherland::get_neutral_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa)
{
    BL_PROFILE("Sutherland::get_neutral_coeffs");
    State &istate = GD::get_state(idx);

    T = istate.get_temperature_from_prim(Q);
    Real alpha = istate.get_alpha_from_prim(Q);
    Real cp = istate.get_cp(alpha);

    Real T_ = T/T0;
    mu = mu_0*T_*sqrt(T_)*(T0+S)/(T+S);
    kappa = mu*cp/Prandtl;

    return;
}


Real Sutherland::get_max_speed(const Vector<Vector<amrex::Real> > &U)
{
    BL_PROFILE("Sutherland::get_max_speed");
    State &istate = GD::get_state(idx);

    Real rho = istate.get_density_from_cons(U[0]);
    Real T = istate.get_temperature_from_cons(U[0]);
    Real gamma = istate.get_gamma(U[0]);


    Real T_ = T/T0;
    Real mu = mu_0*T_*sqrt(T_)*(T0+S)/(T+S);

    return (4*mu*gamma/(Prandtl*rho))/cfl;
}

bool Sutherland::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isHydro) {
        return false;
    }
    return true;
}

// ====================================================================================

PowerLaw::PowerLaw(){}
PowerLaw::~PowerLaw(){}

std::string PowerLaw::tag = "PowerLaw";
bool PowerLaw::registered = GetViscousFactory().Register(PowerLaw::tag, ViscousBuilder<PowerLaw>);

int PowerLaw::get_type(){return Neutral;}
int PowerLaw::get_num(){return NUM_NEUTRAL_DIFF_COEFFS;}

PowerLaw::PowerLaw(const int global_idx, const sol::table& def)
{

    // check valid state
    if (!valid_state(global_idx))
        Abort("Power Law viscosity is not valid for "+GD::get_state(global_idx).name);

    Real mu_ref = def["mu0"];
    Real T_ref = def["T0"];
    n = def["n"];
    Prandtl = def["Pr"];
    cfl = def.get_or("cfl",1.0);

    idx = global_idx;

    Real r_r = GD::rho_ref;
    Real u_r = GD::u_ref;
    Real x_r = GD::x_ref;
    Real T_r = GD::T_ref;

    mu_0 = mu_ref/(r_r*x_r*u_r);
    T0 = T_ref/T_r;


    if ((mu_0 <= 0) || (T0 <= 0) || (n <= 0)) {
        amrex::Abort("Power Law input coefficients non-physical");
    }
}

void PowerLaw::get_neutral_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa)
{
    BL_PROFILE("PowerLaw::get_neutral_coeffs");
    State &istate = GD::get_state(idx);

    T = istate.get_temperature_from_prim(Q);
    Real alpha = istate.get_alpha_from_prim(Q);
    Real cp = istate.get_cp(alpha);


    mu = mu_0*pow(T/T0,n);
    kappa = mu*cp/Prandtl;

    return;
}

Real PowerLaw::get_max_speed(const Vector<Vector<amrex::Real> > &U)
{
    BL_PROFILE("PowerLaw::get_max_speed");
    State &istate = GD::get_state(idx);

    Real rho = istate.get_density_from_cons(U[0]);
    Real T = istate.get_temperature_from_cons(U[0]);
    Real mu = mu_0*pow(T/T0,n);
    Real gamma = istate.get_gamma(U[0]);

    return (4*mu*gamma/(Prandtl*rho))/cfl;
}

bool PowerLaw::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isHydro) {
        return false;
    }
    return true;
}

// ====================================================================================

UserDefinedViscosity::UserDefinedViscosity(){}
UserDefinedViscosity::~UserDefinedViscosity(){}

std::string UserDefinedViscosity::tag = "UserDefined";
bool UserDefinedViscosity::registered = GetViscousFactory().Register(UserDefinedViscosity::tag, ViscousBuilder<UserDefinedViscosity>);

int UserDefinedViscosity::get_type(){return Neutral;}
int UserDefinedViscosity::get_num(){return NUM_NEUTRAL_DIFF_COEFFS;}

UserDefinedViscosity::UserDefinedViscosity(const int global_idx, const sol::table& def)
{

    // check valid state
    if (!valid_state(global_idx))
        Abort("Constant viscosity is not valid for "+GD::get_state(global_idx).name);

    mu_0 = def["mu0"];
    Prandtl = def["Pr"];
    cfl = def.get_or("cfl",1.0);

    idx = global_idx;

    if (mu_0 <= 0) {
        amrex::Abort("Constant viscosity input coefficients non-physical");
    }
}

void UserDefinedViscosity::get_neutral_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa)
{
    BL_PROFILE("UserDefinedViscosity::get_neutral_coeffs");
    State &istate = GD::get_state(idx);

    T = istate.get_temperature_from_prim(Q);
    Real alpha = istate.get_alpha_from_prim(Q);
    Real cp = istate.get_cp(alpha);


    mu = mu_0;
    kappa = mu*cp/Prandtl;

    return;
}

Real UserDefinedViscosity::get_max_speed(const Vector<Vector<amrex::Real> >& U)
{
    BL_PROFILE("UserDefinedViscosity::get_max_speed");
    State &istate = GD::get_state(idx);

    Real rho = istate.get_density_from_cons(U[0]);
    Real mu = mu_0;
    Real gamma = istate.get_gamma(U[0]);

    return (4*mu*gamma/(Prandtl*rho))/cfl;
}

bool UserDefinedViscosity::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isHydro) {
        return false;
    }
    return true;
}

