#include "MFP_hydro_viscous.H"

#include "MFP.H"
#include "MFP_state.H"
#include "MFP_hydro.H"

// ====================================================================================

HydroViscous::HydroViscous()
{
    active = false;
}

std::string HydroViscous::str() const
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

Real HydroViscous::get_min_dt(MFP *mfp) const
{
    BL_PROFILE("get_max_speed");
    const HydroState& istate = HydroState::get_state_global(idx);

    MultiFab& data = mfp->get_new_data(istate.data_idx);

    const size_t n_cons = istate.n_cons();
    Vector<Real> U(n_cons);

    Real max_speed = 0;

#ifdef AMREX_USE_EB
    EBData& eb = mfp->get_eb_data(idx);
#endif

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
        // get the EB data required for later calls

        const FArrayBox& vfrac = eb.volfrac[mfi];

        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();
#endif

        Array4<Real const> const dat = data.const_array(mfi);

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) {
                        continue;
                    }
#endif

                    for (int n=0; n<n_cons; ++n) {
                        U[n] = dat(i,j,k,n);
                    }

                    max_speed = amrex::max(max_speed, calc_stability_condition(U));

                }
            }
        }
    }

    Real min_dt = std::numeric_limits<Real>::max();

    const Real* dx = mfp->Geom().CellSize();

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        min_dt = std::min(min_dt, dx[d]*dx[d]/max_speed);
    }

    return min_dt;
}


ClassFactory<HydroViscous>& GetHydroViscousFactory() {
    static ClassFactory<HydroViscous> F;
    return F;
}

// ====================================================================================

Sutherland::Sutherland()
{
    active = true;
}
Sutherland::~Sutherland(){}

std::string Sutherland::tag = "Sutherland";
bool Sutherland::registered = GetHydroViscousFactory().Register(Sutherland::tag, HydroViscousBuilder<Sutherland>);

Sutherland::Sutherland(const int global_idx, const sol::table& def)
{

    Real mu_ref = def["mu0"];
    Real T_ref = def["T0"];
    Real S_ref = def["S"];
    Prandtl = def["Pr"];

    idx = global_idx;

    Real r_r = MFP::rho_ref;
    Real u_r = MFP::u_ref;
    Real x_r = MFP::x_ref;
    Real T_r = MFP::T_ref;

    mu_0 = mu_ref/(r_r*x_r*u_r);
    T0 = T_ref/T_r;
    S = S_ref/T_r;

    if ((mu_0 <= 0) || (T0 <= 0) || (S <= 0)) {
        amrex::Abort("Sutherland input coefficients non-physical");
    }
}

void Sutherland::get_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa) const
{
    BL_PROFILE("Sutherland::get_neutral_coeffs");
    const HydroState& istate = HydroState::get_state_global(idx);

    T = istate.get_temperature_from_prim(Q);
    Real cp = istate.get_cp_from_prim(Q);

    Real T_ = T/T0;
    mu = mu_0*T_*sqrt(T_)*(T0+S)/(T+S);
    kappa = mu*cp/Prandtl;

    return;
}

Real Sutherland::calc_stability_condition(Vector<Real>& U) const
{
    const HydroState& istate = HydroState::get_state_global(idx);

    const Real rho = U[+HydroDef::ConsIdx::Density];
    const Real T = istate.get_temperature_from_cons(U);
    const Real gamma = istate.get_gamma_from_cons(U);


    const Real T_ = T/T0;
    const Real mu = mu_0*T_*sqrt(T_)*(T0+S)/(T+S);

    return (4*mu*gamma/(Prandtl*rho));
}

// ====================================================================================

PowerLaw::PowerLaw()
{
    active = true;
}
PowerLaw::~PowerLaw(){}

std::string PowerLaw::tag = "PowerLaw";
bool PowerLaw::registered = GetHydroViscousFactory().Register(PowerLaw::tag, HydroViscousBuilder<PowerLaw>);

PowerLaw::PowerLaw(const int global_idx, const sol::table& def)
{

    Real mu_ref = def["mu0"];
    Real T_ref = def["T0"];
    n = def["n"];
    Prandtl = def["Pr"];

    idx = global_idx;

    Real r_r = MFP::rho_ref;
    Real u_r = MFP::u_ref;
    Real x_r = MFP::x_ref;
    Real T_r = MFP::T_ref;

    mu_0 = mu_ref/(r_r*x_r*u_r);
    T0 = T_ref/T_r;


    if ((mu_0 <= 0) || (T0 <= 0) || (n <= 0)) {
        amrex::Abort("Power Law input coefficients non-physical");
    }
}

void PowerLaw::get_coeffs(const Vector<Real>& Q, Real &T, Real &mu, Real &kappa) const
{
    BL_PROFILE("PowerLaw::get_neutral_coeffs");
    const HydroState& istate = HydroState::get_state_global(idx);

    T = istate.get_temperature_from_prim(Q);
    Real cp = istate.get_cp_from_prim(Q);


    mu = mu_0*pow(T/T0,n);
    kappa = mu*cp/Prandtl;

    return;
}

Real PowerLaw::calc_stability_condition(Vector<Real>& U) const
{

    const HydroState& istate = HydroState::get_state_global(idx);

    const Real rho = U[+HydroDef::ConsIdx::Density];
    const Real T = istate.get_temperature_from_cons(U);
    const Real gamma = istate.get_gamma_from_cons(U);

    const Real mu = mu_0*pow(T/T0,n);

    return (4*mu*gamma/(Prandtl*rho));
}

// ====================================================================================

UserDefinedViscosity::UserDefinedViscosity()
{
    active = true;
}
UserDefinedViscosity::~UserDefinedViscosity(){}

std::string UserDefinedViscosity::tag = "UserDefined";
bool UserDefinedViscosity::registered = GetHydroViscousFactory().Register(UserDefinedViscosity::tag, HydroViscousBuilder<UserDefinedViscosity>);

UserDefinedViscosity::UserDefinedViscosity(const int global_idx, const sol::table& def)
{
    mu_0 = def["mu0"];
    Prandtl = def["Pr"];

    idx = global_idx;

    if (mu_0 <= 0) {
        amrex::Abort("Constant viscosity input coefficients non-physical");
    }
}

void UserDefinedViscosity::get_coeffs(const Vector<Real>& Q, Real &T, Real &mu, Real &kappa) const
{
    BL_PROFILE("UserDefinedViscosity::get_neutral_coeffs");
    const HydroState& istate = HydroState::get_state_global(idx);

    T = istate.get_temperature_from_prim(Q);
    Real cp = istate.get_cp_from_prim(Q);

    mu = mu_0;
    kappa = mu*cp/Prandtl;

    return;
}

Real UserDefinedViscosity::calc_stability_condition(Vector<Real>& U) const
{

    const HydroState& istate = HydroState::get_state_global(idx);

    const Real rho = U[+HydroDef::ConsIdx::Density];
    const Real gamma = istate.get_gamma_from_cons(U);
    const Real mu = mu_0;

    return (4*mu*gamma/(Prandtl*rho));
}


