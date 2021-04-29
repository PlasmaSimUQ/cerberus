#include "MFP_mhd.H"

#include "MFP_utility.H"
#include "MFP_global.H"
#include "MFP_lua.H"
#include "MFP_Riemann_solvers.H"

using GD = GlobalData;

Vector<std::string> MhdState::cons_names = {
    "rho",
    "x_mom",
    "y_mom",
    "z_mom",
    "nrg",
    "tracer",
    "x_B",
    "y_B",
    "z_B",
    "psi"
};

Vector<std::string> MhdState::prim_names = {
    "rho",
    "x_vel",
    "y_vel",
    "z_vel",
    "p",
    "alpha",
    "x_B",
    "y_B",
    "z_B",
    "psi"
};

Vector<set_bc> MhdState::bc_set = {
    &set_scalar_bc,
    &set_x_vel_bc,
    &set_y_vel_bc,
    &set_z_vel_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_x_B_bc,
    &set_y_B_bc,
    &set_z_B_bc,
    &set_scalar_bc
};

Vector<int> MhdState::flux_vector_idx = {+MhdState::FluxIdx::Xvel, +MhdState::FluxIdx::Bx};
Vector<int> MhdState::cons_vector_idx = {+MhdState::ConsIdx::Xmom, +MhdState::ConsIdx::Bx};
Vector<int> MhdState::prim_vector_idx = {+MhdState::PrimIdx::Xvel, +MhdState::PrimIdx::Bx};

std::string MhdState::tag = "mhd";
bool MhdState::registered = GetStateFactory().Register(MhdState::tag, StateBuilder<MhdState>);

MhdState::MhdState() : State(){}
MhdState::MhdState(const sol::table &def)
{
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
    sum_cons.resize(+ConsIdx::NUM);
}
MhdState::~MhdState(){}



void MhdState::init_from_lua()
{
    BL_PROFILE("MhdState::init_from_lua");
    State::init_from_lua();

    sol::state& lua = GD::lua;
    const sol::table state_def = lua["states"][name];

    GD::num_fluid += 1;

    //
    // get mass, and gamma
    //
    expand(state_def["mass"], mass);
    expand(state_def["gamma"], gamma);

    mass_const = mass[0] == mass[1];
    gamma_const = gamma[0] == gamma[1];

    //
    // boundary conditions
    //

    const Vector<std::string> dir_name = {"x", "y", "z"};
    const Vector<std::string> side_name = {"lo", "hi"};

    // mapping between name and index for various groupings
    std::map<std::string,std::map<std::string, int>> bc2idx;
    for (int i=+MhdState::PrimIdx::Density; i<=+MhdState::PrimIdx::Alpha; ++i) {
        bc2idx["fill_hydro_bc"][prim_names[i]] = i;
    }
    for (int i=+MhdState::PrimIdx::Bx; i<=+MhdState::PrimIdx::Bz; ++i) {
        bc2idx["fill_B_bc"][prim_names[i]] = i;
    }
    bc2idx["fill_psi_bc"] = {{"psi",+MhdState::PrimIdx::psi}};

    BoundaryState &bs = boundary_conditions;
    bs.phys_fill_bc.resize(+MhdState::PrimIdx::NUM);

    for (int ax = 0; ax < AMREX_SPACEDIM; ++ax) {
        for (int lh=0; lh<2; ++lh) {
            for (const auto &bc : bc2idx) {

                // get the base boundary condition for cell centered values
                std::string side_bc = state_def["bc"][dir_name[ax]][side_name[lh]][bc.first].get_or<std::string>("outflow");
                int i_side_bc = bc_names.at(side_bc);


                // fill in the bc list for AMReX as well as gather any custom values/functions
                for (const auto &var : bc.second) {

                    if (lh==0) {
                        bs.phys_fill_bc[var.second].setLo(ax,i_side_bc);
                    } else {
                        bs.phys_fill_bc[var.second].setHi(ax,i_side_bc);
                    }

                    const sol::object bcv = state_def["bc"][dir_name[ax]][side_name[lh]][var.first].get_or(sol::object());

                    Optional3D1VFunction v;

                    // special case for phi and psi (set to zero in boundary)
                    if (var.second == +PrimIdx::psi) {
                        get_udf(bcv,v,0.0);
                    } else {
                        v = get_udf(bcv);
                    }

                    bs.set(ax, prim_names[var.second],lh,v);

                    // special case for inflow condition
                    if (i_side_bc == PhysBCType::inflow && !v.is_valid()) {
                        Abort("Setting '"+bc.first+" = inflow' requires all primitive variables to be defined, '" + var.first + "' is not defined");
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
    }

    // check validity of inflow bc
    boundary_conditions.post_init();

    //
    // divergence handling
    //

    relative_div_speed = state_def["div_transport"].get_or(0.0);

    // or we can use the projectuon method for divergence error control
    project_divergence = state_def["project_divergence"].get_or(0);



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

const Vector<std::string>& MhdState::get_cons_names() const
{
    return cons_names;
}

const Vector<std::string>& MhdState::get_prim_names() const
{
    return prim_names;
}

const Vector<set_bc>& MhdState::get_bc_set() const
{
    return bc_set;
}

Real MhdState::get_mass(Real alpha) const
{
    BL_PROFILE("MhdState::get_mass");

    if (mass_const) return mass[0];

    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    return (m0*m1)/(m0*alpha + m1*(1-alpha));
}

Real MhdState::get_mass(const Vector<Real> &U) const
{
    BL_PROFILE("MhdState::get_mass");

    if (mass_const) return mass[0];

    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_mass(alpha);
}

Real MhdState::get_gamma(Real alpha) const
{
    BL_PROFILE("MhdState::get_gamma");

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

Real MhdState::get_gamma(const Vector<Real> &U) const
{
    BL_PROFILE("MhdState::get_gamma");

    if (gamma_const) return gamma[0];

    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_gamma(alpha);
}

Real MhdState::get_cp(Real alpha) const
{
    BL_PROFILE("MhdState::get_cp");

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

Real MhdState::get_cp(const Vector<Real> &U) const
{
    BL_PROFILE("MhdState::get_cp");

    if (gamma_const && mass_const) return gamma[0]/(mass[0]*(gamma[0]-1));

    Real alpha = get_alpha_from_cons(U);
    return get_cp(alpha);
}



// in place conversion from conserved to primitive
bool MhdState::cons2prim(Vector<Real>& U) const
{
    BL_PROFILE("MhdState::cons2prim");

    Real rho = U[+ConsIdx::Density];
    Real rhoinv = 1/rho;
    Real u = U[+ConsIdx::Xmom]*rhoinv;
    Real v = U[+ConsIdx::Ymom]*rhoinv;
    Real w = U[+ConsIdx::Zmom]*rhoinv;

    Real alpha = U[+ConsIdx::Tracer]*rhoinv;

    Real g = get_gamma(alpha);

    Real Bx = U[+ConsIdx::Bx];
    Real By = U[+ConsIdx::By];
    Real Bz = U[+ConsIdx::Bz];

    Real nrg = U[+ConsIdx::Eden] - 0.5*(rho*(u*u + v*v + w*w) + Bx*Bx + By*By + Bz*Bz);

    U[+PrimIdx::Xvel] = u;
    U[+PrimIdx::Yvel] = v;
    U[+PrimIdx::Zvel] = w;
    U[+PrimIdx::Prs] = nrg*(g-1);
    U[+PrimIdx::Alpha] = alpha;

    return prim_valid(U);
}

// in-place conversion from primitive to conserved variables
void MhdState::prim2cons(Vector<Real>& U) const
{
    BL_PROFILE("MhdState::prim2cons");

    Real rho = U[+PrimIdx::Density];
    Real mx = rho*U[+PrimIdx::Xvel];
    Real my = rho*U[+PrimIdx::Yvel];
    Real mz = rho*U[+PrimIdx::Zvel];
    Real p = U[+PrimIdx::Prs];
    Real alpha = U[+PrimIdx::Alpha];

    Real Bx = U[+PrimIdx::Bx];
    Real By = U[+PrimIdx::By];
    Real Bz = U[+PrimIdx::Bz];

    Real g = get_gamma(alpha);

    U[+ConsIdx::Xmom] = mx;
    U[+ConsIdx::Ymom] = my;
    U[+ConsIdx::Zmom] = mz;
    U[+ConsIdx::Eden] = p/(g - 1) + 0.5*(mx*mx + my*my + mz*mz)/rho + 0.5*(Bx*Bx + By*By + Bz*Bz);
    U[+ConsIdx::Tracer] *= rho;
}


bool MhdState::prim_valid(Vector<Real> &Q) const
{
    BL_PROFILE("MhdState::prim_valid");
    if ((Q[+PrimIdx::Density] <= 0.0) ||  (Q[+PrimIdx::Prs] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

bool MhdState::cons_valid(Vector<Real> &U) const
{
    BL_PROFILE("MhdState::cons_valid");
    if ((U[+ConsIdx::Density] <= 0.0) ||  (U[+ConsIdx::Eden] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

Real MhdState::get_energy_from_cons(const Vector<Real> &U) const
{
    return U[+ConsIdx::Eden];
}

Real MhdState::get_temperature_from_cons(const Vector<Real> &U) const
{
    Vector<Real> Q = U;

    cons2prim(Q);

    return get_temperature_from_prim(Q);
}

Real MhdState::get_temperature_from_prim(const Vector<Real> &Q) const
{
    Real alpha = Q[+PrimIdx::Alpha];
    Real rho = Q[+PrimIdx::Density];
    Real prs = Q[+PrimIdx::Prs];

    Real m = get_mass(alpha);

    return prs/(rho/m);

}

RealArray MhdState::get_speed_from_cons(const Vector<Real>& U) const
{
    Vector<Real> Q = U;

    cons2prim(Q);

    return get_speed_from_prim(Q);

}

RealArray MhdState::get_speed_from_prim(const Vector<Real>& Q) const
{
    BL_PROFILE("MhdState::get_speed_from_prim");
    Real g = get_gamma(Q[+PrimIdx::Alpha]);

    Real rho = Q[+PrimIdx::Density];
    Real p = Q[+PrimIdx::Prs];

    Real Bx = Q[+PrimIdx::Bx];
    Real By = Q[+PrimIdx::By];
    Real Bz = Q[+PrimIdx::Bz];

    Real B = std::sqrt(Bx*Bx + By*By + Bz*Bz);

    Real cf = std::sqrt( (g*p + B) / rho);

    RealArray s;
    for (int d = 0; d<AMREX_SPACEDIM; ++d) {
        s[d] = cf + std::abs(Q[+PrimIdx::Xvel+d]);
    }

    return s ;
}

Real MhdState::local_shock_detector(const Vector<Real> &L,
                                    const Vector<Real> &R) const
{
    BL_PROFILE("MhdState::local_shock_detector");
    Real pL = L[+PrimIdx::Prs];
    Real pR = R[+PrimIdx::Prs];
    Real varphi = std::abs(pR - pL)/(pR + pL);

    Real shk = 0.5 + 0.5*std::tanh(5*(varphi - 0.75*shock_threshold)/shock_threshold);

    return shk;
}

void MhdState::get_state_values(const Box& box,
                                const FArrayBox& src,
                                std::map<std::string,FArrayBox>& out,
                                Vector<std::string>& updated
                                EB_OPTIONAL(,const FArrayBox& vfrac)
                                ) const
{
    BL_PROFILE("MhdState::get_state_values");
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

void MhdState::calc_velocity(const Box& box,
                             FArrayBox& src,
                             FArrayBox& vel
                             EB_OPTIONAL(,const EBCellFlagFab& flag)
                             ) const
{
    BL_PROFILE("MhdState::calc_velocity");
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
Vector<Real> MhdState::load_state_for_flux(Vector<Array4<const Real>> &face,
                                               int i, int j, int k) const
{
    BL_PROFILE("MhdState::load_state_for_flux");

    const int nf = +FluxIdx::NUM;
    Vector<Real> S(nf);

    // first get the primitives of this state
    Array4<const Real> const &f4 = face[global_idx];

    S[+FluxIdx::Density] = f4(i,j,k,+PrimIdx::Density);
    S[+FluxIdx::Xvel] = f4(i,j,k,+PrimIdx::Xvel);
    S[+FluxIdx::Yvel] = f4(i,j,k,+PrimIdx::Yvel);
    S[+FluxIdx::Zvel] = f4(i,j,k,+PrimIdx::Zvel);
    S[+FluxIdx::Prs] = f4(i,j,k,+PrimIdx::Prs);
    S[+FluxIdx::Alpha] = f4(i,j,k,+PrimIdx::Alpha);
    S[+FluxIdx::Gamma] = get_gamma(S[+FluxIdx::Alpha]);
    S[+FluxIdx::Mass] = get_mass(S[+FluxIdx::Alpha]);
    S[+FluxIdx::Bx] = f4(i,j,k,+PrimIdx::Bx);
    S[+FluxIdx::By] = f4(i,j,k,+PrimIdx::By);
    S[+FluxIdx::Bz] = f4(i,j,k,+PrimIdx::Bz);
    S[+FluxIdx::psi] = f4(i,j,k,+PrimIdx::psi);

    return S;
}

void MhdState::update_div_clean(const Real* dx)
{

    if (relative_div_speed <= 0) {
        div_speed = -relative_div_speed; // hack for testing purposes
        return;
    }

    // get the maximum speed

    Real ch_B = GD::cfl*dx[0]/dt;
    for (int d=1; d<AMREX_SPACEDIM; ++d) {
        ch_B = std::min(ch_B, GD::cfl*dx[d]/dt);
    }

    div_speed = relative_div_speed*ch_B;

    return;
}

void MhdState::write_info(nlohmann::json& js) const
{

    State::write_info(js);


    js["mass"] = mass;
    js["gamma"] = gamma;

    js["div_transport"] = relative_div_speed;

}
