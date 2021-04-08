#include "MFP_hydro.H"

#include <functional>

#include "sol.hpp"
#include "MFP_utility.H"
#include "MFP_global.H"
#include "MFP_lua.H"
#include "MFP_Riemann_solvers.H"
#include "MFP_hydro_bc.H"
#include "Eigen"

#ifdef PYTHON
#include "MFP_diagnostics.H"
#endif

using GD = GlobalData;

//------------
// Hydro State

Vector<std::string> HydroState::cons_names = {
    "rho",
    "x_mom",
    "y_mom",
    "z_mom",
    "nrg",
    "tracer"
};

Vector<std::string> HydroState::prim_names = {
    "rho",
    "x_vel",
    "y_vel",
    "z_vel",
    "p",
    "T",
    "alpha"
};

Vector<set_bc> HydroState::bc_set = {
    &set_scalar_bc,
    &set_x_vel_bc,
    &set_y_vel_bc,
    &set_z_vel_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc
};

//-----------------------------------------------------------------------------

std::string HydroState::tag = "hydro";
bool HydroState::registered = GetStateFactory().Register(HydroState::tag, StateBuilder<HydroState>);

HydroState::HydroState() : State(){}
HydroState::HydroState(const sol::table& def)
{
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
    sum_cons.resize(+ConsIdx::NUM);
}
HydroState::~HydroState(){}

void HydroState::init_from_lua()
{

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

    //
    // viscous terms coefficients
    //
    set_viscosity();

    //
    // domain boundary conditions
    //

    const Vector<std::string> dir_name = {"x", "y", "z"};
    const Vector<std::string> side_name = {"lo", "hi"};
    const Vector<std::string>& hydro_var = get_prim_names();
    const int N = hydro_var.size();

    BoundaryState &bs = boundary_conditions;
    bs.phys_fill_bc.resize(+HydroState::PrimIdx::NUM);

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
                bs.set(ax,hydro_var[j],lh,f);

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

    extra_slope_limits = state_def["extra_slope_limits"].get_or(1);

}



#ifdef AMREX_USE_EB
void HydroState::set_eb_bc(const sol::table &bc_def)
{

    std::string bc_type = bc_def.get<std::string>("type");

    if (bc_type == HydroSlipWall::tag) {
        eb_bcs.push_back(std::unique_ptr<BoundaryEB>(new HydroSlipWall(flux_solver.get())));
    } else if (bc_type == HydroNoSlipWall::tag) {
        if (!viscous) {
            Abort("Requested EB bc of type '" + bc_type + "' without defining 'viscosity' for state '" + name + "'");
        }
        eb_bcs.push_back(std::unique_ptr<BoundaryEB>(new HydroNoSlipWall(flux_solver.get(), viscous.get(), bc_def)));
    } else if (bc_type == DirichletWall::tag) {
        eb_bcs.push_back(std::unique_ptr<BoundaryEB>(new DirichletWall(flux_solver.get(), get_prim_names(), get_prim_vector_idx(), bc_def)));
    } else {
        Abort("Requested EB bc of type '" + bc_type + "' which is not compatible with state '" + name + "'");
    }
}
#endif

Real HydroState::init_from_number_density(std::map<std::string, Real> data)
{
    Real nd = functions["nd"](data);
    Real alpha = functions["alpha"](data);
    Real m = get_mass(alpha);

    return nd*m;
}

void HydroState::set_udf()
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

                rho.set_func(std::bind(&HydroState::init_from_number_density, this, _1));

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

AssociatedType HydroState::get_association_type() const
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

const Vector<std::string>& HydroState::get_cons_names() const
{
    return cons_names;
}

const Vector<std::string>& HydroState::get_prim_names() const
{
    return prim_names;
}

const Vector<set_bc>& HydroState::get_bc_set() const
{
    return bc_set;
}


Real HydroState::get_mass(Real alpha) const
{
    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    return (m0*m1)/(m0*alpha + m1*(1-alpha));
}

Real HydroState::get_mass(const Vector<Real> &U) const
{
    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_mass(alpha);
}

Real HydroState::get_charge(Real alpha) const
{
    // clamp alpha
    alpha = clamp(alpha, 0.0, 1.0);

    Real m0 = mass[0];
    Real m1 = mass[1];

    Real q0 = charge[0];
    Real q1 = charge[1];

    return (alpha*m0*q1 + (1-alpha)*m1*q0)/(m0*alpha + m1*(1-alpha));
}

Real HydroState::get_charge(const Vector<Real> &U) const
{
    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_charge(alpha);
}

Real HydroState::get_gamma(Real alpha) const
{
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

Real HydroState::get_gamma(const Vector<Real> &U) const
{
    // clamp alpha
    Real alpha = get_alpha_from_cons(U);
    return get_gamma(alpha);
}

Real HydroState::get_cp(Real alpha) const
{
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

Real HydroState::get_cp(const Vector<Real> &U) const
{
    Real alpha = get_alpha_from_cons(U);
    return get_cp(alpha);
}


// in place conversion from conserved to primitive
bool HydroState::cons2prim(Vector<Real>& U) const
{
    U.resize(+PrimIdx::NUM);

    // move tracer

    Real rhoinv = 1/U[+ConsIdx::Density];
    U[+PrimIdx::Xvel] = U[+ConsIdx::Xmom]*rhoinv;
    U[+PrimIdx::Yvel] = U[+ConsIdx::Ymom]*rhoinv;
    U[+PrimIdx::Zvel] = U[+ConsIdx::Zmom]*rhoinv;
    Real kineng = 0.5*U[+ConsIdx::Density]*(
                      U[+PrimIdx::Xvel]*U[+PrimIdx::Xvel]
                  + U[+PrimIdx::Yvel]*U[+PrimIdx::Yvel]
                  + U[+PrimIdx::Zvel]*U[+PrimIdx::Zvel]);

    U[+PrimIdx::Alpha] = U[+ConsIdx::Tracer]*rhoinv;

    Real g = get_gamma(U[+PrimIdx::Alpha]);

    U[+PrimIdx::Prs] = (U[+ConsIdx::Eden] - kineng)*(g - 1);

    Real m = get_mass(U[+PrimIdx::Alpha]);
    U[+PrimIdx::Temp] = U[+PrimIdx::Prs]*rhoinv*m;

    return prim_valid(U);
}

// in place conversion from conserved to primitive
bool HydroState::cons2prim(Vector<autodiff::dual>& U) const
{
    U.resize(+PrimIdx::NUM); // expand

    autodiff::dual rhoinv = 1/U[+ConsIdx::Density];
    U[+PrimIdx::Xvel] = U[+ConsIdx::Xmom]*rhoinv;
    U[+PrimIdx::Yvel] = U[+ConsIdx::Ymom]*rhoinv;
    U[+PrimIdx::Zvel] = U[+ConsIdx::Zmom]*rhoinv;
    autodiff::dual kineng = 0.5*U[+ConsIdx::Density]*(
                                U[+PrimIdx::Xvel]*U[+PrimIdx::Xvel]
                            + U[+PrimIdx::Yvel]*U[+PrimIdx::Yvel]
                            + U[+PrimIdx::Zvel]*U[+PrimIdx::Zvel]);

    U[+PrimIdx::Alpha] =  U[+ConsIdx::Tracer]*rhoinv;

    Real g = get_gamma(U[+PrimIdx::Alpha].val);

    U[+PrimIdx::Prs] = (U[+ConsIdx::Eden] - kineng)*(g - 1);


    Real m = get_mass(U[+PrimIdx::Alpha].val);
    U[+PrimIdx::Temp] = U[+PrimIdx::Prs]*rhoinv*m;

    if ((U[+PrimIdx::Density].val <= 0.0) ||  (U[+PrimIdx::Prs].val <= 0.0)) {
        return false;
    }
    return true;
}

// in-place conversion from primitive to conserved variables
void HydroState::prim2cons(Vector<Real>& U) const
{

    Real kineng = 0.5*U[+PrimIdx::Density]*(
                      U[+PrimIdx::Xvel]*U[+PrimIdx::Xvel]
                  + U[+PrimIdx::Yvel]*U[+PrimIdx::Yvel]
                  + U[+PrimIdx::Zvel]*U[+PrimIdx::Zvel]);

    U[+ConsIdx::Xmom] = U[+PrimIdx::Xvel]*U[+PrimIdx::Density];
    U[+ConsIdx::Ymom] = U[+PrimIdx::Yvel]*U[+PrimIdx::Density];
    U[+ConsIdx::Zmom] = U[+PrimIdx::Zvel]*U[+PrimIdx::Density];

    Real g = get_gamma(U[+PrimIdx::Alpha]);

    U[+ConsIdx::Eden] = U[+PrimIdx::Prs]/(g-1) + kineng;

    U[+ConsIdx::Tracer] = U[+PrimIdx::Alpha]*U[+PrimIdx::Density];

    U.resize(+ConsIdx::NUM);

}


bool HydroState::prim_valid(Vector<Real> &Q) const
{
    if ((Q[+PrimIdx::Density] <= 0.0) ||  (Q[+PrimIdx::Prs] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

bool HydroState::cons_valid(Vector<Real> &U) const
{
    if ((U[+ConsIdx::Density] <= 0.0) ||  (U[+ConsIdx::Eden] <= 0.0)
            ) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

Real HydroState::get_energy_from_cons(const Vector<Real> &U) const
{
    return U[+ConsIdx::Eden];
}

Real HydroState::get_temperature_from_cons(const Vector<Real> &U) const
{
    Vector<Real> Q = U;

    cons2prim(Q);

    return get_temperature_from_prim(Q);

}

Real HydroState::get_temperature_from_prim(const Vector<Real> &Q) const
{
    return Q[+PrimIdx::Temp];
}

dual HydroState::get_temperature_from_prim(const Vector<dual> &Q) const
{
    dual T = Q[+PrimIdx::Temp];
    return T;
}

RealArray HydroState::get_speed_from_cons(const Vector<Real>& U) const
{
    Vector<Real> Q = U;

    cons2prim(Q);

    return get_speed_from_prim(Q);

}

RealArray HydroState::get_speed_from_prim(const Vector<Real>& Q) const
{

    Real g = get_gamma(Q[+PrimIdx::Alpha]);

    Real a = std::sqrt(g*Q[+PrimIdx::Prs]/Q[+PrimIdx::Density]);

    RealArray s;
    for (int d = 0; d<AMREX_SPACEDIM; ++d) {
        s[d] = a + std::abs(Q[+PrimIdx::Xvel+d]);
    }

    return s;

}

RealArray HydroState::get_current_from_cons(const Vector<Real> &U) const
{
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



void HydroState::calc_reconstruction(const Box& box,
                                     FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                     EB_OPTIONAL(,const EBCellFlagFab &flag)
                                     EB_OPTIONAL(,const FArrayBox &vfrac)
                                     ) const
{

    // if we don't want to apply extra limiting on the slopes (forced to 2nd order)
    // we can use the default reconstruction scheme
    if (!extra_slope_limits) {
        State::calc_reconstruction(box, prim, rlo, rhi EB_OPTIONAL(,flag,vfrac));
    }

    // convert pressure
    const Box &pbox = prim.box();
    const Dim3 p_lo = amrex::lbound(pbox);
    const Dim3 p_hi = amrex::ubound(pbox);

    FArrayBox gamma_minus_one(pbox);
    Array4<Real> const& src4 = prim.array();
    Array4<Real> const& gam4 = gamma_minus_one.array();

#ifdef AMREX_USE_EB
    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);

    Array4<const EBCellFlag> const& f4 = flag.array();
    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    int N = prim.nComp();

    Vector<Real> stencil(reconstruction->stencil_length);
    int offset = reconstruction->stencil_length/2;
    Array<int,3> stencil_index;
    Vector<Real> cell_value(N), cell_slope(N);

    Real rho_lo, rho_hi;
    Real alpha_lo, alpha_hi;
    Real abs_phi, phi_scale, coeff_eps;
    Real gam_lo, gam_hi;

    // make sure our arrays for putting lo and hi reconstructed values into
    // are the corect size
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        rlo[d].resize(box, N);
        rhi[d].resize(box, N);

#ifdef AMREX_USE_EB
        if (check_eb) {
            rlo[d].copy(prim,box);
            rhi[d].copy(prim,box);
        }
#endif
    }

    // change pressure to internal energy
    for     (int k = p_lo.z; k <= p_hi.z; ++k) {
        for   (int j = p_lo.y; j <= p_hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = p_lo.x; i <= p_hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    continue;
                }
#endif

                gam4(i,j,k) = get_gamma(src4(i,j,k,+PrimIdx::Alpha)) - 1.0;

                src4(i,j,k,+PrimIdx::Prs) /= gam4(i,j,k);

            }
        }
    }

    // now do reconstruction

    // cycle over dimensions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        Array4<Real> const& lo4 = rlo[d].array();
        Array4<Real> const& hi4 = rhi[d].array();

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (check_eb) {

                        // covered cell doesn't need calculating
                        if (f4(i,j,k).isCovered()) {
                            continue;
                        }

                        // cell that references a covered cell doesn't need calculating
                        bool skip = false;
                        stencil_index.fill(0);
                        for (int s=0; s<reconstruction->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            // check if any of the stencil values are from a covered cell
                            if (f4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2]).isCovered()) {
                                skip = true;
                                break;
                            }
                        }

                        if (skip) {
                            continue;
                        }

                    }
#endif

                    // cycle over all components
                    for (int n = 0; n<N; ++n) {

                        // fill in the stencil along dimension index
                        stencil_index.fill(0);
                        for (int s=0; s<reconstruction->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            stencil[s] = src4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], n);
                        }

                        // perform reconstruction
                        cell_slope[n] = reconstruction->get_slope(stencil);
                        cell_value[n] = stencil[offset];

                    }


                    // apply corrections to slopes
                    // J. Sci. Comput. (2014) 60:584-611
                    // Robust Finite Volume Schemes for Two-Fluid Plasma Equations

                    Real &rho     = cell_value[+PrimIdx::Density];
                    Real &phi_rho = cell_slope[+PrimIdx::Density];

                    Real &u     = cell_value[+PrimIdx::Xvel];
                    Real &phi_u = cell_slope[+PrimIdx::Xvel];

                    Real &v     = cell_value[+PrimIdx::Yvel];
                    Real &phi_v = cell_slope[+PrimIdx::Yvel];

                    Real &w     = cell_value[+PrimIdx::Zvel];
                    Real &phi_w = cell_slope[+PrimIdx::Zvel];

                    Real &eps     = cell_value[+PrimIdx::Prs];
                    Real &phi_eps = cell_slope[+PrimIdx::Prs];

                    Real &alpha     = cell_value[+PrimIdx::Alpha];
                    Real &phi_alpha = cell_slope[+PrimIdx::Alpha];



                    // correct density slope
                    if (std::abs(phi_rho) > 2*rho) {
                        phi_rho = 2*sign(phi_rho, 0.0)*rho;
                    }

                    // get some face values
                    rho_lo = rho - 0.5*phi_rho;
                    rho_hi = rho + 0.5*phi_rho;

                    alpha_lo = alpha - 0.5*phi_alpha;
                    alpha_hi = alpha + 0.5*phi_alpha;

                    gam_lo = get_gamma(alpha_lo);
                    gam_hi = get_gamma(alpha_hi);

                    abs_phi = phi_u*phi_u + phi_v*phi_v + phi_w*phi_w;

                    // correct velocity slope
                    Real eps_face = eps - 0.5*std::abs(phi_eps);

                    if (eps_face <= 0.0) {
                        // if the reconstructed face value goes non-physical
                        // just set back to first order with zero slope
                        phi_u = 0.0;
                        phi_v = 0.0;
                        phi_w = 0.0;
                        phi_eps = 0.0;
                    } else {
                        coeff_eps = (rho/(rho_lo*rho_hi))*eps_face;
                        if ((0.125*abs_phi) > coeff_eps) {
                            phi_scale = sqrt(abs_phi);
                            coeff_eps = sqrt(8*coeff_eps);
                            phi_u = (phi_u/phi_scale)*coeff_eps;
                            phi_v = (phi_v/phi_scale)*coeff_eps;
                            phi_w = (phi_w/phi_scale)*coeff_eps;
                        }
                        // update eps
                        abs_phi = phi_u*phi_u + phi_v*phi_v + phi_w*phi_w;
                        eps -= (rho_lo*rho_hi/rho)*0.125*abs_phi;
                    }



                    // density
                    lo4(i,j,k,+PrimIdx::Density) = rho_lo;
                    hi4(i,j,k,+PrimIdx::Density) = rho_hi;

                    // x - velocity
                    lo4(i,j,k,+PrimIdx::Xvel) = u - 0.5*(rho_hi/rho)*phi_u;
                    hi4(i,j,k,+PrimIdx::Xvel) = u + 0.5*(rho_lo/rho)*phi_u;

                    // y - velocity
                    lo4(i,j,k,+PrimIdx::Yvel) = v - 0.5*(rho_hi/rho)*phi_v;
                    hi4(i,j,k,+PrimIdx::Yvel) = v + 0.5*(rho_lo/rho)*phi_v;

                    // z - velocity
                    lo4(i,j,k,+PrimIdx::Zvel) = w - 0.5*(rho_hi/rho)*phi_w;
                    hi4(i,j,k,+PrimIdx::Zvel) = w + 0.5*(rho_lo/rho)*phi_w;

                    // epsilon -> pressure
                    lo4(i,j,k,+PrimIdx::Prs) = (eps - 0.5*phi_eps)*(gam_lo - 1.0);
                    hi4(i,j,k,+PrimIdx::Prs) = (eps + 0.5*phi_eps)*(gam_hi - 1.0);

                    Real prs = lo4(i,j,k,+PrimIdx::Prs);

                    // tracer
                    lo4(i,j,k,+PrimIdx::Alpha) = alpha_lo;
                    hi4(i,j,k,+PrimIdx::Alpha) = alpha_hi;

                    // Temperature (calculate from pressure and density)
                    lo4(i,j,k,+PrimIdx::Temp) = lo4(i,j,k,+PrimIdx::Prs)/(rho_lo/get_mass(alpha_lo));
                    hi4(i,j,k,+PrimIdx::Temp) = hi4(i,j,k,+PrimIdx::Prs)/(rho_hi/get_mass(alpha_hi));

                }
            }
        }
    }


    // convert back to pressure
    for     (int k = p_lo.z; k <= p_hi.z; ++k) {
        for   (int j = p_lo.y; j <= p_hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = p_lo.x; i <= p_hi.x; ++i) {
#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif
                src4(i,j,k,+PrimIdx::Prs) *= gam4(i,j,k);



            }
        }
    }

    return;
}


Real HydroState::local_shock_detector(const Vector<Real> &L,
                                      const Vector<Real> &R) const
{
    Real pL = L[+PrimIdx::Prs];
    Real pR = R[+PrimIdx::Prs];
    Real varphi = std::abs(pR - pL)/(pR + pL);

    Real shk = 0.5 + 0.5*std::tanh(5*(varphi - 0.75*shock_threshold)/shock_threshold);

    return shk;
}

void HydroState::get_state_values(const Box& box,
                                  const FArrayBox& src,
                                  std::map<std::string,FArrayBox>& out,
                                  Vector<std::string>& updated
                                  EB_OPTIONAL(,const FArrayBox& vfrac)
                                  ) const
{
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
        if (s == cons_names[0]) continue;
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
        out[s].setVal(0.0);
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

void HydroState::calc_velocity(const Box& box,
                               FArrayBox& src,
                               FArrayBox& vel
                               EB_OPTIONAL(,const EBCellFlagFab& flag)
                               ) const
{
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
                    continue;
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

void HydroState::calc_viscous_fluxes(const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const Box& pbox,
                                     const Vector<FArrayBox> &prim,
                                     #ifdef AMREX_USE_EB
                                     const EBCellFlagFab& flag,
                                     #endif
                                     const Real* dx) const
{

    // now calculate viscous fluxes and load them into the flux arrays
    if (viscous) {

        Viscous &V = *viscous;
        switch (V.get_type()) {
            case Viscous::Neutral :

//                plot_FAB_2d(flag, "flag", false);
//                plot_FAB_2d(prim[global_idx], 0, "prim density", false, true);

                calc_neutral_viscous_fluxes(box,
                                            fluxes,
                                            pbox,
                                            prim[global_idx],
                                            EB_OPTIONAL(flag,)
                                            dx);
                break;
            case Viscous::Ion :
                calc_ion_viscous_fluxes(box,
                                        fluxes,
                                        pbox,
                                        prim,
                                        EB_OPTIONAL(flag,)
                                        dx);
                break;
            case Viscous::Electron :
                calc_electron_viscous_fluxes(box,
                                             fluxes,
                                             pbox,
                                             prim,
                                             EB_OPTIONAL(flag,)
                                             dx);
                break;
            default :
                break;
        }
    }

}

void HydroState::calc_neutral_diffusion_terms(const Box& box,
                                              const FArrayBox& prim,
                                              FArrayBox& diff
                                              EB_OPTIONAL(,const EBCellFlagFab& flag)
                                              ) const
{

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& prim4 = prim.array();
    Array4<Real> const& d4 = diff.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Viscous &V = *viscous;

    int np = n_prim();

    Vector<Real> Q(np);
    Real T, mu, kappa;

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                for (int n=0; n<np; ++n) {
                    Q[n] = prim4(i,j,k,n);
                }

                V.get_neutral_coeffs(Q, T, mu, kappa);


                d4(i,j,k,Viscous::NeutralTemp) = T;
                d4(i,j,k,Viscous::NeutralKappa) = kappa;
                d4(i,j,k,Viscous::NeutralMu) = mu;
            }
        }
    }

    return;
}

#ifdef AMREX_USE_EB
void HydroState::calc_neutral_viscous_fluxes_eb(const Box& box, Array<FArrayBox,
                                                AMREX_SPACEDIM> &fluxes,
                                                const Box& pbox,
                                                const FArrayBox &prim,
                                                EB_OPTIONAL(const EBCellFlagFab& flag,)
                                                const Real* dx) const
{
    FArrayBox diff(pbox, Viscous::NUM_NEUTRAL_DIFF_COEFFS);
    calc_neutral_diffusion_terms(pbox,
                                 prim,
                                 diff
                                 EB_OPTIONAL(,flag)
                                 );

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<Real, AMREX_SPACEDIM> dxinv;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        dxinv[d] = 1/dx[d];
    }

    Array4<const Real> const& p4 = prim.array();
    Array4<const Real> const& d4 = diff.array();

    Array4<const EBCellFlag> const& f4 = flag.array();

    Real dudx=0, dudy=0, dudz=0, dvdx=0, dvdy=0, dwdx=0, dwdz=0, divu=0;
    const Real two_thirds = 2/3;

    const Array<Real,3>  weights = {0.0, 1.0, 0.5};
    Real whi, wlo;

    Real tauxx, tauxy, tauxz, dTdx, muf;

    int iTemp = Viscous::NeutralTemp;
    int iMu = Viscous::NeutralMu;
    int iKappa = Viscous::NeutralKappa;

    Vector<int> prim_vel_id = get_prim_vector_idx();

    int Xvel = prim_vel_id[0] + 0;
    int Yvel = prim_vel_id[0] + 1;
    int Zvel = prim_vel_id[0] + 2;

    Vector<int> cons_vel_id = get_cons_vector_idx();

    int Xmom = cons_vel_id[0] + 0;
    int Ymom = cons_vel_id[0] + 1;
    int Zmom = cons_vel_id[0] + 2;

    Vector<int> nrg_id = get_nrg_idx();
    int Eden = nrg_id[0];


    // X - direction
    Array4<Real> const& fluxX = fluxes[0].array();
    for     (int k = lo.z-AMREX_D_PICK(0,0,1); k <= hi.z+AMREX_D_PICK(0,0,1); ++k) {
        for   (int j = lo.y-AMREX_D_PICK(0,1,1); j <= hi.y+AMREX_D_PICK(0,1,1); ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x + 1; ++i) {

                bool covered = f4(i,j,k).isCovered();
                bool connected = f4(i,j,k).isConnected(-1,0,0);
                bool other_covered = f4(i-1,j,k).isCovered();

                // only calculate fluxes for fluid cells and between cells that are connected
                if (covered || other_covered || !connected)
                    continue;

                dTdx = (d4(i,j,k,iTemp) - d4(i-1,j,k,iTemp))*dxinv[0];

                dudx = (p4(i,j,k,Xvel) - p4(i-1,j,k,Xvel))*dxinv[0];
                dvdx = (p4(i,j,k,Yvel) - p4(i-1,j,k,Yvel))*dxinv[0];
                dwdx = (p4(i,j,k,Zvel) - p4(i-1,j,k,Zvel))*dxinv[0];

#if AMREX_SPACEDIM >= 2

                const int jhip = j + (int)f4(i,j,k).isConnected(0, 1,0);
                const int jhim = j - (int)f4(i,j,k).isConnected(0,-1,0);
                const int jlop = j + (int)f4(i-1,j,k).isConnected(0, 1,0);
                const int jlom = j - (int)f4(i-1,j,k).isConnected(0,-1,0);
                whi = weights[jhip-jhim];
                wlo = weights[jlop-jlom];
                dudy = (0.5*dxinv[1]) * ((p4(i  ,jhip,k,Xvel)-p4(i  ,jhim,k,Xvel))*whi+(p4(i-1,jlop,k,Xvel)-p4(i-1,jlom,k,Xvel))*wlo);
                dvdy = (0.50*dxinv[1]) * ((p4(i  ,jhip,k,Yvel)-p4(i  ,jhim,k,Yvel))*whi+(p4(i-1,jlop,k,Yvel)-p4(i-1,jlom,k,Yvel))*wlo);

#endif
#if AMREX_SPACEDIM == 3

                const int khip = k + (int)f4(i,j,k).isConnected(0,0, 1);
                const int khim = k - (int)f4(i,j,k).isConnected(0,0,-1);
                const int klop = k + (int)f4(i-1,j,k).isConnected(0,0, 1);
                const int klom = k - (int)f4(i-1,j,k).isConnected(0,0,-1);
                whi = weights[khip-khim];
                wlo = weights[klop-klom];
                dudz = (0.5*dxinv[2]) * ((p4(i  ,j,khip,Xvel)-p4(i  ,j,khim,Xvel))*whi + (p4(i-1,j,klop,Xvel)-p4(i-1,j,klom,Xvel))*wlo);
                dwdz = (0.5*dxinv[2]) * ((p4(i  ,j,khip,Zvel)-p4(i  ,j,khim,Zvel))*whi + (p4(i-1,j,klop,Zvel)-p4(i-1,j,klom,Zvel))*wlo);

#endif
                divu = dudx + dvdy + dwdz;

                muf = 0.5*(d4(i,j,k,iMu)+d4(i-1,j,k,iMu));
                tauxx = muf*(2*dudx-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauxz = muf*(dudz+dwdx);

                fluxX(i,j,k,Xmom) -= tauxx;
                fluxX(i,j,k,Ymom) -= tauxy;
                fluxX(i,j,k,Zmom) -= tauxz;
                fluxX(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel) +  p4(i-1,j,k,Xvel))*tauxx
                                          +(p4(i,j,k,Yvel) + p4(i-1,j,k,Yvel))*tauxy
                                          +(p4(i,j,k,Zvel) + p4(i-1,j,k,Zvel))*tauxz
                                          +(d4(i,j,k,iKappa)+d4(i-1,j,k,iKappa))*dTdx);
            }
        }
    }

    // Y - direction
#if AMREX_SPACEDIM >= 2
    Real tauyy, tauyz, dTdy;
    Real dvdz=0, dwdy=0;
    Array4<Real> const& fluxY = fluxes[1].array();
    for     (int k = lo.z-AMREX_D_PICK(0,0,1); k <= hi.z+AMREX_D_PICK(0,0,1); ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x-1; i <= hi.x+1; ++i) {

                bool covered = f4(i,j,k).isCovered();
                bool connected = f4(i,j,k).isConnected(0,-1,0);
                bool other_covered = f4(i,j-1,k).isCovered();

                // only calculate fluxes for fluid cells and between cells that are connected
                if (covered || other_covered || !connected)
                    continue;

                dTdy = (d4(i,j,k,iTemp)-d4(i,j-1,k,iTemp))*dxinv[1];
                dudy = (p4(i,j,k,Xvel)-p4(i,j-1,k,Xvel))*dxinv[1];
                dvdy = (p4(i,j,k,Yvel)-p4(i,j-1,k,Yvel))*dxinv[1];
                dwdy = (p4(i,j,k,Zvel)-p4(i,j-1,k,Zvel))*dxinv[1];

                const int ihip = i + (int)f4(i,j,  k).isConnected( 1,0,0);
                const int ihim = i - (int)f4(i,j,  k).isConnected(-1,0,0);
                const int ilop = i + (int)f4(i,j-1,k).isConnected( 1,0,0);
                const int ilom = i - (int)f4(i,j-1,k).isConnected(-1,0,0);
                whi = weights[ihip-ihim];
                wlo = weights[ilop-ilom];

                dudx = (0.5*dxinv[0]) * ((p4(ihip,j  ,k,Xvel)-p4(ihim,j  ,k,Xvel))*whi + (p4(ilop,j-1,k,Xvel)-p4(ilom,j-1,k,Xvel))*wlo);
                dvdx = (0.5*dxinv[0]) * ((p4(ihip,j  ,k,Yvel)-p4(ihim,j  ,k,Yvel))*whi + (p4(ilop,j-1,k,Yvel)-p4(ilom,j-1,k,Yvel))*wlo);

#if AMREX_SPACEDIM == 3

                const int khip = k + (int)f4(i,j,  k).isConnected(0,0, 1);
                const int khim = k - (int)f4(i,j,  k).isConnected(0,0,-1);
                const int klop = k + (int)f4(i,j-1,k).isConnected(0,0, 1);
                const int klom = k - (int)f4(i,j-1,k).isConnected(0,0,-1);
                whi = weights[khip-khim];
                wlo = weights[klop-klom];

                dvdz = (0.5*dxinv[2]) * ((p4(i,j  ,khip,Yvel)-p4(i,j  ,khim,Yvel))*whi + (p4(i,j-1,klop,Yvel)-p4(i,j-1,klom,Yvel))*wlo);
                dwdz = (0.5*dxinv[2]) * ((p4(i,j  ,khip,Zvel)-p4(i,j  ,khim,Zvel))*whi + (p4(i,j-1,klop,Zvel)-p4(i,j-1,klom,Zvel))*wlo);

#endif
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,iMu)+d4(i,j-1,k,iMu));
                tauyy = muf*(2*dvdy-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauyz = muf*(dwdy+dvdz);

                fluxY(i,j,k,Xmom) -= tauxy;
                fluxY(i,j,k,Ymom) -= tauyy;
                fluxY(i,j,k,Zmom) -= tauyz;
                fluxY(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel)+p4(i,j-1,k,Xvel))*tauxy
                                          +(p4(i,j,k,Yvel)+p4(i,j-1,k,Yvel))*tauyy
                                          +(p4(i,j,k,Zvel)+p4(i,j-1,k,Zvel))*tauyz
                                          +(d4(i,j,k,iKappa) + d4(i,j-1,k,iKappa))*dTdy);

            }
        }
    }
#endif

    // Z - direction
#if AMREX_SPACEDIM == 3
    Real tauzz, dTdz;
    Array4<Real> const& fluxZ = fluxes[2].array();
    for     (int k = lo.z; k <= hi.z+1; ++k) {
        for   (int j = lo.y-1; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x-1; i <= hi.x+1; ++i) {

                bool covered = f4(i,j,k).isCovered();
                bool connected = f4(i,j,k).isConnected(0,0,-1);
                bool other_covered = f4(i,j,k-1).isCovered();

                // only calculate fluxes for fluid cells and between cells that are connected
                if (covered || other_covered || !connected)
                    continue;

                dTdz = (d4(i,j,k,iTemp)-d4(i,j,k-1,iTemp))*dxinv[2];
                dudz = (p4(i,j,k,Xvel)-p4(i,j,k-1,Xvel))*dxinv[2];
                dvdz = (p4(i,j,k,Yvel)-p4(i,j,k-1,Yvel))*dxinv[2];
                dwdz = (p4(i,j,k,Zvel)-p4(i,j,k-1,Zvel))*dxinv[2];

                const int ihip = i + (int)f4(i,j,k  ).isConnected( 1,0,0);
                const int ihim = i - (int)f4(i,j,k  ).isConnected(-1,0,0);
                const int ilop = i + (int)f4(i,j,k-1).isConnected( 1,0,0);
                const int ilom = i - (int)f4(i,j,k-1).isConnected(-1,0,0);
                whi = weights[ihip-ihim];
                wlo = weights[ilop-ilom];

                dudx = (0.5*dxinv[0]) * ((p4(ihip,j,k  ,Xvel)-p4(ihim,j,k  ,Xvel))*whi + (p4(ilop,j,k-1,Xvel)-p4(ilom,j,k-1,Xvel))*wlo);
                dwdx = (0.5*dxinv[0]) * ((p4(ihip,j,k  ,Zvel)-p4(ihim,j,k  ,Zvel))*whi + (p4(ilop,j,k-1,Zvel)-p4(ilom,j,k-1,Zvel))*wlo);

                const int jhip = j + (int)f4(i,j,k  ).isConnected(0 ,1,0);
                const int jhim = j - (int)f4(i,j,k  ).isConnected(0,-1,0);
                const int jlop = j + (int)f4(i,j,k-1).isConnected(0 ,1,0);
                const int jlom = j - (int)f4(i,j,k-1).isConnected(0,-1,0);
                whi = weights[jhip-jhim];
                wlo = weights[jlop-jlom];

                dvdy = (0.5*dxinv[1]) * ((p4(i,jhip,k  ,Yvel)-p4(i,jhim,k  ,Yvel))*whi + (p4(i,jlop,k-1,Yvel)-p4(i,jlom,k-1,Yvel))*wlo);
                dwdy = (0.5*dxinv[1]) * ((p4(i,jhip,k  ,Zvel)-p4(i,jhim,k  ,Zvel))*whi + (p4(i,jlop,k-1,Zvel)-p4(i,jlom,k-1,Zvel))*wlo);

                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,iMu)+d4(i,j,k-1,iMu));
                tauxz = muf*(dudz+dwdx);
                tauyz = muf*(dvdz+dwdy);
                tauzz = muf*(2.*dwdz-two_thirds*divu);

                fluxZ(i,j,k,Xmom) -= tauxz;
                fluxZ(i,j,k,Ymom) -= tauyz;
                fluxZ(i,j,k,Zmom) -= tauzz;
                fluxZ(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel)+p4(i,j,k-1,Xvel))*tauxz
                                          +(p4(i,j,k,Yvel)+p4(i,j,k-1,Yvel))*tauyz
                                          +(p4(i,j,k,Zvel)+p4(i,j,k-1,Zvel))*tauzz
                                          +(d4(i,j,k,iKappa) +d4(i,j,k-1,iKappa))*dTdz);

            }
        }
    }

#endif

}

#endif

void HydroState::calc_neutral_viscous_fluxes(const Box& box, Array<FArrayBox,
                                             AMREX_SPACEDIM> &fluxes,
                                             const Box& pbox,
                                             const FArrayBox &prim,
                                             EB_OPTIONAL(const EBCellFlagFab& flag,)
                                             const Real* dx) const
{

#ifdef AMREX_USE_EB
    if (flag.getType() != FabType::regular) {
        calc_neutral_viscous_fluxes_eb(box, fluxes, pbox, prim, flag, dx);
        return;
    }
#endif

    FArrayBox diff(pbox, Viscous::NUM_NEUTRAL_DIFF_COEFFS);
    calc_neutral_diffusion_terms(pbox,
                                 prim,
                                 diff
                             #ifdef AMREX_USE_EB
                                 ,flag
                             #endif
                                 );

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<Real, AMREX_SPACEDIM> dxinv;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        dxinv[d] = 1/dx[d];
    }

    Array4<const Real> const& p4 = prim.array();
    Array4<const Real> const& d4 = diff.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Real dudx=0, dudy=0, dudz=0, dvdx=0, dvdy=0, dwdx=0, dwdz=0, divu=0;
    const Real two_thirds = 2/3;

    Real tauxx, tauxy, tauxz, dTdx, muf;

    int iTemp = Viscous::NeutralTemp;
    int iMu = Viscous::NeutralMu;
    int iKappa = Viscous::NeutralKappa;

    Vector<int> prim_vel_id = get_prim_vector_idx();

    int Xvel = prim_vel_id[0] + 0;
    int Yvel = prim_vel_id[0] + 1;
    int Zvel = prim_vel_id[0] + 2;

    Vector<int> cons_vel_id = get_cons_vector_idx();

    int Xmom = cons_vel_id[0] + 0;
    int Ymom = cons_vel_id[0] + 1;
    int Zmom = cons_vel_id[0] + 2;

    Vector<int> nrg_id = get_nrg_idx();
    int Eden = nrg_id[0];

    Array4<Real> const& fluxX = fluxes[0].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x + 1; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                dTdx = (d4(i,j,k,iTemp) - d4(i-1,j,k,iTemp))*dxinv[0];

                dudx = (p4(i,j,k,Xvel) - p4(i-1,j,k,Xvel))*dxinv[0];
                dvdx = (p4(i,j,k,Yvel) - p4(i-1,j,k,Yvel))*dxinv[0];
                dwdx = (p4(i,j,k,Zvel) - p4(i-1,j,k,Zvel))*dxinv[0];

#if AMREX_SPACEDIM >= 2
                dudy = (p4(i,j+1,k,Xvel)+p4(i-1,j+1,k,Xvel)-p4(i,j-1,k,Xvel)-p4(i-1,j-1,k,Xvel))*(0.25*dxinv[1]);
                dvdy = (p4(i,j+1,k,Yvel)+p4(i-1,j+1,k,Yvel)-p4(i,j-1,k,Yvel)-p4(i-1,j-1,k,Yvel))*(0.25*dxinv[1]);
#endif
#if AMREX_SPACEDIM == 3
                dudz = (p4(i,j,k+1,Xvel)+p4(i-1,j,k+1,Xvel)-p4(i,j,k-1,Xvel)-p4(i-1,j,k-1,Xvel))*(0.25*dxinv[2]);
                dwdz = (p4(i,j,k+1,Zvel)+p4(i-1,j,k+1,Zvel)-p4(i,j,k-1,Zvel)-p4(i-1,j,k-1,Zvel))*(0.25*dxinv[2]);
#endif
                divu = dudx + dvdy + dwdz;

                muf = 0.5*(d4(i,j,k,iMu)+d4(i-1,j,k,iMu));
                tauxx = muf*(2*dudx-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauxz = muf*(dudz+dwdx);

                fluxX(i,j,k,Xmom) -= tauxx;
                fluxX(i,j,k,Ymom) -= tauxy;
                fluxX(i,j,k,Zmom) -= tauxz;
                fluxX(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel) +  p4(i-1,j,k,Xvel))*tauxx
                                          +(p4(i,j,k,Yvel) + p4(i-1,j,k,Yvel))*tauxy
                                          +(p4(i,j,k,Zvel) + p4(i-1,j,k,Zvel))*tauxz
                                          +(d4(i,j,k,iKappa)+d4(i-1,j,k,iKappa))*dTdx);

            }
        }
    }

#if AMREX_SPACEDIM >= 2
    Real tauyy, tauyz, dTdy;
    Real dvdz=0, dwdy=0;
    Array4<Real> const& fluxY = fluxes[1].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                dTdy = (d4(i,j,k,iTemp)-d4(i,j-1,k,iTemp))*dxinv[1];
                dudy = (p4(i,j,k,Xvel)-p4(i,j-1,k,Xvel))*dxinv[1];
                dvdy = (p4(i,j,k,Yvel)-p4(i,j-1,k,Yvel))*dxinv[1];
                dwdy = (p4(i,j,k,Zvel)-p4(i,j-1,k,Zvel))*dxinv[1];
                dudx = (p4(i+1,j,k,Xvel)+p4(i+1,j-1,k,Xvel)-p4(i-1,j,k,Xvel)-p4(i-1,j-1,k,Xvel))*(0.25*dxinv[0]);
                dvdx = (p4(i+1,j,k,Yvel)+p4(i+1,j-1,k,Yvel)-p4(i-1,j,k,Yvel)-p4(i-1,j-1,k,Yvel))*(0.25*dxinv[0]);
#if AMREX_SPACEDIM == 3
                dvdz = (p4(i,j,k+1,Yvel)+p4(i,j-1,k+1,Yvel)-p4(i,j,k-1,Yvel)-p4(i,j-1,k-1,Yvel))*(0.25*dxinv[2]);
                dwdz = (p4(i,j,k+1,Zvel)+p4(i,j-1,k+1,Zvel)-p4(i,j,k-1,Zvel)-p4(i,j-1,k-1,Zvel))*(0.25*dxinv[2]);
#endif
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,iMu)+d4(i,j-1,k,iMu));
                tauyy = muf*(2*dvdy-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauyz = muf*(dwdy+dvdz);

                fluxY(i,j,k,Xmom) -= tauxy;
                fluxY(i,j,k,Ymom) -= tauyy;
                fluxY(i,j,k,Zmom) -= tauyz;
                fluxY(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel)+p4(i,j-1,k,Xvel))*tauxy
                                          +(p4(i,j,k,Yvel)+p4(i,j-1,k,Yvel))*tauyy
                                          +(p4(i,j,k,Zvel)+p4(i,j-1,k,Zvel))*tauyz
                                          +(d4(i,j,k,iKappa) + d4(i,j-1,k,iKappa))*dTdy);

            }
        }
    }


#endif
#if AMREX_SPACEDIM == 3
    Real tauzz, dTdz;
    Array4<Real> const& fluxZ = fluxes[2].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                dTdz = (d4(i,j,k,iTemp)-d4(i,j,k-1,iTemp))*dxinv[2];
                dudz = (p4(i,j,k,Xvel)-p4(i,j,k-1,Xvel))*dxinv[2];
                dvdz = (p4(i,j,k,Yvel)-p4(i,j,k-1,Yvel))*dxinv[2];
                dwdz = (p4(i,j,k,Zvel)-p4(i,j,k-1,Zvel))*dxinv[2];
                dudx = (p4(i+1,j,k,Xvel)+p4(i+1,j,k-1,Xvel)-p4(i-1,j,k,Xvel)-p4(i-1,j,k-1,Xvel))*(0.25*dxinv[0]);
                dwdx = (p4(i+1,j,k,Zvel)+p4(i+1,j,k-1,Zvel)-p4(i-1,j,k,Zvel)-p4(i-1,j,k-1,Zvel))*(0.25*dxinv[0]);
                dvdy = (p4(i,j+1,k,Yvel)+p4(i,j+1,k-1,Yvel)-p4(i,j-1,k,Yvel)-p4(i,j-1,k-1,Yvel))*(0.25*dxinv[1]);
                dwdy = (p4(i,j+1,k,Zvel)+p4(i,j+1,k-1,Zvel)-p4(i,j-1,k,Zvel)-p4(i,j-1,k-1,Zvel))*(0.25*dxinv[1]);
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,iMu)+d4(i,j,k-1,iMu));
                tauxz = muf*(dudz+dwdx);
                tauyz = muf*(dvdz+dwdy);
                tauzz = muf*(2.*dwdz-two_thirds*divu);

                fluxZ(i,j,k,Xmom) -= tauxz;
                fluxZ(i,j,k,Ymom) -= tauyz;
                fluxZ(i,j,k,Zmom) -= tauzz;
                fluxZ(i,j,k,Eden) -= 0.5*((p4(i,j,k,Xvel)+p4(i,j,k-1,Xvel))*tauxz
                                          +(p4(i,j,k,Yvel)+p4(i,j,k-1,Yvel))*tauyz
                                          +(p4(i,j,k,Zvel)+p4(i,j,k-1,Zvel))*tauzz
                                          +(d4(i,j,k,iKappa) +d4(i,j,k-1,iKappa))*dTdz);

            }
        }
    }

#endif


    return;
}

// ====================================================================================
void HydroState::calc_ion_diffusion_terms(const Box& box,const Vector<FArrayBox>& prim,
                                          State& EMstate,Array4<const Real> const& prim_EM4,
                                          State& ELEstate,Array4<const Real> const& prim_ELE4,
                                          FArrayBox& diff
                                          EB_OPTIONAL(,const EBCellFlagFab& flag)
                                          ) const {

    return;
}

void HydroState::calc_ion_viscous_fluxes(const Box& box, 
                                         Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                         const Box& pbox, const Vector<FArrayBox>& prim,
                                         EB_OPTIONAL(const EBCellFlagFab& flag,)
                                         const Real* dx) const {

    return;
}

// ====================================================================================

void HydroState::calc_electron_diffusion_terms(const Box& box,const Vector<FArrayBox>& prim,
                                               State& EMstate,
                                               Array4<const Real> const& prim_EM4,
                                               State& IONstate,
                                               Array4<const Real> const& prim_ION4,
                                               FArrayBox& diff
                                               EB_OPTIONAL(,const EBCellFlagFab& flag)
                                               ) const {
    return;
}

void HydroState::calc_electron_viscous_fluxes(const Box& box, 
                                              Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                              const Box& pbox, const Vector<FArrayBox>& prim,
                                              EB_OPTIONAL(const EBCellFlagFab& flag,)
                                              const Real* dx) const {
    return;
}

// ====================================================================================

void HydroState::calc_charged_viscous_fluxes(int passed_idx,
                                             int ion_idx,
                                             int electron_idx,
                                             int em_idx,
                                             const Box& box,
                                             Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                             const Box& pbox,
                                             const Vector<FArrayBox>& prim,
                                             EB_OPTIONAL(const EBCellFlagFab& flag,)
                                             const Real* dx, FArrayBox& diff) const {
    return;
}


// ====================================================================================
void HydroState::write_info(nlohmann::json& js) const
{

    State::write_info(js);

    js["mass"] = mass;
    js["charge"] = charge;
    js["gamma"] = gamma;

    if (viscous) {

        auto& grp = js["viscosity"];

        int tp = viscous->get_type();
        grp["type"] = tp;

        const auto coeffs = viscous->get_refs();

        for (const auto& cf : coeffs) {
            grp[cf.first] = cf.second;
        }
    }
}

std::string HydroState::str() const
{
    std::stringstream msg;

    msg << State::str();

    if (viscous) {

        msg << "    viscosity : " << viscous->str() << "\n";

    }

    return msg.str();
}
