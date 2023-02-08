#include "MFP_hydro.H"
#include "MFP_lua.H"
#include "MFP.H"
#include "MFP_diagnostics.H"
#include "MFP_transforms.H"
#include "MFP_hydro_refine.H"

#include "Eigen"

std::string HydroState::multicomp_prim_name = "alpha";
std::string HydroState::multicomp_cons_name = "tracer";

Vector<std::string> HydroState::cons_names = {
    "rho",
    "x_mom",
    "y_mom",
    "z_mom",
    "nrg",
};

Vector<std::string> HydroState::prim_names = {
    "rho",
    "x_vel",
    "y_vel",
    "z_vel",
    "p",
    "T",
    "gamma",
    "cp",
};

Array<int,1> HydroState::cons_vector_idx = {+HydroDef::ConsIdx::Xmom};
Array<int,1> HydroState::prim_vector_idx = {+HydroDef::PrimIdx::Xvel};

std::map<std::string, int> HydroState::bc_names = {{"interior",  PhysBCType::interior},
                                                   {"inflow",    PhysBCType::inflow},
                                                   {"outflow",   PhysBCType::outflow},
                                                   {"symmetry",  PhysBCType::symmetry},
                                                   {"slipwall",  PhysBCType::slipwall},
                                                   {"noslipwall",PhysBCType::noslipwall}};

Vector<set_bc> HydroState::bc_set = {
    &set_scalar_bc,
    &set_x_vel_bc,
    &set_y_vel_bc,
    &set_z_vel_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc
};

std::string HydroState::tag = "hydro";
bool HydroState::registered = GetStateFactory().Register(HydroState::tag, StateBuilder<HydroState>);

HydroState::HydroState(){}

HydroState::HydroState(const sol::table& def)
{
    num_grow = 0;
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
}

HydroState::~HydroState(){}

void HydroState::set_gas()
{

    //
    // gas model
    //

    ClassFactory<HydroGas> gfact = GetHydroGasFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string gas_type = state_def["gas"]["type"].get_or<std::string>("");

    gas = gfact.Build(gas_type, state_def);

    if (!gas_type.empty() && !gas)
        Abort("Invalid gas option '"+gas_type+"'. Options are "+vec2str(gfact.getKeys()));
}

void HydroState::set_viscosity()
{

    //
    // viscous terms coefficients
    //

    ClassFactory<HydroViscous> vfact = GetHydroViscousFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string visc = state_def["viscosity"]["type"].get_or<std::string>("");

    viscous = vfact.Build(visc, state_def);

    if (!visc.empty() && !viscous)
        Abort("Invalid viscosity option '"+visc+"'. Options are "+vec2str(vfact.getKeys()));

    if (viscous) {
        set_num_grow(2); // needs at least two cells
    }
}

Real HydroState::init_from_number_density(std::map<std::string, Real> data)
{

    Vector<Real> Q(n_prim());

    Real nd = other_functions["nd"](data);
    if (gas->mass_const)
        return nd * gas->mass[0];

    for (int i=0; i < n_tracers; ++i) {
        Q[+HydroDef::PrimIdx::NUM + i] = functions[+HydroDef::PrimIdx::NUM + i](data);
    }

    return gas->get_density_from_number_density(nd, Q);
}

void HydroState::set_udf()
{
    using namespace std::placeholders;

    sol::state& lua = MFP::lua;

    sol::table state_def = lua["states"][name];

    // check if we have 'value' defined
    const sol::table value = state_def["value"].get_or(sol::table());

    if (!value.valid())
        Abort("State "+name+" does not have 'value' defined for initial conditions");

    bool check, success;

    // are there any alpha values?
    sol::object alpha = value[multicomp_prim_name].get_or(sol::object());
    const int i_start = +HydroDef::PrimIdx::NUM;
    if (alpha.valid()) {

        // turn it into a list if it isn't one already
        sol::table alpha_list;
        if (alpha.get_type() == sol::type::table) {
            alpha_list = alpha;
        } else {
            alpha_list = MFP::lua.create_table_with(1,alpha);
        }

        // handle the situation where we have only one component but still want to have a tracer
        if (gas->n_species() == 1) {
            n_tracers = alpha_list.size();

            // expand the components to suit
            gas->mass.resize(n_tracers + 1,gas->mass[0]);
            gas->charge.resize(n_tracers + 1,gas->charge[0]);
        }

        // are there enough tracer functions?
        if (alpha_list.size() != n_tracers) Abort("Incorrect number of 'alpha' defined for state '"+name+"' given the number of components, need "+num2str(n_tracers));

        const int i_stop = i_start + alpha_list.size();

        functions.resize(i_stop);

        for (int i = i_start; i<i_stop; ++i) {
            Optional3D1VFunction& v = functions[i];
            if (get_udf(alpha_list[i-i_start+1], v, 0.0)) {
                functions[i] = v;
            }
        }
    } else {
        const int i_stop = i_start + n_tracers;
        functions.resize(i_stop);
        for (int i = i_start; i<i_stop; ++i) {
            functions[i].set_value(0.0);
        }
    }

    // now for the primitives
    for (int i = 0; i<+HydroDef::PrimIdx::NUM; ++i) {

        std::string comp = prim_names[i];

        // is there a variable with this name?
        success = false;
        check = value[comp].valid();

        // doesn't exist, is there an alternative?
        if (!check) {

            // use number density instead of density
            if (i == +HydroDef::PrimIdx::Density) {

                check = value["nd"].valid();

                Optional3D1VFunction nd;
                success = get_udf(value["nd"], nd, 0.0);

                other_functions["nd"] = nd;

                Optional3D1VFunction rho;

                rho.set_func(std::bind(&HydroState::init_from_number_density, this, _1));

                functions[i] = rho;
            }

        }

        if (!success) {

            Optional3D1VFunction v;

            success = get_udf(value[comp], v, 0.0);

            functions[i] = v;
        }
    }

    return;
}

void HydroState::set_flux()
{

    if (!is_transported()) return;

    ClassFactory<HydroRiemannSolver> rfact = GetHydroRiemannSolverFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;
    state_def["n_prim"] = n_prim();
    state_def["n_cons"] = n_cons();
    state_def["n_tracer"] = n_tracers;

    std::string flux = state_def["flux"].get_or<std::string>("null");

    if (flux == "null")
        Abort("Flux option required for state '"+name+"'. Options are "+vec2str(rfact.getKeys()));

    flux_solver = rfact.Build(flux, state_def);

    if (!flux_solver)
        Abort("Invalid flux solver option '"+flux+"'. Options are "+vec2str(rfact.getKeys()));


    return;

}

void HydroState::set_shock_detector()
{

    ClassFactory<HydroShockDetector> sdfact = GetHydroShockDetectorFactory();

    sol::table sd_def = MFP::lua["states"][name]["shock_detector"].get_or(sol::table());

    if (!sd_def.valid()) return;

    sd_def["global_idx"] = global_idx;

    std::string sd_name = sd_def["name"].get_or<std::string>("");

    shock_detector = sdfact.Build(sd_name, sd_def);

    if (!sd_name.empty() && !shock_detector)
        Abort("Invalid shock_detector option '"+sd_name+"'. Options are "+vec2str(sdfact.getKeys()));
}

void HydroState::set_refinement()
{

    ClassFactory<Refinement> rfact = GetHydroRefinementFactory();

    sol::table r_def = MFP::lua["states"][name]["refinement"].get_or(sol::table());

    if (!r_def.valid()) return;

    r_def["global_idx"] = global_idx;

    std::string r_name = r_def["name"].get_or<std::string>("");

    refine = rfact.Build(r_name, r_def);

    if (!r_name.empty() && !refine)
        Abort("Invalid refinement option '"+r_name+"'. Options are "+vec2str(rfact.getKeys()));
}

void HydroState::init_from_lua()
{
    BL_PROFILE("HydroState::init_from_lua");

    EulerianState::init_from_lua();

    sol::state& lua = MFP::lua;

    const sol::table state_def = lua["states"][name];


    //
    // gas definition
    //
    set_gas();

    n_tracers = gas->n_tracers();

    //
    // viscous definition
    //
    set_viscosity();

    //
    // user defined functions
    //
    set_udf();

    for (int i = 0; i < n_tracers; ++i) {
        bc_set.push_back(&set_scalar_bc);
    }


    //
    // domain boundary conditions
    //

    const Vector<std::string> dir_name = {"x", "y", "z"};
    const Vector<std::string> side_name = {"lo", "hi"};
    const Vector<std::string>& hydro_var = prim_names;
    const int N = hydro_var.size();

    BoundaryState &bs = boundary_conditions;
    bs.phys_fill_bc.resize(+HydroDef::PrimIdx::NUM+n_tracers);

    const int j_start = +HydroDef::PrimIdx::NUM;
    const int j_stop = +HydroDef::PrimIdx::NUM+n_tracers;

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

            // now handle the tracers

            const sol::object alpha_funcs = state_def["bc"][dir_name[ax]][side_name[lh]][multicomp_prim_name].get_or(sol::object());

            if (alpha_funcs.valid()) {
                if (alpha_funcs.get_type() == sol::type::table) {
                    if (alpha_funcs.as<sol::table>().size() < n_tracers) {
                        Abort("Not enough boundary conditions specified for alpha in state '"+name+"'");
                    }
                }
            }

            for (int j=j_start; j<j_stop; ++j) {
                if (lh==0) {
                    bs.phys_fill_bc[j].setLo(ax,i_side_bc);
                } else {
                    bs.phys_fill_bc[j].setHi(ax,i_side_bc);
                }

                sol::object v;

                if (alpha_funcs.valid()) {
                    if (alpha_funcs.get_type() == sol::type::table) {
                        sol::table alpha_funcs_table = alpha_funcs.as<sol::table>();
                        v = alpha_funcs_table[j - j_start + 1];
                    } else {
                        v = alpha_funcs;
                    }
                }

                Optional3D1VFunction f = get_udf(v);
                // special case for inflow condition
                if (i_side_bc == PhysBCType::inflow && !f.is_valid()) {
                    Abort("Setting 'fill_hydro_bc = inflow' requires all primitive variables to be defined, 'alpha' is not defined");
                }

                bs.set(ax,get_multicomp_name(multicomp_prim_name, j-j_start),lh,f);
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
    // riemann solver
    //
    set_flux();

    //
    // shock detector
    //
    set_shock_detector();

    //
    // refinement
    //

    set_refinement();

}

void HydroState::variable_setup(Vector<int> periodic)
{

    boundary_conditions.fill_bc.resize(n_prim());

    for (int icomp=0; icomp < n_cons(); ++icomp) {
        set_bc s = bc_set[icomp]; // the function that sets the bc

        // make sure our periodicity isn't being overwritten
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            if (periodic[d]) {
                boundary_conditions.phys_fill_bc[icomp].setLo(d, PhysBCType::interior);
                boundary_conditions.phys_fill_bc[icomp].setHi(d, PhysBCType::interior);
            }
        }

        // grab the per component BCRec and apply the set bc function to it
        (*s)(boundary_conditions.fill_bc[icomp], boundary_conditions.phys_fill_bc[icomp]);
    }

    Vector<std::string> comp_names(n_cons());
    for (int icomp=0; icomp < +HydroDef::ConsIdx::NUM; ++icomp) {
        comp_names[icomp] = cons_names[icomp] + "-" + name;
    }

    for (int icomp=+HydroDef::ConsIdx::NUM; icomp < n_cons(); ++icomp) {
        comp_names[icomp] = multicomp_cons_name + "_" + num2str(icomp-+HydroDef::ConsIdx::NUM) + "-" + name;
    }

    int ng = num_grow;

#ifdef AMREX_USE_EB
    Interpolater* interp = &eb_cell_cons_interp;
#else
    Interpolater* interp = &cell_cons_interp;
#endif

    bool state_data_extrap = false;
    bool store_in_checkpoint = true;

    DescriptorList& desc_lst = MFP::get_desc_lst();

    data_idx = desc_lst.size();

    desc_lst.addDescriptor(data_idx, IndexType::TheCellType(),
                           StateDescriptor::Point, ng, n_cons(),
                           interp, state_data_extrap,
                           store_in_checkpoint);

    desc_lst.setComponent(
                data_idx, 0, comp_names, boundary_conditions.fill_bc,
                FillBC());


    if (MFP::verbosity >= 1) {
        Print() << str();
    }
}

void HydroState::init_data(MFP* mfp, const Real time)
{

    const Real* dx = mfp->Geom().CellSize();
    const Real* prob_lo = mfp->Geom().ProbLo();

    MultiFab& S_new = mfp->get_data(data_idx, time);

    Vector<Real> U(n_cons());
    Vector<Real> Q(n_prim());

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        FArrayBox& S_data = S_new[mfi];

        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);
        Array4<Real> const& S4 = S_data.array();

#ifdef AMREX_USE_EB
        FabArray<EBCellFlagFab>& flags = mfp->get_eb_data(global_idx).flags;
        Array4<const EBCellFlag> const& flag4 = flags.array(mfi);
#endif

        Real x, y, z;
        for     (int k = lo.z; k <= hi.z; ++k) {
            z = prob_lo[2] + (k + 0.5)*dx[2];
            for   (int j = lo.y; j <= hi.y; ++j) {
                y = prob_lo[1] + (j + 0.5)*dx[1];
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {
                    x = prob_lo[0] + (i + 0.5)*dx[0];

#ifdef AMREX_USE_EB
                    const EBCellFlag &cflag = flag4(i,j,k);

                    if (cflag.isCovered()) {
                        for (int n=0; n<n_cons(); ++n) {
                            S4(i,j,k,n) = 0.0;
                        }
                        continue;
                    }
#endif

                    // grab the primitive variables as defined by the user functions
                    for (int icomp=0; icomp<n_prim(); ++icomp) {
                        const auto& f = functions[icomp];

                        Q[icomp] = f(x, y, z);

                    }

                    gas->define_rho_p_T(Q);

                    gas->prim_valid(Q);

                    // convert primitive to conserved
                    gas->prim2cons(Q, U);

                    // copy into array
                    for (int n=0; n<n_cons(); ++n) {
                        S4(i,j,k,n) = U[n];
                    }
                }
            }
        }
    }
}



Vector<std::string> HydroState::get_plot_output_names() const
{
    Vector<std::string> out;
    out.insert(out.end(), cons_names.begin(), cons_names.end());
    out.insert(out.end(), prim_names.begin(), prim_names.end());
    out.push_back("charge");
    out.push_back("mass");
    out.push_back("gamma");
#ifdef AMREX_USE_EB
    out.push_back("vfrac");
#endif

    for (int i=0; i<n_tracers; ++i) {
        out.push_back(get_multicomp_name(multicomp_cons_name, i));
        out.push_back(get_multicomp_name(multicomp_prim_name, i));
    }

    return out;
}

void HydroState::get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string,FArrayBox>& out,
                                 Vector<std::string>& updated
                                 #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                                 #endif
                                 ) const
{
    BL_PROFILE("HydroState::get_state_values");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
    Array4<const Real> const& vf4 = vfrac.array();
#endif

    updated.resize(0);

    // check conserved variables
    std::map<std::string,int> cons_tags;
    for (int i=0; i<+HydroDef::ConsIdx::NUM; ++i) {
        const std::string s = cons_names[i];
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        cons_tags[var_name] = i;
        updated.push_back(var_name);
    }

    for (int i=0; i<n_tracers; ++i) {
        const std::string var_name = get_multicomp_name(multicomp_cons_name, i)+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        cons_tags[var_name] = i + +HydroDef::ConsIdx::NUM;
        updated.push_back(var_name);
    }



    // check primitive variables
    std::map<std::string,int> prim_tags;
    for (int i=0; i<+HydroDef::PrimIdx::NUM; ++i) {
        const std::string s = prim_names[i];
        if (s == cons_names[0]) continue;
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        prim_tags[var_name] = i;
        updated.push_back(var_name);
    }

    for (int i=0; i<n_tracers; ++i) {
        const std::string var_name = get_multicomp_name(multicomp_prim_name, i)+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        prim_tags[var_name] = i + +HydroDef::PrimIdx::NUM;
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

    //    const std::string gamma_name = "gamma-"+name;
    //    bool load_gamma = out.find(gamma_name) != out.end();
    //    if (load_gamma) other.push_back(gamma_name);

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
    Vector<Real> S(n_cons());
    Vector<Real> Q(n_prim());

    Array4<const Real> const& src4 = src.array();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (vf4(i,j,k) == 0.0) continue;
#endif

                for (int n=0; n<n_cons(); ++n) {
                    S[n] = src4(i,j,k,n);

                    //                    if (std::isnan(S[n])) {
                    //                        Abort();
                    //                    }

                }

                if (S[+HydroDef::ConsIdx::Density] < effective_zero) continue;


                if (load_charge) out4[charge_name](i,j,k) = gas->get_charge_from_cons(S);
                if (load_mass)   out4[mass_name](i,j,k)   = gas->get_mass_from_cons(S);
                //                if (load_gamma)  out4[gamma_name](i,j,k)  = get_gamma_from_cons(S);
#ifdef AMREX_USE_EB
                if (load_vfrac)  out4[vfrac_name](i,j,k)  = vf4(i,j,k);
#endif

                if (!cons_tags.empty()) {
                    for (const auto& var : cons_tags) {
                        out4[var.first](i,j,k) = S[var.second];
                    }
                }

                if (!prim_tags.empty()) {
                    gas->cons2prim(S, Q);

                    for (const auto& var : prim_tags) {
                        out4[var.first](i,j,k) = Q[var.second];
                    }
                }
            }
        }
    }


    return;
}

Real HydroState::get_allowed_time_step(MFP* mfp) const
{


    Real dt = EulerianState::get_allowed_time_step(mfp);

    // check for any viscous time step limitation
    if (viscous) {
        dt = std::min(dt, viscous->get_min_dt(mfp));
    }

    return dt;
}

void HydroState::calc_velocity(const Box& box,
                               FArrayBox& cons,
                               FArrayBox &prim
                               #ifdef AMREX_USE_EB
                               ,const FArrayBox& vfrac
                               #endif
                               ) const
{
    BL_PROFILE("HydroState::calc_velocity");


    prim.resize(box, AMREX_SPACEDIM);

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);
    Array4<Real> const& s4 = cons.array();
    Array4<Real> const& p4 = prim.array();

#ifdef AMREX_USE_EB
    Array4<const Real> const& vfrac4 = vfrac.array();
#endif

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (vfrac4(i,j,k) == 0.0) continue;
#endif

                const Real irho = 1/s4(i,j,k,+HydroDef::ConsIdx::Density);

                for (int n = 0; n<AMREX_SPACEDIM; ++n) {
                    p4(i,j,k,n) =  irho*s4(i,j,k,+HydroDef::ConsIdx::Xmom+n);
                }

            }
        }
    }

    return;
}

void HydroState::calc_reconstruction(const Box& box,
                                     FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                     #ifdef AMREX_USE_EB
                                     ,const EBCellFlagFab &flag
                                     ,const FArrayBox &vfrac
                                     #endif
                                     ) const
{
    BL_PROFILE("HydroState::calc_reconstruction");
    // if we don't want to apply extra limiting on the slopes (forced to 2nd order)
    // we can use the default reconstruction scheme

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

    Vector<Real> stencil(reconstructor->stencil_length);
    int offset = reconstructor->stencil_length/2;
    Array<int,3> stencil_index;
    Vector<Real> Q(n_prim()), cell_slope(n_prim());

    Real rho_lo, rho_hi;
    Real alpha_lo, alpha_hi;
    Real abs_phi, phi_scale, coeff_eps;
    Real gam_lo, gam_hi;

    Vector<Real> alphas_lo(n_tracers), alphas_hi(n_tracers);

    // make sure our arrays for putting lo and hi reconstructed values into
    // are the corect size
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        rlo[d].resize(box, n_prim());
        rhi[d].resize(box, n_prim());

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

                for (int n=0; n<n_prim(); ++n) {
                    Q[n] = src4(i,j,k,n);
                }

                gam4(i,j,k) = gas->get_gamma_from_prim(Q) - 1.0;

                src4(i,j,k,+HydroDef::PrimIdx::Prs) /= gam4(i,j,k);

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
                        for (int s=0; s<reconstructor->stencil_length; ++s) {
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
                    for (int n = 0; n<n_prim(); ++n) {

                        // fill in the stencil along dimension index
                        stencil_index.fill(0);
                        for (int s=0; s<reconstructor->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            stencil[s] = src4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], n);
                        }

                        // perform reconstruction
                        cell_slope[n] = reconstructor->get_slope(stencil);
                        Q[n] = stencil[offset];

                    }


                    // apply corrections to slopes
                    // J. Sci. Comput. (2014) 60:584-611
                    // Robust Finite Volume Schemes for Two-Fluid Plasma Equations

                    Real &rho     = Q[+HydroDef::PrimIdx::Density];
                    Real &phi_rho = cell_slope[+HydroDef::PrimIdx::Density];

                    Real &u     = Q[+HydroDef::PrimIdx::Xvel];
                    Real &phi_u = cell_slope[+HydroDef::PrimIdx::Xvel];

                    Real &v     = Q[+HydroDef::PrimIdx::Yvel];
                    Real &phi_v = cell_slope[+HydroDef::PrimIdx::Yvel];

                    Real &w     = Q[+HydroDef::PrimIdx::Zvel];
                    Real &phi_w = cell_slope[+HydroDef::PrimIdx::Zvel];

                    Real &eps     = Q[+HydroDef::PrimIdx::Prs];
                    Real &phi_eps = cell_slope[+HydroDef::PrimIdx::Prs];

                    // correct density slope
                    if (std::abs(phi_rho) > 2*rho) {
                        phi_rho = 2*sign(phi_rho, 0.0)*rho;
                    }

                    // get some face values
                    rho_lo = rho - 0.5*phi_rho;
                    rho_hi = rho + 0.5*phi_rho;

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
                    lo4(i,j,k,+HydroDef::PrimIdx::Density) = rho_lo;
                    hi4(i,j,k,+HydroDef::PrimIdx::Density) = rho_hi;

                    // x - velocity
                    lo4(i,j,k,+HydroDef::PrimIdx::Xvel) = u - 0.5*(rho_hi/rho)*phi_u;
                    hi4(i,j,k,+HydroDef::PrimIdx::Xvel) = u + 0.5*(rho_lo/rho)*phi_u;

                    // y - velocity
                    lo4(i,j,k,+HydroDef::PrimIdx::Yvel) = v - 0.5*(rho_hi/rho)*phi_v;
                    hi4(i,j,k,+HydroDef::PrimIdx::Yvel) = v + 0.5*(rho_lo/rho)*phi_v;

                    // z - velocity
                    lo4(i,j,k,+HydroDef::PrimIdx::Zvel) = w - 0.5*(rho_hi/rho)*phi_w;
                    hi4(i,j,k,+HydroDef::PrimIdx::Zvel) = w + 0.5*(rho_lo/rho)*phi_w;


                    for (int n=0; n<n_tracers; ++n) {
                        Real &alpha     = Q[+HydroDef::PrimIdx::NUM + n];
                        Real &phi_alpha = cell_slope[+HydroDef::PrimIdx::NUM + n];

                        alpha_lo = alpha - 0.5*phi_alpha;
                        alpha_hi = alpha + 0.5*phi_alpha;

                        // tracer
                        lo4(i,j,k,+HydroDef::PrimIdx::NUM + n) = alpha_lo;
                        hi4(i,j,k,+HydroDef::PrimIdx::NUM + n) = alpha_hi;

                        alphas_lo[n] = alpha_lo;
                        alphas_hi[n] = alpha_hi;
                    }

                    gam_lo = gas->get_gamma_from_prim(alphas_lo,0);
                    gam_hi = gas->get_gamma_from_prim(alphas_hi,0);

                    // epsilon -> pressure
                    lo4(i,j,k,+HydroDef::PrimIdx::Prs) = (eps - 0.5*phi_eps)*(gam_lo - 1.0);
                    hi4(i,j,k,+HydroDef::PrimIdx::Prs) = (eps + 0.5*phi_eps)*(gam_hi - 1.0);


                    // Temperature (calculate from pressure and density)
                    lo4(i,j,k,+HydroDef::PrimIdx::Temp) = lo4(i,j,k,+HydroDef::PrimIdx::Prs)/(rho_lo/gas->get_mass_from_prim(alphas_lo,0));
                    hi4(i,j,k,+HydroDef::PrimIdx::Temp) = hi4(i,j,k,+HydroDef::PrimIdx::Prs)/(rho_hi/gas->get_mass_from_prim(alphas_hi,0));

                    // gamma
                    lo4(i,j,k,+HydroDef::PrimIdx::Gamma) = gam_lo;
                    hi4(i,j,k,+HydroDef::PrimIdx::Gamma) = gam_hi;

                    // specific heat
                    lo4(i,j,k,+HydroDef::PrimIdx::SpHeat) = gas->get_cp_from_prim(alphas_lo,0);
                    hi4(i,j,k,+HydroDef::PrimIdx::SpHeat) = gas->get_cp_from_prim(alphas_hi,0);



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
                src4(i,j,k,+HydroDef::PrimIdx::Prs) *= gam4(i,j,k);



            }
        }
    }

    return;
}


void HydroState::calc_diffusion_terms(const FArrayBox& prim,
                                      FArrayBox& diff
                                      #ifdef AMREX_USE_EB
                                      ,const EBCellFlagFab& flag
                                      #endif
                                      ) const
{
    BL_PROFILE("HydroState::calc_neutral_diffusion_terms");
    const Dim3 lo = amrex::lbound(prim.box());
    const Dim3 hi = amrex::ubound(prim.box());

    Array4<const Real> const& prim4 = prim.array();
    Array4<Real> const& d4 = diff.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Vector<Real> Q(n_prim());
    Real T, mu, kappa;

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                for (int n=0; n<n_prim(); ++n) {
                    Q[n] = prim4(i,j,k,n);
                }

                viscous->get_coeffs(Q, T, mu, kappa);


                d4(i,j,k,+HydroViscous::CoeffIdx::Temp) = T;
                d4(i,j,k,+HydroViscous::CoeffIdx::Kappa) = kappa;
                d4(i,j,k,+HydroViscous::CoeffIdx::Mu) = mu;
            }
        }
    }

    return;
}


void HydroState::calc_viscous_fluxes(const Box& box, Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const FArrayBox &prim,
                                     #ifdef AMREX_USE_EB
                                     const EBCellFlagFab& flag,
                                     #endif
                                     const Real* dx) const
{
    BL_PROFILE("HydroState::calc_viscous_fluxes");
#ifdef AMREX_USE_EB
    if (flag.getType() != FabType::regular) {
        calc_viscous_fluxes_eb(box, fluxes, prim, flag, dx);
        return;
    }
#endif

    const Box pbox = prim.box();

    FArrayBox diff(pbox, +HydroViscous::CoeffIdx::NUM);
    calc_diffusion_terms(prim,
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

    Array4<Real> const& fluxX = fluxes[0].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x + 1; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered())
                    continue;
#endif

                dTdx = (d4(i,j,k,+HydroViscous::CoeffIdx::Temp) - d4(i-1,j,k,+HydroViscous::CoeffIdx::Temp))*dxinv[0];

                dudx = (p4(i,j,k,+HydroDef::PrimIdx::Xvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Xvel))*dxinv[0];
                dvdx = (p4(i,j,k,+HydroDef::PrimIdx::Yvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Yvel))*dxinv[0];
                dwdx = (p4(i,j,k,+HydroDef::PrimIdx::Zvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Zvel))*dxinv[0];

#if AMREX_SPACEDIM >= 2
                dudy = (p4(i,j+1,k,+HydroDef::PrimIdx::Xvel)+p4(i-1,j+1,k,+HydroDef::PrimIdx::Xvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Xvel)-p4(i-1,j-1,k,+HydroDef::PrimIdx::Xvel))*(0.25*dxinv[1]);
                dvdy = (p4(i,j+1,k,+HydroDef::PrimIdx::Yvel)+p4(i-1,j+1,k,+HydroDef::PrimIdx::Yvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Yvel)-p4(i-1,j-1,k,+HydroDef::PrimIdx::Yvel))*(0.25*dxinv[1]);
#endif
#if AMREX_SPACEDIM == 3
                dudz = (p4(i,j,k+1,Xvel)+p4(i-1,j,k+1,Xvel)-p4(i,j,k-1,Xvel)-p4(i-1,j,k-1,Xvel))*(0.25*dxinv[2]);
                dwdz = (p4(i,j,k+1,Zvel)+p4(i-1,j,k+1,Zvel)-p4(i,j,k-1,Zvel)-p4(i-1,j,k-1,Zvel))*(0.25*dxinv[2]);
#endif
                divu = dudx + dvdy + dwdz;

                muf = 0.5*(d4(i,j,k,+HydroViscous::CoeffIdx::Mu)+d4(i-1,j,k,+HydroViscous::CoeffIdx::Mu));
                tauxx = muf*(2*dudx-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauxz = muf*(dudz+dwdx);

                fluxX(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxx;
                fluxX(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauxy;
                fluxX(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauxz;
                fluxX(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel) +  p4(i-1,j,k,+HydroDef::PrimIdx::Xvel))*tauxx
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel) + p4(i-1,j,k,+HydroDef::PrimIdx::Yvel))*tauxy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel) + p4(i-1,j,k,+HydroDef::PrimIdx::Zvel))*tauxz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa)+d4(i-1,j,k,+HydroViscous::CoeffIdx::Kappa))*dTdx);

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

                dTdy = (d4(i,j,k,+HydroViscous::CoeffIdx::Temp)-d4(i,j-1,k,+HydroViscous::CoeffIdx::Temp))*dxinv[1];
                dudy = (p4(i,j,k,+HydroDef::PrimIdx::Xvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Xvel))*dxinv[1];
                dvdy = (p4(i,j,k,+HydroDef::PrimIdx::Yvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Yvel))*dxinv[1];
                dwdy = (p4(i,j,k,+HydroDef::PrimIdx::Zvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Zvel))*dxinv[1];
                dudx = (p4(i+1,j,k,+HydroDef::PrimIdx::Xvel)+p4(i+1,j-1,k,+HydroDef::PrimIdx::Xvel)-p4(i-1,j,k,+HydroDef::PrimIdx::Xvel)-p4(i-1,j-1,k,+HydroDef::PrimIdx::Xvel))*(0.25*dxinv[0]);
                dvdx = (p4(i+1,j,k,+HydroDef::PrimIdx::Yvel)+p4(i+1,j-1,k,+HydroDef::PrimIdx::Yvel)-p4(i-1,j,k,+HydroDef::PrimIdx::Yvel)-p4(i-1,j-1,k,+HydroDef::PrimIdx::Yvel))*(0.25*dxinv[0]);
#if AMREX_SPACEDIM == 3
                dvdz = (p4(i,j,k+1,Yvel)+p4(i,j-1,k+1,Yvel)-p4(i,j,k-1,Yvel)-p4(i,j-1,k-1,Yvel))*(0.25*dxinv[2]);
                dwdz = (p4(i,j,k+1,Zvel)+p4(i,j-1,k+1,Zvel)-p4(i,j,k-1,Zvel)-p4(i,j-1,k-1,Zvel))*(0.25*dxinv[2]);
#endif
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,+HydroViscous::CoeffIdx::Mu)+d4(i,j-1,k,+HydroViscous::CoeffIdx::Mu));
                tauyy = muf*(2*dvdy-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauyz = muf*(dwdy+dvdz);

                fluxY(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxy;
                fluxY(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauyy;
                fluxY(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauyz;
                fluxY(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Xvel))*tauxy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Yvel))*tauyy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Zvel))*tauyz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa) + d4(i,j-1,k,+HydroViscous::CoeffIdx::Kappa))*dTdy);

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

                fluxZ(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauyz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauzz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Xvel))*tauxz
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Yvel))*tauyz
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Zvel))*tauzz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa) +d4(i,j,k-1,+HydroViscous::CoeffIdx::Kappa))*dTdz);

            }
        }
    }

#endif


    return;
}

#ifdef AMREX_USE_EB
void HydroState::calc_viscous_fluxes_eb(const Box& box, Array<FArrayBox,
                                        AMREX_SPACEDIM> &fluxes,
                                        const FArrayBox &prim,
                                        #ifdef AMREX_USE_EB
                                        const EBCellFlagFab& flag,
                                        #endif
                                        const Real* dx) const
{
    BL_PROFILE("HydroState::calc_neutral_viscous_fluxes_eb");
    FArrayBox diff(prim.box(), +HydroViscous::CoeffIdx::NUM);
    calc_diffusion_terms(prim,
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

    Array4<const EBCellFlag> const& f4 = flag.array();

    Real dudx=0, dudy=0, dudz=0, dvdx=0, dvdy=0, dwdx=0, dwdz=0, divu=0;
    const Real two_thirds = 2/3;

    const Array<Real,3>  weights = {0.0, 1.0, 0.5};
    Real whi, wlo;

    Real tauxx, tauxy, tauxz, dTdx, muf;


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

                dTdx = (d4(i,j,k,+HydroViscous::CoeffIdx::Temp) - d4(i-1,j,k,+HydroViscous::CoeffIdx::Temp))*dxinv[0];

                dudx = (p4(i,j,k,+HydroDef::PrimIdx::Xvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Xvel))*dxinv[0];
                dvdx = (p4(i,j,k,+HydroDef::PrimIdx::Yvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Yvel))*dxinv[0];
                dwdx = (p4(i,j,k,+HydroDef::PrimIdx::Zvel) - p4(i-1,j,k,+HydroDef::PrimIdx::Zvel))*dxinv[0];

#if AMREX_SPACEDIM >= 2

                const int jhip = j + (int)f4(i,j,k).isConnected(0, 1,0);
                const int jhim = j - (int)f4(i,j,k).isConnected(0,-1,0);
                const int jlop = j + (int)f4(i-1,j,k).isConnected(0, 1,0);
                const int jlom = j - (int)f4(i-1,j,k).isConnected(0,-1,0);
                whi = weights[jhip-jhim];
                wlo = weights[jlop-jlom];
                dudy = (0.5*dxinv[1]) * ((p4(i  ,jhip,k,+HydroDef::PrimIdx::Xvel)-p4(i  ,jhim,k,+HydroDef::PrimIdx::Xvel))*whi+(p4(i-1,jlop,k,+HydroDef::PrimIdx::Xvel)-p4(i-1,jlom,k,+HydroDef::PrimIdx::Xvel))*wlo);
                dvdy = (0.50*dxinv[1]) * ((p4(i  ,jhip,k,+HydroDef::PrimIdx::Yvel)-p4(i  ,jhim,k,+HydroDef::PrimIdx::Yvel))*whi+(p4(i-1,jlop,k,+HydroDef::PrimIdx::Yvel)-p4(i-1,jlom,k,+HydroDef::PrimIdx::Yvel))*wlo);

#endif
#if AMREX_SPACEDIM == 3

                const int khip = k + (int)f4(i,j,k).isConnected(0,0, 1);
                const int khim = k - (int)f4(i,j,k).isConnected(0,0,-1);
                const int klop = k + (int)f4(i-1,j,k).isConnected(0,0, 1);
                const int klom = k - (int)f4(i-1,j,k).isConnected(0,0,-1);
                whi = weights[khip-khim];
                wlo = weights[klop-klom];
                dudz = (0.5*dxinv[2]) * ((p4(i  ,j,khip,+HydroDef::PrimIdx::Xvel)-p4(i  ,j,khim,+HydroDef::PrimIdx::Xvel))*whi + (p4(i-1,j,klop,+HydroDef::PrimIdx::Xvel)-p4(i-1,j,klom,+HydroDef::PrimIdx::Xvel))*wlo);
                dwdz = (0.5*dxinv[2]) * ((p4(i  ,j,khip,+HydroDef::PrimIdx::Zvel)-p4(i  ,j,khim,+HydroDef::PrimIdx::Zvel))*whi + (p4(i-1,j,klop,+HydroDef::PrimIdx::Zvel)-p4(i-1,j,klom,+HydroDef::PrimIdx::Zvel))*wlo);

#endif
                divu = dudx + dvdy + dwdz;

                muf = 0.5*(d4(i,j,k,+HydroViscous::CoeffIdx::Mu)+d4(i-1,j,k,+HydroViscous::CoeffIdx::Mu));
                tauxx = muf*(2*dudx-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauxz = muf*(dudz+dwdx);

                fluxX(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxx;
                fluxX(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauxy;
                fluxX(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauxz;
                fluxX(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel) +  p4(i-1,j,k,+HydroDef::PrimIdx::Xvel))*tauxx
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel) + p4(i-1,j,k,+HydroDef::PrimIdx::Yvel))*tauxy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel) + p4(i-1,j,k,+HydroDef::PrimIdx::Zvel))*tauxz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa)+d4(i-1,j,k,+HydroViscous::CoeffIdx::Kappa))*dTdx);
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

                dTdy = (d4(i,j,k,+HydroViscous::CoeffIdx::Temp)-d4(i,j-1,k,+HydroViscous::CoeffIdx::Temp))*dxinv[1];
                dudy = (p4(i,j,k,+HydroDef::PrimIdx::Xvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Xvel))*dxinv[1];
                dvdy = (p4(i,j,k,+HydroDef::PrimIdx::Yvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Yvel))*dxinv[1];
                dwdy = (p4(i,j,k,+HydroDef::PrimIdx::Zvel)-p4(i,j-1,k,+HydroDef::PrimIdx::Zvel))*dxinv[1];

                const int ihip = i + (int)f4(i,j,  k).isConnected( 1,0,0);
                const int ihim = i - (int)f4(i,j,  k).isConnected(-1,0,0);
                const int ilop = i + (int)f4(i,j-1,k).isConnected( 1,0,0);
                const int ilom = i - (int)f4(i,j-1,k).isConnected(-1,0,0);
                whi = weights[ihip-ihim];
                wlo = weights[ilop-ilom];

                dudx = (0.5*dxinv[0]) * ((p4(ihip,j  ,k,+HydroDef::PrimIdx::Xvel)-p4(ihim,j  ,k,+HydroDef::PrimIdx::Xvel))*whi + (p4(ilop,j-1,k,+HydroDef::PrimIdx::Xvel)-p4(ilom,j-1,k,+HydroDef::PrimIdx::Xvel))*wlo);
                dvdx = (0.5*dxinv[0]) * ((p4(ihip,j  ,k,+HydroDef::PrimIdx::Yvel)-p4(ihim,j  ,k,+HydroDef::PrimIdx::Yvel))*whi + (p4(ilop,j-1,k,+HydroDef::PrimIdx::Yvel)-p4(ilom,j-1,k,+HydroDef::PrimIdx::Yvel))*wlo);

#if AMREX_SPACEDIM == 3

                const int khip = k + (int)f4(i,j,  k).isConnected(0,0, 1);
                const int khim = k - (int)f4(i,j,  k).isConnected(0,0,-1);
                const int klop = k + (int)f4(i,j-1,k).isConnected(0,0, 1);
                const int klom = k - (int)f4(i,j-1,k).isConnected(0,0,-1);
                whi = weights[khip-khim];
                wlo = weights[klop-klom];

                dvdz = (0.5*dxinv[2]) * ((p4(i,j  ,khip,+HydroDef::PrimIdx::Yvel)-p4(i,j  ,khim,+HydroDef::PrimIdx::Yvel))*whi + (p4(i,j-1,klop,+HydroDef::PrimIdx::Yvel)-p4(i,j-1,klom,+HydroDef::PrimIdx::Yvel))*wlo);
                dwdz = (0.5*dxinv[2]) * ((p4(i,j  ,khip,+HydroDef::PrimIdx::Zvel)-p4(i,j  ,khim,+HydroDef::PrimIdx::Zvel))*whi + (p4(i,j-1,klop,+HydroDef::PrimIdx::Zvel)-p4(i,j-1,klom,+HydroDef::PrimIdx::Zvel))*wlo);

#endif
                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,+HydroViscous::CoeffIdx::Mu)+d4(i,j-1,k,+HydroViscous::CoeffIdx::Mu));
                tauyy = muf*(2*dvdy-two_thirds*divu);
                tauxy = muf*(dudy+dvdx);
                tauyz = muf*(dwdy+dvdz);

                fluxY(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxy;
                fluxY(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauyy;
                fluxY(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauyz;
                fluxY(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Xvel))*tauxy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Yvel))*tauyy
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel)+p4(i,j-1,k,+HydroDef::PrimIdx::Zvel))*tauyz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa) + d4(i,j-1,k,+HydroViscous::CoeffIdx::Kappa))*dTdy);

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

                dTdz = (d4(i,j,k,+HydroViscous::CoeffIdx::Temp)-d4(i,j,k-1,+HydroViscous::CoeffIdx::Temp))*dxinv[2];
                dudz = (p4(i,j,k,+HydroDef::PrimIdx::Xvel)-p4(i,j,k-1,+HydroDef::PrimIdx::Xvel))*dxinv[2];
                dvdz = (p4(i,j,k,+HydroDef::PrimIdx::Yvel)-p4(i,j,k-1,+HydroDef::PrimIdx::Yvel))*dxinv[2];
                dwdz = (p4(i,j,k,+HydroDef::PrimIdx::Zvel)-p4(i,j,k-1,+HydroDef::PrimIdx::Zvel))*dxinv[2];

                const int ihip = i + (int)f4(i,j,k  ).isConnected( 1,0,0);
                const int ihim = i - (int)f4(i,j,k  ).isConnected(-1,0,0);
                const int ilop = i + (int)f4(i,j,k-1).isConnected( 1,0,0);
                const int ilom = i - (int)f4(i,j,k-1).isConnected(-1,0,0);
                whi = weights[ihip-ihim];
                wlo = weights[ilop-ilom];

                dudx = (0.5*dxinv[0]) * ((p4(ihip,j,k  ,+HydroDef::PrimIdx::Xvel)-p4(ihim,j,k  ,+HydroDef::PrimIdx::Xvel))*whi + (p4(ilop,j,k-1,+HydroDef::PrimIdx::Xvel)-p4(ilom,j,k-1,+HydroDef::PrimIdx::Xvel))*wlo);
                dwdx = (0.5*dxinv[0]) * ((p4(ihip,j,k  ,+HydroDef::PrimIdx::Zvel)-p4(ihim,j,k  ,+HydroDef::PrimIdx::Zvel))*whi + (p4(ilop,j,k-1,+HydroDef::PrimIdx::Zvel)-p4(ilom,j,k-1,+HydroDef::PrimIdx::Zvel))*wlo);

                const int jhip = j + (int)f4(i,j,k  ).isConnected(0 ,1,0);
                const int jhim = j - (int)f4(i,j,k  ).isConnected(0,-1,0);
                const int jlop = j + (int)f4(i,j,k-1).isConnected(0 ,1,0);
                const int jlom = j - (int)f4(i,j,k-1).isConnected(0,-1,0);
                whi = weights[jhip-jhim];
                wlo = weights[jlop-jlom];

                dvdy = (0.5*dxinv[1]) * ((p4(i,jhip,k  ,+HydroDef::PrimIdx::Yvel)-p4(i,jhim,k  ,+HydroDef::PrimIdx::Yvel))*whi + (p4(i,jlop,k-1,+HydroDef::PrimIdx::Yvel)-p4(i,jlom,k-1,+HydroDef::PrimIdx::Yvel))*wlo);
                dwdy = (0.5*dxinv[1]) * ((p4(i,jhip,k  ,+HydroDef::PrimIdx::Zvel)-p4(i,jhim,k  ,+HydroDef::PrimIdx::Zvel))*whi + (p4(i,jlop,k-1,+HydroDef::PrimIdx::Zvel)-p4(i,jlom,k-1,+HydroDef::PrimIdx::Zvel))*wlo);

                divu = dudx + dvdy + dwdz;
                muf = 0.5*(d4(i,j,k,+HydroViscous::CoeffIdx::Mu)+d4(i,j,k-1,+HydroViscous::CoeffIdx::Mu));
                tauxz = muf*(dudz+dwdx);
                tauyz = muf*(dvdz+dwdy);
                tauzz = muf*(2.*dwdz-two_thirds*divu);

                fluxZ(i,j,k,+HydroDef::ConsIdx::Xmom) -= tauxz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Ymom) -= tauyz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Zmom) -= tauzz;
                fluxZ(i,j,k,+HydroDef::ConsIdx::Eden) -= 0.5*((p4(i,j,k,+HydroDef::PrimIdx::Xvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Xvel))*tauxz
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Yvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Yvel))*tauyz
                                                              +(p4(i,j,k,+HydroDef::PrimIdx::Zvel)+p4(i,j,k-1,+HydroDef::PrimIdx::Zvel))*tauzz
                                                              +(d4(i,j,k,+HydroViscous::CoeffIdx::Kappa) +d4(i,j,k-1,+HydroViscous::CoeffIdx::Kappa))*dTdz);

            }
        }
    }

#endif

}


#endif

void HydroState::calc_current_and_charge(const Box& box,
                                         const FArrayBox& cons,
                                         FArrayBox* cd,
                                         FArrayBox* J
                                         #ifdef AMREX_USE_EB
                                         ,const FArrayBox& vfrac
                                         #endif
                                         ) const
{
    BL_PROFILE("HydroState::calc_current_and_charge");

    Vector<Real> U(n_cons());

    const bool get_current = (J != nullptr);
    const bool get_charge = (cd != nullptr);

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);
    Array4<const Real> const& s4 = cons.array();


    Array4<Real> J4, cd4;

    if (get_charge) cd4 = cd->array();
    if (get_current) J4 = J->array();

#ifdef AMREX_USE_EB
    Array4<const Real> const& vfrac4 = vfrac.array();

    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);
#endif

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (vfrac4(i,j,k) == 0.0) {
                    continue;
                } else {
#endif
                    // grab the conserved variables
                    for (int n=0; n<n_cons(); ++n) {
                        U[n] = s4(i,j,k,n);
                    }
#ifdef AMREX_USE_EB
                }
#endif

                // calculate the mass and charge
                const Real m = gas->get_mass_from_cons(U);
                const Real q = gas->get_charge_from_cons(U);

                // calculate the charge density and current
                if (get_charge) {
                    cd4(i,j,k) += U[+HydroDef::ConsIdx::Density]*q/m;
                }

                if (get_current) {
                    J4(i,j,k,0) += U[+HydroDef::ConsIdx::Xmom]*q/m;
                    J4(i,j,k,1) += U[+HydroDef::ConsIdx::Ymom]*q/m;
                    J4(i,j,k,2) += U[+HydroDef::ConsIdx::Zmom]*q/m;
                }
            }
        }
    }

    return;
}

void HydroState::write_info(nlohmann::json &js) const
{

    EulerianState::write_info(js);

    // write out stuff that is common to all states

    js["type"] = tag;
    js["type_idx"] = +get_type();

    gas->write_info(js["gas"]);

    if (viscous) {

        auto& grp = js["viscosity"];

        grp["type"] = viscous->get_tag();

        const auto coeffs = viscous->get_refs();

        for (const auto& cf : coeffs) {
            grp[cf.first] = cf.second;
        }
    }

}
