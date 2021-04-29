#include "MFP_global.H"

#include "MFP_source.H"
#include "MFP_ode_solver.h"
#include "MFP_utility.H"
#include "MFP_eb_sdf.H"
#include "MFP_read_geom.h"
#include "MFP_polyspline.H"

using namespace amrex;

sol::state GlobalData::lua;
std::string GlobalData::lua_script;  // lua script that defines our problem

int GlobalData::Cost_Idx = -1;
int GlobalData::Shock_Idx = -1;
int GlobalData::num_states = 0;
int GlobalData::num_solve_state = 0;
int GlobalData::full_state_size = 0;

int GlobalData::zero_dimensional = false;

Array<bool, AMREX_SPACEDIM> GlobalData::periodic = {AMREX_D_DECL(false, false, false)};

Vector<std::string> GlobalData::state_names;
std::map<std::string, int> GlobalData::state_index;
Vector<std::string> GlobalData::state_tags;
bool GlobalData::detect_shocks = false;
bool GlobalData::plot_shock_detector = false;

int GlobalData::num_fluid = 0; // how many fluid type states are there

int GlobalData::num_shock_detector = 0;

int GlobalData::do_CTU = 1;

int GlobalData::do_face_src = 1;

Real GlobalData::force_dt = 0;

// Parameters
int GlobalData::verbose = 0;
int GlobalData::linear_solver_verbosity = 0;

IntVect GlobalData::tile_size;
Real GlobalData::cfl = 0.25;

int GlobalData::refine_cutcells = 1;

Vector<amrex::RealBox> GlobalData::refine_boxes;
Vector<RefineBoxType> GlobalData::refine_box_type;
bool GlobalData::only_refine_in_box = false;
bool GlobalData::derefine_box = false;

#ifdef AMREX_PARTICLES
int GlobalData::do_tracer_particles = 0;
int GlobalData::particle_verbose = 0;
#endif

//int GlobalData::do_tracer_particles = 0;
//int GlobalData::particle_verbose = 0;
//std::string GlobalData::particle_init_file;

Real GlobalData::x_ref = 1.0;
Real GlobalData::n_ref = 1.0;
Real GlobalData::m_ref = 1.0;
Real GlobalData::rho_ref = 1.0;
Real GlobalData::T_ref = 1.0;
Real GlobalData::u_ref = 1.0;
Real GlobalData::n0 = 1.0;
Real GlobalData::lightspeed = 1.0;
Real GlobalData::beta = 1.0;
Real GlobalData::skin_depth = 1.0;
Real GlobalData::Larmor = 1.0;
Real GlobalData::Debye = 1.0;

bool GlobalData::plasma_params_set = false;

Real GlobalData::effective_zero = 1e-14;

int GlobalData::msg = 0;

Vector<std::unique_ptr<State>> GlobalData::states;
Vector<std::unique_ptr<ODESystem>> GlobalData::ode_source_terms;

#ifdef AMREX_USE_EB
Vector<DataEB> GlobalData::eb_def;
Vector<Vector<EBData>> GlobalData::eb_data;
#endif

std::map<std::string,Array<int,AMREX_SPACEDIM+1>> GlobalData::plot_variables;
Vector<std::pair<std::string, Optional3D1VFunction>> GlobalData::plot_functions;

// pointer to state initialization routine
typedef void (*state_config)(const std::string &name, const int idx);

// check that the supplied index gives a state that is valid for the thing that will be
// applied to it
typedef bool (*state_valid)(const int idx);


//------------
// Global data


GlobalData::GlobalData()
{
}

GlobalData::~GlobalData()
{
}

void GlobalData::set_num_levels(int n)
{
#ifdef AMREX_USE_EB
    eb_data.resize(n);
#endif
}

void GlobalData::set_lua_script(const std::string &script)
{
    BL_PROFILE("GlobalData::set_lua_script");
    lua_script =
        #include "default.lua"
            ;

    lua_script += script;

    // insert call to preprocessor defined in default.lua
    lua_script += "\npreprocess()\n";
}

void GlobalData::read_config(const Vector<int> &is_periodic, const bool plot_output)
{
    BL_PROFILE("GlobalData::read_config");
    // periodicity
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        periodic[d] = (bool) is_periodic[d];
    }

    lua.open_libraries(sol::lib::base,
                       sol::lib::package,
                       sol::lib::math,
                       sol::lib::table,
                       sol::lib::string,
                       sol::lib::io);

#ifdef AMREX_USE_EB
    // register geometry stuff
    LuaSplineIF::register_with_lua(lua);
    ReadSTL::register_with_lua(lua);
#if AMREX_SPACEDIM == 2
    PolySpline::register_with_lua(lua);
    ReadSVG::register_with_lua(lua);
#endif
#endif

    // parse the input script
    lua.script(lua_script);

    // just in case there is any io going on wait until all have finished
    // this is a bit of a hack as we currently do not handle i/o collisions
    // between different lua threads
    ParallelDescriptor::Barrier();

    verbose = lua["verbosity"];

    linear_solver_verbosity = lua["linear_solver_verbosity"];

    //
    // reference variables
    //

    update_ref();

    //
    // generate states
    //

    // what states do we have?
    sol::table names = lua.script("return get_sorted_keys(states)");

    resize(names.size());

    PhysicsFactory<State> sfact = GetStateFactory();

    // get a list of all the tags
    for (const auto& S : sfact.getRegistered()) {
        state_tags.push_back(S.first);
    }


    int global_idx = 0;
    int type_idx;
    for (auto& item : names) {
        std::string name = item.second.as<std::string>();

        sol::table state_def = lua["states"][name];

        state_def["name"] = name;
        state_def["global_idx"] = global_idx;

        std::string state_type = state_def["type"].get<std::string>();

        std::unique_ptr<State> istate = sfact.Build(state_type, state_def);

        if (!istate)
            Abort("Failed to read state "+name+", must be one of "+vec2str(state_tags));

        states.push_back(std::move(istate));

        state_names[global_idx] = name;
        state_index[name] = global_idx;

        full_state_size += states[global_idx]->n_cons();

        global_idx++;
    }

    //
    // initialize states
    //

    for (auto &istate : states) {
        istate->init_from_lua();
    }

    //
    // forced refinement
    //

    sol::table ref_boxes = lua["refine_boxes"];
    Array<Real,AMREX_SPACEDIM> boxlo, boxhi;
    for (const auto& ubox : ref_boxes) {
        const sol::table &pair = ubox.second;
        const sol::table &lo = pair[1];
        const sol::table &hi = pair[2];
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            boxlo[i] = lo[i + 1];
            boxhi[i] = hi[i + 1];
        }
        refine_boxes.push_back(RealBox(boxlo, boxhi));

        std::string type = pair.get_or<std::string>("type","force_refine");
        if (type == "force_refine") {
            refine_box_type.push_back(RefineBoxType::ForceRefine);
        } else if (type == "no_refine") {
            refine_box_type.push_back(RefineBoxType::NoRefine);
            derefine_box = !only_refine_in_box && true;
        } else if (type == "only_refine") {
            refine_box_type.push_back(RefineBoxType::OnlyRefine);
            only_refine_in_box = true;
        } else {
            Abort("Refine box has undefined type '"+type+"', valid types are 'force_refine', 'no_refine', 'only_refine'");
        }
    }

    num_shock_detector = 0;
    for (global_idx=0; global_idx<num_solve_state; ++global_idx) {
        State &istate = get_state(global_idx);
        if (istate.shock_idx > -1) {
            istate.shock_idx = num_shock_detector; ++num_shock_detector;
        }
    }

    //
    // miscellaneous
    //

    zero_dimensional = lua["zero_dimensional"];

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        tile_size[i] = lua["tile_size"][i + 1];
    }

    do_CTU = lua["do_CTU"];

    do_face_src = lua["do_face_sources"];

    force_dt = lua["force_dt"];

    effective_zero = lua["effective_zero"];

    cfl = lua["cfl"];

    //
    // geometry
    //

#ifdef AMREX_USE_EB
    refine_cutcells = lua["refine_cutcells"];

    // get a sorted list of keys so that we have consistent tagging of boundaries
    sol::table eb_keys = lua.script("return get_sorted_keys(embedded_boundaries)");
    sol::table eb_list = lua["embedded_boundaries"];

    std::size_t i = 0;
    for (const auto& eb_key : eb_keys) {
        //        const sol::object &eb_name = eb_item.first;
        const sol::table &eb_desc = eb_list[eb_key.second.as<std::string>()];

        DataEB eb_dat;

        // get the geometry function
        sol::object eb_geom = eb_desc["geom"];

        if (eb_geom.get_type() == sol::type::function) {
            eb_dat.geom_func = eb_geom.as<lua_IF>();
        } else {
            Abort("EB geometry must be defined as a signed distance function");
        }

        // inside/outside
        eb_dat.inside = eb_desc.get_or("inside", 1) > 0;

        // how to insert it if there are multiples
        // default to complement
        const std::string boolean_op = eb_desc.get_or<std::string>("boolean_operation", "or");
        if (boolean_op == "or") {
            eb_dat.insertion_type = 1;
        } else if (boolean_op == "and"){
            eb_dat.insertion_type = -1;
        } else {
            Abort("Unknown option for 'boolean_operation', try 'and' or 'or'");
        }

        eb_dat.index = i;

        // get the names of the states that this geometry is interacting with
        // and the type of boundary condition that state uses to interact with
        // the embedded boundary
        std::size_t j = 0;
        const sol::table &bcs = eb_desc["bcs"];
        for (const auto& bc : bcs) {
            const sol::object &state = bc.first;
            const sol::object &bc_def = bc.second;
            std::string state_name = state.as<std::string>();
            if (state_name != "func") {
                State &istate = get_state(state_name);
                eb_dat.states.push_back({istate.global_idx,istate.eb_bcs.size()});

                // define bc
                istate.set_eb_bc(bc_def);
            }
            ++j;
        }

        eb_def.push_back(eb_dat);

        ++i;
    }
#endif

    // what source term collections have been specified
    sol::table ode_systems = lua["sources"];

    for (const auto& lua_ode : ode_systems) {

        const std::string system_name = lua_ode.first.as<std::string>();
        const sol::table system_def = lua_ode.second;

        // create a new collection
        ode_source_terms.push_back(std::move(std::unique_ptr<ODESystem>(new ODESystem())));
        ODESystem &local_ode = *ode_source_terms.back();

        local_ode.name = system_name;

        //---------------------------------------------------------------------
        // get the solver

        PhysicsFactory<SolveODE> solvefact = GetSolveODEFactory();

        std::string solve_name = system_def["solver"].get<std::string>();

        std::unique_ptr<SolveODE> solver = solvefact.Build(solve_name, system_def);

        if (!solver)
            Abort("Source solver option not recognised, must be one of "+vec2str(solvefact.getKeys()));

        local_ode.set_solver(std::move(solver));



        //---------------------------------------------------------------------

        // get all of the source terms in this collection
        sol::table sources = system_def["sources"];

        for (auto& srcs : sources) {

            const std::string tag = srcs.first.as<std::string>();
            sol::table source_def = srcs.second;

            PhysicsFactory<SourceTerm> sfact = GetSourceTermFactory();

            source_def["name"] = tag;
            source_def["solver"] = +local_ode.solver->get_type();

            std::string key = source_def["type"].get<std::string>();

            if (key.empty())
                Abort("Invalid source type provided for source '"+tag+"', must be one of "+vec2str(sfact.getKeys()));

            std::unique_ptr<SourceTerm> src = sfact.Build(key, source_def);

            if (!src)
                Abort("Source option '"+key+"' not recognised, must be one of "+vec2str(sfact.getKeys()));

            // let the state know what sources it is a part of
            for (const auto& oi : src->offsets) {
                State& istate = get_state(oi.global);
                istate.associated_sources.push_back({ode_source_terms.size()-1, local_ode.sources.size()});
            }

            local_ode.add_source(std::move(src));
        }

        if (local_ode.empty()) {
            Abort("Error: ODE system '"+system_name+"' inactive due to no active/valid sources");
        }

        if (verbose >= 1) {
            Print() << local_ode.print();
        }

    }

    // update anything in the states that requires source terms to be defined
    for (auto &istate : states) {
        istate->post_init_from_lua();
    }

    // plot functions
    const sol::table plot = lua["plot"];

    const Array<std::string, AMREX_SPACEDIM> grad = {AMREX_D_DECL("-dx","-dy","-dz")};

    const sol::table plot_vars = plot["variables"];


    // first get all the variables we want in our output
    for (const auto& key_value_pair : plot_vars ) {
         sol::object value = key_value_pair.second;
         std::string name = value.as<std::string>();
         plot_variables[name] = {0, AMREX_D_DECL(0,0,0)};
    }

    // now check if we have any request for gradients
    for (const auto& var : plot_variables ) {
        std::string name = var.first;
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            size_t found = name.find(grad[i]);
            if (found != std::string::npos) {
                const std::string base_name = name.erase(found,found+grad[i].size());
                if (plot_variables.find(base_name) == plot_variables.end())
                    plot_variables[base_name] = {-1, AMREX_D_DECL(0,0,0)};
                plot_variables[base_name][i+1] = 1;
                break;
            }
        }
    }

    // get everything???
    if (plot_variables.find("all") != plot_variables.end()) {

        plot_variables["cost"][0] = 0;

        for (int idx=0; idx<num_solve_state; ++idx) {
            State &istate = get_state(idx);
            for (auto& name : istate.get_prim_names()) {
                plot_variables[name+"-"+istate.name][0] = 0;
            }

#ifdef AMREX_USE_EB

            plot_variables["vfrac-"+istate.name][0] = 0;
#endif

            if (istate.shock_idx > -1) {
                plot_shock_detector = true;
                plot_variables["shock-"+istate.name][0] = 0;
            }
        }
        plot_variables.erase("all");
    }

    int count = 0;
    for (auto& var : plot_variables) {
        if (var.second[0] > -1) {
            var.second[0] = count;
            count++;
        }
    }

//    for (const auto& var : plot_variables) {
//        Print() << var.first << " : " << var.second[0] << std::endl;
//    }

    const sol::table plot_funcs = plot["functions"];
    // need to get sorted function names to ensure consistency across processors
    const sol::table func_names = lua.script("return get_sorted_keys(plot['functions'])");
    for (const auto& key_value_pair : func_names ) {
         const std::string name = key_value_pair.second.as<std::string>();
         plot_functions.push_back(std::make_pair(name, get_udf(plot_funcs[name])));
    }
}

void GlobalData::update_ref()
{
    BL_PROFILE("GlobalData::update_ref");
    x_ref = lua["ref_length"];
    rho_ref = lua["ref_density"];
    m_ref = lua["ref_mass"];
    T_ref = lua["ref_temp"];
    Real c_ref = lua["ref_lightspeed"].get_or(299792458.0);

    Real c = lua["lightspeed"];
    beta = lua["beta"];
    skin_depth = lua["skin_depth"];
    Larmor = lua["Larmor"];
    Debye = lua["Debye"];

    n_ref = rho_ref/m_ref;

    n0 = n_ref*x_ref*x_ref*x_ref;

    const Real kB = 1.38064852e-23; // m^2.kg.s^-2.K^-1
    const Real mu0= 1.25663706e-6 ; // m kg s-2 A-2
    const Real ep0= 8.85418782e-12; // m-3 kg-1 s4 A2
    // have choice of T or c
    if ((T_ref > 0) && (c > 0)) {
        Abort("ref_temp and lightspeed both set, only set one");
    } else if (T_ref > 0) {
        //Print() << "T_ref >0, set c from this.";
        T_ref = T_ref;
        u_ref = sqrt(kB*T_ref/m_ref);
        lightspeed = c_ref/u_ref;
    } else if (c > 0) { // default
        lightspeed = c;
        //Print() << "c>0, set T_ref from this.";
        u_ref = c_ref/lightspeed;
        T_ref = ((c_ref*c_ref)/(lightspeed*lightspeed))*(m_ref/kB);
        //Print() << "T_ref = " + std::to_string(T_ref);
    } else {
        Abort("Either ref_temp or lightspeed need to be set, only set one");
    }

    // need Larmor and Debye, can set these explicitly or via dS and beta
    int set = 0;
    if ((Larmor > 0) && (Debye > 0)) {
        set += 1;
    } else if ((beta > 0) && (skin_depth > 0)) {
        set += 2;
    }

    switch (set) {
    case 0:
        Abort("plasma parameters not set, define skin_depth & beta or Larmor & Debye, choose one pair");
    case 1:
        skin_depth = Debye*lightspeed;
        beta = 2*(Larmor/skin_depth)*(Larmor/skin_depth);
        plasma_params_set = true;
        break;
    case 2:
        Debye =  skin_depth/lightspeed;
        Larmor = sqrt(beta/2.0)*skin_depth;
        break;
    case 3:
        Abort("skin_depth+beta as well as Larmor+Debye defined, choose one pair");
    }
    Real B_ref = std::sqrt(2*mu0*n_ref*m_ref*u_ref*u_ref/beta), E_ref = c_ref*B_ref;
  
    plasma_params_set = true;

    Print() << "\n===Simulation Reference parameters===\n  Mechanical\n  x_ref = "
            << std::scientific << x_ref <<"\n  rho_ref = " << std::scientific
            << rho_ref<< "\n  m_ref = " << std::scientific << m_ref <<  +"\n  T_ref = "
            << std::scientific << T_ref << "\n  n_ref = " << std::scientific << n_ref
            << "\n  n0_ref = " << std::scientific << n0 <<  "\n  t_ref = "
            << std::scientific << x_ref/u_ref << "\n\n  Electromagnetic\n  c_ref = "
            << std::scientific << c_ref << "\n  B_ref = " << std::scientific
            << B_ref << "\n  E_ref = " << std::scientific << E_ref
            << "\n\n  NonDimensional c = " << std::scientific << c
            << "\n  NonDimensional beta = " << std::scientific << beta
            << "\n  NonDimensional skin_depth = " << std::scientific << skin_depth
            << "\n  NonDimensional Larmor = " << std::scientific << Larmor
            << "\n  NonDimensional Debye = " << std::scientific << Debye << "\n";
}

void GlobalData::resize(const int size)
{
    BL_PROFILE("GlobalData::resize");
    num_solve_state = size;

    state_names.resize(num_solve_state);
}



State& GlobalData::get_state(const int idx)
{
    BL_PROFILE("GlobalData::get_state(idx)");
    return *states[idx];
}

State& GlobalData::get_state(const std::string& name)
{
    BL_PROFILE("GlobalData::get_state(name)");
    if ( state_index.find(name) == state_index.end() ) {
        Abort("Attempting to reference a state that doesn't exist");
    }

    return *states[state_index[name]];
}

Vector<int> GlobalData::get_states_index(const Vector<std::string>& names)
{
    BL_PROFILE("GlobalData::get_states_index");
    Vector<int> index;
    for (const auto& name : names) {
        if ( state_index.find(name) == state_index.end() ) {
            Abort("Attempting to reference a state that doesn't exist");
        }
        index.push_back(state_index[name]);
    }

    return index;
}


void GlobalData::write_info(nlohmann::json &js) const
{
    BL_PROFILE("GlobalData::write_info");
    // write out globally defined data

    js["num_state"] = num_solve_state;
    js["x_ref"] = x_ref;
    js["n_ref"] = n_ref;
    js["m_ref"] = m_ref;
    js["rho_ref"] = rho_ref;
    js["T_ref"] = T_ref;
    js["u_ref"] = u_ref;
    js["n0"] = n0;
    js["lightspeed"] = lightspeed;
    js["beta"] = beta;
    js["skin_depth"] = skin_depth;
    js["Larmor"] = Larmor;
    js["Debye"] = Debye;
    js["effective_zero"] = effective_zero;

}

void GlobalData::clean_up()
{
    BL_PROFILE("GlobalData::clean_up");
#ifdef AMREX_USE_EB
    eb_data.clear();
#endif
}
