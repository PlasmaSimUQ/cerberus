#include <AMReX_ParmParse.H>

#include "MFP.H"
#include "MFP_state.H"
#include "MFP_eb_sdf.H"
#include "MFP_read_geom.h"
#include "MFP_utility.H"

#include "MFP_state.H"
#include "MFP_action.H"


void MFP::set_lua_script(const std::string &script)
{
    BL_PROFILE("MFP::set_lua_script");
    lua_script =
        #include "default.lua"
            ;

    lua_script += script;

    // insert call to preprocessor defined in default.lua
    lua_script += "\npreprocess()\n";
}

void MFP::save_lua_script() {

    if (ParallelDescriptor::IOProcessor()) {

        std::ofstream ofs;

        ofs.open("MFP_Lua_initialisation.lua", std::ofstream::out | std::ofstream::trunc);
        ofs << lua_script;
        ofs.close();
        Print() << "'MFP_Lua_initialisation.lua' written to file\n";

    }

    return;
}

void MFP::read_config()
{
    BL_PROFILE("MFP::read_config");


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

    verbosity = lua["verbosity"];

    std::string time_integrator = lua["time_integration_scheme"];

    if (time_integrator == "one_step") {
        time_integration_scheme = TimeIntegrator::OneStep;
        need_scratch_space = true;
    } else if (time_integrator == "strang") {
        time_integration_scheme = TimeIntegrator::StrangSplitting;
        need_scratch_space = true;
    } else if (time_integrator == "symplectic") {
        time_integration_scheme = TimeIntegrator::Symplectic;
        need_scratch_space = false;
    } else {
        Abort("Time integration scheme '"+time_integrator+"' is not recognised, try ['one_step', 'strang', 'symplectic']");
    }

    linear_solver_verbosity = lua["linear_solver_verbosity"];

    //
    // reference variables
    //

    update_ref();

    //
    // generate states
    //

    // what states do we have?
    sol::table state_def_names = lua.script("return get_sorted_keys(states)");

    ClassFactory<State> sfact = GetStateFactory();

    // get a list of all the tags
    Vector<std::string> state_tags;
    for (const auto& S : sfact.getRegistered()) {
        state_tags.push_back(S.first);
    }


    int global_idx = 0;
    for (auto& item : state_def_names) {
        std::string name = item.second.as<std::string>();

        sol::table state_def = lua["states"][name];

        state_def["name"] = name;
        state_def["global_idx"] = global_idx;

        std::string state_type = state_def["type"].get<std::string>();

        std::unique_ptr<State> istate = sfact.Build(state_type, state_def);

        if (!istate)
            Abort("Failed to read state "+name+", must be one of "+vec2str(state_tags));

        states.push_back(std::move(istate));
        state_names.push_back(name);
        state_index[name] = global_idx;

        State& S = get_state(name);

        switch (S.get_classification()) {
        case State::StateClassification::Eulerian:
            eulerian_states.push_back(global_idx);
            break;
        case State::StateClassification::Lagrangian:
            lagrangian_states.push_back(global_idx);
            break;
        default:
            Abort("How did we get here?");
        }



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

    //
    // miscellaneous
    //

    zero_dimensional = lua["zero_dimensional"];

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        tile_size[i] = lua["tile_size"][i + 1];
    }

    force_dt = lua["force_dt"];

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
                eb_dat.states.push_back({istate.global_idx,istate.get_eb_bc_size()});

                // define bc
                istate.set_eb_bc(bc_def);
            }
            ++j;
        }

        eb_def.push_back(eb_dat);

        ++i;
    }
#endif


    //
    // generate actions
    //

    // what actions do we have?
    sol::table act_def_names = lua.script("return get_sorted_keys(actions)");

    ClassFactory<Action> act_fact = GetActionFactory();

    // get a list of all the tags
    Vector<std::string> act_tags;
    for (const auto& S : act_fact.getRegistered()) {
        act_tags.push_back(S.first);
    }


    int act_idx = 0;
    for (auto& item : act_def_names) {
        std::string act_name = item.second.as<std::string>();

        sol::table act_def = lua["actions"][act_name];

        act_def["name"] = act_name;
        act_def["src_idx"] = act_idx;

        std::string src_type = act_def["type"].get<std::string>();

        std::unique_ptr<Action> act = act_fact.Build(src_type, act_def);

        if (!act)
            Abort("Failed to read source "+act_name+", must be one of "+vec2str(act_tags));

        actions.push_back(std::move(act));
        source_names.push_back(act_name);
        source_index[act_name] = act_idx;

        act_idx++;
    }



    // update anything in the states that requires actions to be defined
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

        for (const auto& state : states) {
            int idx = state->data_idx;
            for (auto& name : state->get_plot_output_names()) {
                plot_variables[name+"-"+state->name][0] = 0;
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

void MFP::update_ref()
{
    BL_PROFILE("MFP::update_ref");
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

// pointer to state initialization routine
typedef void (*state_config)(const std::string &name, const int idx);

void MFP::read_params()
{
    BL_PROFILE("MFP::read_params()");

    ParmParse mfp("mfp");
    std::string lua_script;

    if (mfp.contains("lua")) {
        mfp.get("lua", lua_script);
    } else {
        amrex::Abort("Exiting as we have not defined the problem (mfp.lua = '')");
    }

    set_lua_script(lua_script);

#ifdef AMREX_DEBUG
    save_lua_script();
#endif

    ParmParse amr("amr");

    bool plot_output = false;
    amr.query("plot_files_output", plot_output);

    amr.query("check_archive", archive_checkpoint);

    // pass off to global data routine
    read_config();

}
