#include <MFP.H>
#include <MFP_utility.H>
#include <AMReX_ParmParse.H>

#include <tuple>

#include "MFP_lua.H"
#include "MFP_solvers.H"
#include "MFP_source.H"
#include "MFP_sources.H"
#include "MFP_interp.H"

using namespace amrex;

// pointer to state initialization routine
typedef void (*state_config)(const std::string &name, const int idx);

void MFP::read_params() {
    BL_PROFILE("MFP::read_params()");

    ParmParse mfp("mfp");
    std::string lua_script;

    if (mfp.contains("lua")) {
        mfp.get("lua", lua_script);
    } else {
        amrex::Abort("Exiting as we have not defined the problem (mfp.lua = '')");
    }

    gd.set_lua_script(lua_script);

#ifdef AMREX_DEBUG
    save_lua_script();
#endif

    ParmParse pgm("geometry");

    Vector<int> is_per(AMREX_SPACEDIM,0);
    pgm.queryarr("is_periodic",is_per);

    ParmParse amr("amr");

    bool plot_output = false;
    amr.query("plot_files_output", plot_output);

    amr.query("check_archive", archive_checkpoint);

    // pass off to global data routine
    gd.read_config(is_per, plot_output);

    // === OTHER ===

    // get particle params
#ifdef AMREX_PARTICLES
    read_particle_params();
#endif

}
