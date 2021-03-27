
#include <MFP.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

using namespace amrex;


// Iterate over all boxes within the level to perform a single time-step update

Real MFP::advance(Real time, Real dt, int iteration, int ncycle) {
    BL_PROFILE("MFP::advance()");

    if (Level() == 0) {
        solve_static_fields(time);
    }

    ParallelDescriptor::Barrier();


    MultiFab& C_new = get_new_data(gd.Cost_Idx);
    C_new.setVal(0.0);


    // ==========================================================================
    // A. source update

    apply_cell_sources(time, dt/2.0);

    // ==========================================================================
    // B. flux update

    if (!gd.zero_dimensional) {
        apply_cell_transport(time, dt);
    }

    // ==========================================================================
    // C. update from cell source terms

    apply_cell_sources(time, dt/2.0);


    // ==========================================================================

    for (int i = 0; i < gd.num_states; ++i) {
        state[i].setNewTimeLevel(time+dt);
    }

    // do some lua garbage collection
    gd.lua.script("collectgarbage('collect')");
    gd.lua.script("collectgarbage('collect')");
    return dt;
}

