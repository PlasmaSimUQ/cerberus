
#include <AMReX.H>

#ifdef AMREX_USE_EB
    #include <AMReX_EB2.H>
#endif
#include "MFP.H"
#include "MFP_eulerian.H"

#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#ifdef PPROF
    #include <gperftools/profiler.h>
#endif

using namespace amrex;

amrex::LevelBld* getLevelBld();

#ifdef AMREX_USE_EB
void initialize_EB2(Amr& amr, int ngrow);
#endif

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    BL_PROFILE_VAR("main()", pmain);

#ifdef PPROF
    ProfilerStart("gperftools_profile.prof");
#endif

    Real timer_tot = ParallelDescriptor::second();
    Real timer_init = 0.;
    Real timer_advance = 0.;

    int max_step;
    Real strt_time;
    Real stop_time;

    {
        ParmParse pp;

        max_step = -1;
        strt_time = 0.0;
        stop_time = -1.0;

        pp.query("max_step", max_step);
        pp.query("strt_time", strt_time);
        pp.query("stop_time", stop_time);

        Vector<float> dx;
        //     if (nx=pp.countval("cell_size")) {
        //        // get nx values starting at index 0 and store in dx.
        //        // dx is automatically resized here.
        //        pp.getarr("cell_size",dx,0,nx);
        //     }
    }

    if (strt_time < 0.0) { amrex::Abort("MUST SPECIFY a non-negative strt_time"); }

    if (max_step < 0 && stop_time < 0.0) {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

#ifdef AMREX_USE_EB
    {
        // check that all refinement ratios are <= 2 for EB
        Vector<int> ratios;
        ParmParse pp("amr");
        pp.queryarr("ref_ratio", ratios);

        for (const int& r : ratios)
            if (r > 2) amrex::Abort("Exiting because EB is enabled and refinement ratio > 2");
    }
#endif

    {
        timer_init = ParallelDescriptor::second();

        Amr amr(getLevelBld());

        int n_grow = 4;

        for (int data_idx = 0; data_idx < MFP::eulerian_states.size(); ++data_idx) {
            EulerianState& istate = EulerianState::get_state(data_idx);
            n_grow = std::max(n_grow, istate.get_num_grow() + 3);
        }

#ifdef AMREX_USE_EB
        AmrLevel::SetEBSupportLevel(EBSupport::full);
        AmrLevel::SetEBMaxGrowCells(n_grow, n_grow, n_grow);

        initialize_EB2(amr, n_grow);
#endif

        // handle if we have archived level data in a restart folder
        Vector<std::string> to_remove;
        if (ParallelDescriptor::IOProcessor()) {
            const std::string restart_chkfile = amr.theRestartFile();
            if (!restart_chkfile.empty() && restart_chkfile != "init") {
                // un-tar all of the levels
                std::string LevelDir, FullPath;
                for (int ilev = 0; ilev <= amr.maxLevel(); ++ilev) {
                    // get the name of the level folder
                    LevelDir = amrex::Concatenate("Level_", ilev, 1);
                    FullPath = restart_chkfile;
                    if (!FullPath.empty() && FullPath.back() != '/') { FullPath += '/'; }
                    FullPath += LevelDir;

                    // check if we actually need to un-tar
                    if (FileExists(FullPath)) continue;

                    // check tar exists
                    if (!FileExists(FullPath + ".tar")) continue;

                    // add to list of things to clean up later
                    to_remove.push_back(FullPath);

                    // perform un-tar
                    std::string cmd = "\\tar -C " + restart_chkfile + " -xf " + FullPath + ".tar ";
                    const char* command = {cmd.c_str()};
                    int retVal = std::system(command);
                    if (retVal == -1 || WEXITSTATUS(retVal) != 0) {
                        Abort("Error: Unable to un-tar '" + FullPath + "'");
                    }
                }
            }
        }

        // make sure any tar operations have completed
        ParallelDescriptor::Barrier();

        // perform initialization
        amr.init(strt_time, stop_time);

        // clean up after un-tar operation
        if (ParallelDescriptor::IOProcessor()) {
            for (const auto& path : to_remove) { FileSystem::RemoveAll(path); }
        }

        timer_init = ParallelDescriptor::second() - timer_init;

        timer_advance = ParallelDescriptor::second();

        while (amr.okToContinue() && (amr.levelSteps(0) < max_step || max_step < 0) &&
               (amr.cumTime() < stop_time || stop_time < 0.0))

        {
            //
            // Do a coarse timestep.  Recursively calls timeStep()
            //
            amr.coarseTimeStep(stop_time);
        }

        timer_advance = ParallelDescriptor::second() - timer_advance;

        // Write final checkpoint and plotfile
        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) { amr.checkPoint(); }

        if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) { amr.writePlotFile(); }
    }

    timer_tot = ParallelDescriptor::second() - timer_tot;

    ParallelDescriptor::ReduceRealMax({timer_tot, timer_init, timer_advance},
                                      ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run Time total        = " << timer_tot << "\n"
                   << "Run Time init         = " << timer_init << "\n"
                   << "Run Time advance      = " << timer_advance << "\n";

    BL_PROFILE_VAR_STOP(pmain);

#ifdef PPROF
    ProfilerStop();
#endif

    amrex::Finalize();

    return 0;
}
