#include "MFP.H"

#include "MFP_diagnostics.H"
#include "MFP_eulerian.H"
#include "MFP_field.H"
#include "MFP_lagrangian.H"
#include "MFP_state.H"

using namespace amrex;

MFP::TimeIntegrator MFP::time_integration_scheme;
int MFP::time_integration_nsteps;

sol::state MFP::lua;
std::string MFP::lua_script;

Real MFP::force_dt = 0.0;
Real MFP::cfl = 0.0;

int MFP::Cost_Idx;
bool MFP::archive_checkpoint = true;
int MFP::verbosity;
int MFP::linear_solver_verbosity;

Vector<MFP::RefineBoxType> MFP::refine_box_type;
Vector<amrex::RealBox> MFP::refine_boxes;
bool MFP::only_refine_in_box;
bool MFP::derefine_box;

bool MFP::refine_cutcells;

bool MFP::zero_dimensional = false;

IntVect MFP::tile_size;

std::map<std::string, Array<int, AMREX_SPACEDIM + 1>> MFP::plot_variables;
Vector<std::pair<std::string, Optional3D1VFunction>> MFP::plot_functions;

constexpr int MFP::level_mask_interior;
constexpr int MFP::level_mask_covered;
constexpr int MFP::level_mask_notcovered;
constexpr int MFP::level_mask_physbnd;

Vector<std::unique_ptr<State>> MFP::states;
Vector<std::string> MFP::state_names;
std::map<std::string, int> MFP::state_index;
Vector<size_t> MFP::eulerian_states;
Vector<size_t> MFP::lagrangian_states;

Vector<std::unique_ptr<Action>> MFP::actions;
Vector<std::string> MFP::source_names;
std::map<std::string, int> MFP::source_index;

#ifdef AMREX_USE_EB
Vector<DataEB> MFP::eb_def;
Vector<Vector<EBData>> MFP::eb_data;
#endif

Real MFP::x_ref = 1.0;
Real MFP::n_ref = 1.0;
Real MFP::m_ref = 1.0;
Real MFP::rho_ref = 1.0;
Real MFP::T_ref = 1.0;
Real MFP::u_ref = 1.0;
Real MFP::n0 = 1.0;
Real MFP::prs_ref = 1.0;
Real MFP::t_ref = 1.0;
Real MFP::lightspeed = 1.0;
Real MFP::beta = 1.0;
Real MFP::skin_depth = 1.0;
Real MFP::Larmor = 1.0;
Real MFP::Debye = 1.0;

bool MFP::plasma_params_set = false;

MFP::MFP() {}

MFP::MFP(Amr& papa,
         int lev,
         const Geometry& level_geom,
         const BoxArray& bl,
         const DistributionMapping& dm,
         Real time) :
    AmrLevel(papa, lev, level_geom, bl, dm, time)
{
    BL_PROFILE("MFP::MFP");

    if (level > 0) {
        flux_reg.resize(eulerian_states.size());

        for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
            const int global_idx = eulerian_states[data_idx];
            flux_reg[data_idx].define(bl,
                                      papa.boxArray(level - 1),
                                      dm,
                                      papa.DistributionMap(level - 1),
                                      level_geom,
                                      papa.Geom(level - 1),
                                      papa.refRatio(level - 1),
                                      level,
                                      desc_lst[data_idx].nComp());
        }
    }

    build_eb();
}

MFP::~MFP() {}

void MFP::init(AmrLevel& old)
{
    BL_PROFILE("MFP::init");

    auto& oldlev = dynamic_cast<MFP&>(old);

    Real dt_new = parent->dtLevel(level);
    Real cur_time = oldlev.state[0].curTime();
    Real prev_time = oldlev.state[0].prevTime();
    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    MultiFab& C_new = get_new_data(Cost_Idx);
    FillPatch(old, C_new, 0, cur_time, Cost_Idx, 0, 1);

    for (const auto& state : states) {
        int idx = state->data_idx;
        if (idx >= 0) {
            MultiFab& S_new = get_new_data(idx);
            FillPatch(old, S_new, 0, cur_time, idx, 0, S_new.nComp());
        }
    }
}

void MFP::init()
{
    BL_PROFILE("MFP::init");

    Real dt = parent->dtLevel(level);
    Real cur_time = getLevel(level - 1).state[0].curTime();
    Real prev_time = getLevel(level - 1).state[0].prevTime();
    Real dt_old = (cur_time - prev_time) / static_cast<Real>(parent->MaxRefRatio(level - 1));
    setTimeLevel(cur_time, dt_old, dt);

    MultiFab& C_new = get_new_data(Cost_Idx);
    FillCoarsePatch(C_new, 0, cur_time, Cost_Idx, 0, 1);

    for (const auto& state : states) {
        int idx = state->data_idx;
        if (idx >= 0) {
            MultiFab& S_new = get_new_data(idx);
            FillCoarsePatch(S_new, 0, cur_time, idx, 0, S_new.nComp());
        }
    }
}

void MFP::initData()
{
    BL_PROFILE("MFP::initData");

    for (const auto& state : states) { state->init_data(this, get_cum_time()); }

    // initialise cost state
    MultiFab& C_new = get_new_data(Cost_Idx);
    C_new.setVal(0.0);
}

void MFP::post_regrid(int lbase, int new_finest)
{
    BL_PROFILE("MFP::post_regrid");

#ifdef AMREX_PARTICLES
    for (const auto& i : lagrangian_states) {
        LagrangianState& istate = LagrangianState::get_state_global(i);
        istate.redistribute(0, -1, 0);
    }
#endif
}

void MFP::post_timestep(int iteration)
{
    BL_PROFILE("MFP::post_timestep");

    // reflux
    if (level < parent->finestLevel()) {
        MFP& fine_level = getLevel(level + 1);

        for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
            EulerianState& istate = EulerianState::get_state(data_idx);
            if (istate.reflux) {
                MultiFab& S_crse = get_new_data(data_idx);
#ifdef AMREX_USE_EB
                MultiFab& S_fine = fine_level.get_new_data(data_idx);
                EBData& eb = get_eb_data(istate.global_idx);
                fine_level.flux_reg[data_idx].Reflux(S_crse, eb.volfrac, S_fine, eb.volfrac);
#else
                fine_level.flux_reg[data_idx].Reflux(S_crse, 0);
#endif
            }
        }
    }

    if (level < parent->finestLevel()) { avgDown(); }

#ifdef AMREX_PARTICLES
    const int ncycle = parent->nCycle(level);
    // Don't redistribute on the final subiteration except on the coarsest grid.
    if (iteration < ncycle || level == 0) {
        int ngrow = (level == 0) ? 0 : iteration;

        for (const auto& i : lagrangian_states) {
            LagrangianState& istate = LagrangianState::get_state_global(i);
            istate.redistribute(level, parent->finestLevel(), ngrow);
        }
    }
#endif
}

void MFP::avgDown()
{
    BL_PROFILE("MFP::avgDown");

    if (level == parent->finestLevel()) return;

    auto& fine_lev = getLevel(level + 1);

    for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
        const int global_idx = eulerian_states[data_idx];

        MultiFab& S_crse = get_new_data(data_idx);
        MultiFab& S_fine = fine_lev.get_new_data(data_idx);

#ifdef AMREX_USE_EB
        MultiFab volume(S_fine.boxArray(), S_fine.DistributionMap(), 1, 0);
        volume.setVal(1.0);

        EBData& eb = get_eb_data(level + 1, global_idx);
        amrex::EB_average_down(S_fine, S_crse, volume, eb.volfrac, 0, S_fine.nComp(), fine_ratio);
#else
        amrex::average_down(S_fine, S_crse, 0, S_fine.nComp(), fine_ratio);
#endif
    }
}

void MFP::postCoarseTimeStep(Real time)
{
    BL_PROFILE("MFP::postCoarseTimeStep");

    if (verbosity >= 1) {
        amrex::Print().SetPrecision(17)
          << "[MFP] : step = " << parent->levelSteps(0) << ", time = " << time << std::endl;
    }
}

void MFP::post_init(Real)
{
    BL_PROFILE("MFP::post_init");

    // check our reconstruction doesn't go out-of-bounds

    const Box& box = geom.Domain();
    const IntVect size = box.bigEnd() - box.smallEnd();
    for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
        EulerianState& istate = EulerianState::get_state(data_idx);

        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (istate.reconstructor->get_num_grow() > size[d])
                Abort("Reconstruction stencil is larger than domain in dim=" + num2str(d));
        }
    }

    if (level > 0) return;
    for (int k = parent->finestLevel() - 1; k >= 0; --k) { getLevel(k).avgDown(); }
}

void MFP::post_restart()
{
    BL_PROFILE("MFP::post_restart");

    if (level > 0) {
        flux_reg.resize(eulerian_states.size());
        for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
            EulerianState& istate = EulerianState::get_state(data_idx);
            if (istate.reflux) {
                flux_reg[data_idx].define(grids,
                                          parent->boxArray(level - 1),
                                          dmap,
                                          parent->DistributionMap(level - 1),
                                          geom,
                                          parent->Geom(level - 1),
                                          parent->refRatio(level - 1),
                                          level,
                                          desc_lst[data_idx].nComp());
            }
        }
    }

    build_eb();

#ifdef AMREX_PARTICLES
    std::string restart_chkfile;
    ParmParse pp("amr");
    pp.query("restart", restart_chkfile);
    ParticlePostRestart(restart_chkfile);
#endif
}

#ifdef AMREX_PARTICLES

void MFP::writeParticles(const std::string& dir)
{
    BL_PROFILE("MFP::writeParticles");

    for (const auto& i : lagrangian_states) {
        LagrangianState& istate = LagrangianState::get_state_global(i);
        istate.checkpoint(dir);
    }
}

void MFP::ParticlePostRestart(const std::string& dir)
{
    BL_PROFILE("MFP::ParticlePostRestart");

    if ((level == 0) && !lagrangian_states.empty()) {
        // handle if we have archived level data in a restart folder
        Vector<std::string> to_remove;
        if (ParallelDescriptor::IOProcessor()) {
            // get the names of the particle folders
            for (const auto& i : lagrangian_states) {
                LagrangianState& istate = LagrangianState::get_state_global(i);
                std::string particle_folder = "Particles_" + istate.name;
                std::string FullPath = dir + "/" + particle_folder;

                // check if we actually need to un-tar
                if (FileExists(FullPath)) continue;

                // check tar exists
                if (!FileExists(FullPath + ".tar")) continue;

                // add to list of things to clean up later
                to_remove.push_back(FullPath);

                // perform un-tar
                std::string cmd = "\\tar -C " + dir + " -xf " + FullPath + ".tar ";
                const char* command = {cmd.c_str()};
                int retVal = std::system(command);
                if (retVal == -1 || WEXITSTATUS(retVal) != 0) {
                    Abort("Error: Unable to un-tar '" + FullPath + "'");
                }
            }
        }

        ParallelDescriptor::Barrier();

        // retrieve particle data
        for (const auto& i : lagrangian_states) {
            LagrangianState& istate = LagrangianState::get_state_global(i);
            istate.init(this, false);
            istate.restart(dir);
        }

        // clean up after un-tar operation
        if (ParallelDescriptor::IOProcessor()) {
            for (const auto& path : to_remove) { FileSystem::RemoveAll(path); }
        }
    }
}

#endif

State& MFP::get_state(const std::string& name)
{
    if (state_index.find(name) == state_index.end()) {
        Abort("Attempting to reference a state that doesn't exist");
    }
    return *states[state_index[name]];
}

Action& MFP::get_source(const std::string& name)
{
    if (source_index.find(name) == source_index.end()) {
        Abort("Attempting to reference a source that doesn't exist");
    }
    return *actions[source_index[name]];
}
