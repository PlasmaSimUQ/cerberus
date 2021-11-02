#ifndef MFP_H_
#define MFP_H_

#ifdef AMREX_USE_EB
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFluxRegister.H>
#else
#include <AMReX_YAFluxRegister.H>
#endif

#include <AMReX_AmrLevel.H>
#include <AMReX_LO_BCTYPES.H>
#include <memory>
#include "sol.hpp"
#include "json.hpp"

#include "MFP_eb_structs.H"
#include "MFP_optional_func.H"
#include "MFP_action.H"
#include "MFP_reconstruction.H"

class State;
class Action;

using namespace amrex;

class MFP : public amrex::AmrLevel {
public:
    MFP();
    MFP(amrex::Amr& papa, int lev, const amrex::Geometry& level_geom,
        const amrex::BoxArray& bl, const amrex::DistributionMapping& dm,
        amrex::Real time);
    virtual ~MFP();

    MFP(const MFP& rhs) = delete;
    MFP& operator=(const MFP& rhs) = delete;

    std::string getVersion();

    virtual std::string thePlotFileType() const override {
        return {"HyperCLaw-V1.1"};
    }

#ifdef AMREX_USE_EB
    bool check_covered_stencil(Array4<const EBCellFlag> const& flag,
                               int i, int j, int k, int d,
                               int stencil_length);
#endif

    void calc_slope(const Box& box,
                    const FArrayBox& src,
                    FArrayBox &slope,
                #ifdef AMREX_USE_EB
                const EBCellFlagFab &flag,
                #endif
                    const Real *dx,
                    int index,
                    int dim,
                    Reconstruction &reco);

    void getPlotData(MultiFab &plot_data, std::vector<std::string> &plot_names);

    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how = VisMF::NFiles) override;

    virtual void writePlotFilePost(const std::string &dir, std::ostream &os) override;

    virtual void checkPointPost (const std::string& dir,std::ostream& os) override;

    void archive_folder(const std::string &dir);

    // Initialize data on this level from another MFP (during regrid).
    virtual void init(amrex::AmrLevel& old) override;

    // Initialize data on this level after regridding if old level did not
    // previously exist
    virtual void init() override;

    // Initialize grid data at problem start-up.
    virtual void initData() override;

    // Advance grids at this level in time.
    virtual amrex::Real advance(amrex::Real time, amrex::Real dt, int iteration,
                                int ncycle) override;

    enum class TimeIntegrator {
        OneStep = 0,
        StrangSplitting,
        Symplectic,
    };

    static TimeIntegrator time_integration_scheme;

    void advance_euler(Real time, Real dt, int iteration, int ncycle);
    void advance_strang(Real time, Real dt, int iteration, int ncycle);
    void apply_derivative(Vector<std::pair<int, MultiFab> > &dU);

    Real estTimeStep();

#ifdef AMREX_PARTICLES
    void init_particles();
    void writeParticles(const std::string& dir);
    void ParticlePostRestart (const std::string &dir);
#endif

    virtual void computeInitialDt(int finest_level, int sub_cycle,
                                  amrex::Vector<int>& n_cycle,
                                  const amrex::Vector<amrex::IntVect>& ref_ratio,
                                  amrex::Vector<amrex::Real>& dt_level,
                                  amrex::Real stop_time) override;

    virtual void computeNewDt(int finest_level, int sub_cycle,
                              amrex::Vector<int>& n_cycle,
                              const amrex::Vector<amrex::IntVect>& ref_ratio,
                              amrex::Vector<amrex::Real>& dt_min,
                              amrex::Vector<amrex::Real>& dt_level,
                              amrex::Real stop_time,
                              int post_regrid_flag) override;

    virtual void post_regrid(int lbase, int new_finest) override;

    // Do work after timestep().
    virtual void post_timestep(int iteration) override;

    // After a full time step
    virtual void postCoarseTimeStep(amrex::Real time) override;

    // Do work after init().
    virtual void post_init(amrex::Real stop_time) override;

    virtual void post_restart() override;

    // Error estimation for regridding.
    virtual void errorEst(amrex::TagBoxArray& tb, int clearval, int tagval,
                          amrex::Real time, int n_error_buf = 0,
                          int ngrow = 0) override;

    virtual int WorkEstType() override { return 0; }

    static void read_params();

    int get_level() {return level;}
    Amr* get_parent() {return parent;}
    Real get_cum_time() {return parent->cumTime();}
    Real get_dt() {return parent->dtLevel(level);}

    // Define data descriptors.
    static DescriptorList& get_desc_lst () noexcept { return desc_lst; }
    static void variableSetUp();

    // Cleanup data descriptors at end of run.
    static void variableCleanUp();

    MFP& getLevel(int lev) { return dynamic_cast<MFP&>(parent->getLevel(lev)); }

    static State* get_state_ptr(const int idx) {return states[idx].get();}
    static State& get_state(const int idx) { return *states[idx]; }
    static State& get_state(const std::string& name);

    static Action* get_source_ptr(const int idx) {return actions[idx].get();}
    static Action& get_source(const int idx) { return *actions[idx]; }
    static Action& get_source(const std::string& name);

    void build_eb();
    void avgDown();

    static void set_lua_script(const std::string &script);
    static void save_lua_script();
    static void read_config();
    static void update_ref();

    void write_info(nlohmann::json &js) const;

    static sol::state lua;
    static std::string lua_script;  // lua script that defines our problem

    static Real cfl;
    static Real force_dt;

    static int Cost_Idx;

    static bool archive_checkpoint;
    static int verbosity;
    static int linear_solver_verbosity;

    enum class RefineBoxType : int {
        ForceRefine=0,
        NoRefine,
        OnlyRefine,
    };

    static Vector<RefineBoxType> refine_box_type;
    static Vector<amrex::RealBox> refine_boxes;
    static bool only_refine_in_box;
    static bool derefine_box;

    static bool refine_cutcells;

    static bool zero_dimensional;

    static IntVect tile_size;

    static Vector<std::unique_ptr<State>> states;
    static Vector<std::string> state_names;
    static std::map<std::string, int> state_index;
    static Vector<size_t> eulerian_states; // mapping from data_idx to global_idx
    static Vector<size_t> lagrangian_states;

    static bool need_scratch_space;
    Vector<std::pair<int,MultiFab>> eulerian_du;

    static Vector<std::unique_ptr<Action>> actions;
    static Vector<std::string> source_names;
    static std::map<std::string, int> source_index;

    static std::map<std::string,Array<int,AMREX_SPACEDIM+1>> plot_variables;
    static Vector<std::pair<std::string,Optional3D1VFunction>> plot_functions;

    static Real x_ref, n_ref, m_ref, rho_ref, T_ref, u_ref, n0;
    static Real lightspeed, beta, skin_depth, Larmor, Debye;
    static bool plasma_params_set;

    enum LevelMask {
        Inside = 0, // valid cells
        Covered,    // ghost cells covered by valid cells of this level
        NotCovered, // ghost cells not covered
        PhysBnd,    // outside domain
    };

    static constexpr int level_mask_interior =   0;  // valid cells
    static constexpr int level_mask_covered =    1;  // ghost cells covered by valid cells of this level
    static constexpr int level_mask_notcovered = 2;  // ghost cells not covered
    static constexpr int level_mask_physbnd =    3;  // outside domain

    iMultiFab level_mask;

#ifdef AMREX_USE_EB
    static Vector<DataEB> eb_def;
    static Vector<Vector<EBData>> eb_data;
    EBData& get_eb_data(const int idx) { return eb_data[level][idx]; }
    EBData& get_eb_data(const int lev, const int idx) { return eb_data[lev][idx]; }
#endif


#ifdef AMREX_USE_EB
    std::vector<amrex::EBFluxRegister> flux_reg;
#else
    std::vector<amrex::YAFluxRegister> flux_reg;
#endif

};

#endif
