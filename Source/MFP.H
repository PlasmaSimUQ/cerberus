#ifndef MFP_H_
#define MFP_H_

#include "MFP_global.H"

#ifdef AMREX_USE_EB
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFluxRegister.H>
#else
#include <AMReX_YAFluxRegister.H>
#endif

#include <AMReX_AmrLevel.H>
#include <AMReX_LO_BCTYPES.H>
#include <memory>

#include "MFP_fillbc.H"
#include "MFP_reconstruction.H"
#include "MFP_Riemann.H"
#include "MFP_source.H"
#include "MFP_ode_solver.h"
#include "MFP_viscous.H"
#include "MFP_diagnostics.H"

#ifdef AMREX_PARTICLES
#include "AMReX_AmrParticles.H"
using ParTileType = ParticleTile<AMREX_SPACEDIM, 0, 0, 0>;
using ParIterType = ParIter<AMREX_SPACEDIM, 0, 0, 0>;
using ParticleType = Particle<AMREX_SPACEDIM, 0>;
#endif

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

    // function for compiling a multifab that holds all of the data that needs
    // to be plotted for this level
    void getPlotData(MultiFab& plot_data, std::vector<std::string>& plot_names);

    virtual void writePlotFile (const std::string& dir,
                                std::ostream&      os,
                                VisMF::How         how = VisMF::NFiles) override;

    virtual void writePlotFilePost(const std::string &dir, std::ostream &os) override;

    virtual void checkPointPost (const std::string& dir,std::ostream&      os) override;

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

    virtual int WorkEstType() override { return gd.Cost_Idx; }

    // Define data descriptors.
    static void variableSetUp();

    // Cleanup data descriptors at end of run.
    static void variableCleanUp();

protected:

    static bool first_step;

    static GlobalData gd;

    static void get_lua_script();
    static void save_lua_script();

    // initialization routines for all of the different types of states

    //    static void init_source_terms();

    static void read_params();

    MFP& getLevel(int lev) { return dynamic_cast<MFP&>(parent->getLevel(lev)); }

    void avgDown();


    amrex::iMultiFab level_mask;
    void buildMetrics();


    amrex::Real estTimeStep();

    // Compute initial time step.
    amrex::Real initialTimeStep();

    void printTotal();

    bool is_vacuum(const Box& bx, const FArrayBox& fab, const int idx);

#ifdef AMREX_PARTICLES
    static int particle_verbosity;
    static Vector<std::string> particle_names;
    static Vector<int> particle_idx;
    static Vector<amrex::AmrTracerParticleContainer*> particles;

    static void read_particle_params ();
    void init_particles();
    void writeParticles(const std::string& dir);
    void ParticlePostRestart (const std::string &dir);
    void advance_particles (int iteration, amrex::Real time, amrex::Real dt);
    void AdvectWithUcc (AmrTracerParticleContainer* pc,
                        const MultiFab& Ucc,
                        int lev,
                        Real dt
                        EB_OPTIONAL(,const FabArray<EBCellFlagFab>& flags)
                        );

    void push_particles (ParTileType& ptile,
                         const FArrayBox& prim,
                         Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                         Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                         const int vel_idx,
                         const Real dt
                         EB_OPTIONAL(,const EBCellFlagFab& flag)
                         );

#endif

#ifdef AMREX_USE_EB
    std::vector<amrex::EBFluxRegister> flux_reg;

    Vector<EBData>& getEBData() {return gd.eb_data[level];}
    EBData& getEBData(const int idx) {return gd.eb_data[level][idx];}
    EBData& getEBData(const int ilev, const int idx) {return gd.eb_data[ilev][idx];}
    EBData& getFineEBData(const int idx) {return gd.eb_data[level+1][idx];}
    FArrayBox get_eb_faces(const MFIter &mfi, const int idx);
    void tag_cut_cells (TagBoxArray& tags, const int idx);
#else
    std::vector<amrex::YAFluxRegister> flux_reg;
#endif

    static constexpr int level_mask_interior =   0;  // valid cells
    static constexpr int level_mask_covered =    1;  // ghost cells covered by valid cells of this level
    static constexpr int level_mask_notcovered = 2;  // ghost cells not covered
    static constexpr int level_mask_physbnd =    3;  // outside domain

    static bool archive_checkpoint;

    // calculate the starting condition for all states via the provided lua
    // functions for each component of the state vector
    void init_data(const Box& box,
                   EB_OPTIONAL(const EBCellFlagFab &flag,)
                   FArrayBox& src,
                   const Real* dx,
                   const Real* prob_lo,
                   const int idx);

    // provide an estimate of the minimum time step required for the stable
    // computation of a single time update
    void estimate_dt (const Box& box,
                      Vector<const FArrayBox*>& src,
                      EB_OPTIONAL(Vector<const FArrayBox*>& vfrac,)
                      const Real* dx,
                      Real& dt);

    // elliptic solutions of field variables

    void solve_static_fields(const Real time);

    void correct_dynamic_fields(const Real time);


    void project_divergence(const Real time);

    void solve_divergence(const int state_idx,
                          const int vector_idx,
                          const int phi_idx,
                          const BCRec &bc,
                          Vector<MultiFab*> S_cc_ptr={});

    // calculate any modifications to face values due to source terms
    void calc_face_source (const Box& box,
                           Vector<FArrayBox*>& src_dat,
                           Vector<Array<FArrayBox, AMREX_SPACEDIM> > &R_lo,
                           Vector<Array<FArrayBox, AMREX_SPACEDIM> > &R_hi,
                           EB_OPTIONAL(const Vector<const EBCellFlagFab*> &flag,)
                           const Real* dx,
                           Real time,
                           Real dt);

    // calculate source terms for cell centred values
    void calc_cell_source (const Box& box,
                           Vector<FArrayBox*>& src_dat,
                           Vector<FArrayBox*>& dst_dat,
                           EB_OPTIONAL(Vector<const EBCellFlagFab*> &flags,)
                           EB_OPTIONAL(Vector<const FArrayBox*> &volfrac,)
                           Real time, Real dt);

    // calculate and apply cell source terms, this acts directly on data at 'time'
    void apply_cell_sources(const Real time,
                            const Real dt);

    // calculate and apply cell transport terms, this acts directly on data at 'time'
    void apply_cell_transport(amrex::Real time,
                              amrex::Real dt);


    // tag cells for refinement according to the passed in configuration parameters
    void tag_refinement(const Box& box,
                        const FArrayBox &src,
                        EB_OPTIONAL(const EBCellFlagFab& flags,)
                        TagBox &tags,
                        const int n,
                        const Real threshold,
                        const Real min_val,
                        const char tagval,
                        const char clearval);

};

#endif
