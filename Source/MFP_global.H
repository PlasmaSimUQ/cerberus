#ifndef MFP_STATE_H
#define MFP_STATE_H

#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <string>

#include "MFP_utility.H"

#include "sol.hpp"

#include "MFP_transforms.H"

//TODO check if could be problematic
#include "MFP_state.H"
#include "MFP_ode_system.h"

#include "MFP_hydro.H"
#include "MFP_hydro2p.H"
#include "MFP_mhd.H"
#include "MFP_field.H"

#include "MFP_eb_sdf.H"

using namespace amrex;

enum LevelMask {
    Inside = 0, // valid cells
    Covered,    // ghost cells covered by valid cells of this level
    NotCovered, // ghost cells not covered
    PhysBnd,    // outside domain
    NUM
};

#ifdef AMREX_USE_EB

struct EBData {
    std::unique_ptr<EBFArrayBoxFactory> ebfactory;
    FabArray<EBCellFlagFab> flags;
    MultiFab volfrac;
    const amrex::MultiCutFab* bndrycent;
    const amrex::MultiCutFab* bndrynorm;
    Array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
    Array<const MultiCutFab*, AMREX_SPACEDIM> facecent;
    MultiCutFab bndryidx;
};

struct DataEB {
    int index;
    const EB2::IndexSpace* index_space;
    std::function<Real (AMREX_D_DECL(Real, Real, Real))> geom_func;
    int insertion_type; // add = 1, subtract = -1
    bool inside;
    Vector<std::pair<int,int>> states;
    bool regular;
};
#endif

enum class RefineBoxType : int {
    ForceRefine=0,
    NoRefine,
    OnlyRefine,
    NUM
};


//-----------------------
// holds everything that is possible without needing any local dependencies
// - acts as the go-between

class GlobalData
{
public:
    GlobalData();
    virtual ~GlobalData();

public:

    static sol::state lua;
    static std::string lua_script;  // lua script that defines our problem

    static int Cost_Idx;
    static int Shock_Idx;
    static int num_states;
    static int num_solve_state;

    static int full_state_size;

    static int zero_dimensional;

    static Array<bool,AMREX_SPACEDIM> periodic;

    static Vector<std::string> state_names;
    static std::map<std::string, int> state_index;
    static Vector<std::string> state_tags;
    static bool detect_shocks;
    static bool plot_shock_detector;
    static int num_shock_detector;

    static int num_fluid; // how many fluid type states are there

    static int do_CTU;

    static int do_face_src;

    static Real force_dt;

    // Parameters
    static int verbose;
    static int linear_solver_verbosity;
    static IntVect tile_size;
    static Real cfl;

    static int refine_cutcells;

    static Vector<RefineBoxType> refine_box_type;
    static Vector<amrex::RealBox> refine_boxes;
    static bool only_refine_in_box;
    static bool derefine_box;

#ifdef AMREX_PARTICLES
    static int do_tracer_particles;
    static int particle_verbose;
#endif

    static Real x_ref, n_ref, m_ref, rho_ref, T_ref, u_ref, n0;
    static Real lightspeed, beta, skin_depth, Larmor, Debye;

    static bool plasma_params_set;

    static Real effective_zero;

    static int msg;


    static Vector<std::unique_ptr<State>> states;
    static Vector<std::unique_ptr<ODESystem>> ode_source_terms;

#ifdef AMREX_USE_EB
    static Vector<DataEB> eb_def;
    static Vector<Vector<EBData>> eb_data;
#endif

    static std::map<std::string,Array<int,AMREX_SPACEDIM+1>> plot_variables;
    static Vector<std::pair<std::string,Optional3D1VFunction>> plot_functions;

    // functions

    static void set_num_levels(int n);

    static void set_lua_script(const std::string &script);
    static void read_config(const Vector<int> &is_periodic, const bool plot_output=false);
    static void init_source_terms();


    static void update_ref();

    static void resize(const int size);

    static State& get_state(const int idx);
    static State& get_state(const std::string &name);
    Vector<int> get_states_index(const Vector<std::string>& names);

    void write_info(nlohmann::json& js) const;

    static void clean_up();

};

#endif // MFP_STATE_H
