#ifndef MFP_GENERIC_H
#define MFP_GENERIC_H

#include <set>

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

#ifdef AMREX_USE_EB
#include <AMReX_MultiCutFab.H>
#include "MFP_eb_bc.H"
#include "MFP_eb_divergence.H"
#endif


#include "forward.hpp"
#include "json.hpp"

#include "MFP_Riemann.H"
#include "MFP_bc.H"
#include "MFP_reconstruction.H"

#include "MFP_Riemann.H"
#include "MFP_viscous.H"
#include "MFP_shockdetector.H"

#include "MFP_fillbc.H"
#include "MFP_optional_func.H"
#include "MFP_utility.H"

#ifdef AMREX_USE_EB
#define EB_OPTIONAL(...) __VA_ARGS__
#else
#define EB_OPTIONAL(...)
#endif

#ifdef PYTHON
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#endif

using namespace amrex;
using dual = autodiff::dual;

enum class StateType {
    isNull = -1,
    isHydro,
    isHydro2P,
    isMHD,
    isField,
    NUM
};

enum class AssociatedType {
    isNull = -1,
    Neutral,
    Ion,
    Electron,
    Field,
    MHD,
};

struct BoundaryInfo
{
        int imin, imax;
        int jmin, jmax;
        int kmin, kmax;
        int dir, lo_hi;
};

/*
 * The State class is intended to be used as a base for the definition of
 * all the various states that will be used to define fluids within the
 * simulation. All of the members are public to allow for simple access
 * to the data within a state.
 */

class State
{
    public:
        State();
        ~State();

        std::string name;
        int global_idx;

        int shock_idx;

        int refine_grad_max_lev;
        Vector<std::tuple<bool,int,Real, Real>> refine_grad_threshold; // {primitive?, idx, threshold, min value}

        int reflux;
        int project_divergence;

        BoundaryState boundary_conditions;

        bool mass_const, charge_const, gamma_const;
        Array<Real,2> mass, charge, gamma;

        int num_grow;

        Real dt; // time step restriction for this state (not the global time step)

        Real pressure_relaxation_rate;

        enum class DivIndex : int {
            DivB=0,
            DivD,
            NUM
        };

        // ratio of cleaning speed to fastest speed in domain
        Real relative_div_speed = 0.0;

        // the actual cleaning speed
        Real div_speed = 0.0;

        bool projection_cleaning = false;

        std::map<std::string,Optional3D1VFunction> functions;
        std::map<int,Optional3D1VFunction*> dynamic_functions; // index is for primitive values

        std::unique_ptr<Reconstruction> reconstruction;
        std::unique_ptr<RiemannSolver> flux_solver;
        std::unique_ptr<ShockDetector> shock_detector;
        std::unique_ptr<Viscous> viscous;

        Vector<std::pair<int,int>> associated_sources; // (system idx, src idx)
        std::map<AssociatedType,Vector<int>> associated_state; // (AssociatedType, {state idx, ...})

        Vector<Real> sum_cons;

        std::string particle_init;
        int particle_index;

        // initialisation routines
        virtual void init_from_lua();
        virtual void post_init_from_lua();
        virtual void set_udf();
        virtual void set_reconstruction();
        virtual void set_flux();
        virtual void set_reflux();
        virtual void set_viscosity();
        virtual void set_shock_detector();
        virtual void set_grad_refinement();
        //virtual void set_num_grow();
        virtual void set_num_grow(int n=0);

        // getters

        virtual const std::string& get_tag() const = 0;
        virtual int get_type() const {return +StateType::isNull;}
        virtual AssociatedType get_association_type() const {return AssociatedType::isNull;}
        virtual int n_prim() const {return 0;}
        virtual int n_cons() const {return 0;}
        virtual int n_flux() const {return 0;}

        virtual bool is_transported() const {return true;}
        virtual bool has_field() const {return false;}

        virtual bool get_enforce_positivity() const {return false;}
        virtual Vector<int> get_enforce_positivity_idx() const {return {-1};}

        int get_idx() const {return global_idx;}
        int get_shock_idx() const {return shock_idx;}

        virtual const Vector<std::string>& get_cons_names() const = 0;
        virtual const Vector<std::string>& get_prim_names() const = 0;
        virtual const Vector<std::pair<int,Real>> init_with_value() const {return {};}
        virtual const Vector<set_bc>& get_bc_set() const = 0;

        virtual Real get_mass(Real alpha) const {return 0.0;}
        virtual Real get_mass(const Vector<Real> &U) const{return 0.0;}

        virtual Real get_charge(Real alpha) const{return 0.0;}
        virtual Real get_charge(const Vector<Real> &U) const{return 0.0;}

        virtual Real get_gamma(Real alpha) const{return 0.0;}
        virtual Real get_gamma(const Vector<Real> &U) const{return 0.0;}

        virtual Real get_cp(Real alpha) const{return 0.0;}
        virtual Real get_cp(const Vector<Real> &U) const{return 0.0;}

        virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const {return false;}
        virtual bool cons2prim(Vector<autodiff::dual>& U, Vector<autodiff::dual>& Q) const {return false;}

        virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const {}
        virtual bool prim_valid(Vector<Real>& Q) const {return false;}
        virtual bool cons_valid(Vector<Real>& U) const {return false;}

        virtual int get_cons_density_idx() const {return -1;}
        virtual int get_prim_density_idx() const {return -1;}

        virtual const Vector<int> get_prim_vector_idx() const {return {-1};}
        virtual const Vector<int> get_flux_vector_idx() const {return get_prim_vector_idx();}
        virtual const Vector<int> get_cons_vector_idx() const {return {-1};}

        virtual const Vector<int> get_prs_idx() const {return {-1};}
        virtual const Vector<int> get_nrg_idx() const {return {-1};}

        virtual const int get_cons_D_idx() const {return -1;}
        virtual const int get_cons_B_idx() const {return -1;}

        virtual const int get_cons_phi_idx() const {return -1;}
        virtual const int get_cons_psi_idx() const {return -1;}

        virtual int get_cons_tracer_idx() const {return -1;}
        virtual int get_prim_alpha_idx() const {return -1;}

        virtual void update_div_clean(const Real* dx) {}

        virtual Real get_density_from_cons(const Vector<Real> &U) const {return -1.0;}
        virtual Real get_density_from_prim(const Vector<Real> &Q) const {return -1.0;}
        virtual dual get_density_from_prim(const Vector<dual> &Q) const {return -1.0;}

        virtual Real get_alpha_from_cons(const Vector<Real> &U) const {return -1.0;}
        virtual Real get_alpha_from_prim(const Vector<Real> &Q) const {return -1.0;}
        virtual dual get_alpha_from_prim(const Vector<dual> &Q) const {return -1.0;}

        virtual Real get_energy_from_cons(const Vector<Real> &U) const {return -1.0;}

        virtual Real get_temperature_from_cons(const Vector<Real> &U) const {return -1.0;}
        virtual Real get_temperature_from_prim(const Vector<Real> &Q) const {return -1.0;}
        virtual dual get_temperature_from_prim(const Vector<dual> &Q) const {return -1.0;}

        virtual RealArray get_speed_from_cons(const Vector<Real> &U) const {return {AMREX_D_DECL(0,0,0)};}
        virtual RealArray get_speed_from_prim(const Vector<Real> &Q) const {return {AMREX_D_DECL(0,0,0)};}

        virtual RealArray get_current_from_cons(const Vector<Real> &U) const {return {AMREX_D_DECL(0,0,0)};}

        // functions that work on Fabs

        static void get_state_vector(const FArrayBox& U, const int i, const int j, const int k, Vector<Real> &S);

        virtual void get_state_values(const Box& box,
                                      const FArrayBox& src,
                                      std::map<std::string,FArrayBox>& out,
                                      Vector<std::string>& updated
                                      EB_OPTIONAL(,const FArrayBox& vfrac)
                                      ) const;

        virtual void calc_velocity(const Box& box,
                                   FArrayBox& src,
                                   FArrayBox& vel
                                   EB_OPTIONAL(,const EBCellFlagFab &flag)
                                   ) const;

        virtual void calc_primitives(const Box& box,
                                     FArrayBox& cons,
                                     FArrayBox &prim,
                                     const Real* dx,
                                     const Real t,
                                     const Real* prob_lo
                                     EB_OPTIONAL(,const FArrayBox& vfrac)
                                     ) const;

        virtual void calc_reconstruction(const Box& box,
                                         FArrayBox& prim,
                                         Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                         Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                         EB_OPTIONAL(,const EBCellFlagFab &flag)
                                         EB_OPTIONAL(,const FArrayBox &vfrac)
                                         ) const;

        static void calc_slope(const Box& box,
                               const FArrayBox& src,
                               FArrayBox &slope,
                               EB_OPTIONAL(const EBCellFlagFab &flag,)
                               const Real *dx,
                               int index,
                               int dim,
                               Reconstruction& reco);


        virtual void calc_time_averaged_faces(const Box& box,
                                              const FArrayBox &prim,
                                              Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                              Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                                              EB_OPTIONAL(const EBCellFlagFab &flag,)
                                              const Real* dx,
                                              Real dt) const;

        Vector<BoundaryInfo> get_bc_limits(const Box& box,
                                           const Geometry &geom) const;

#ifdef AMREX_USE_EB
        void update_eb_vfrac(const Geometry &geom,
                       FArrayBox& vfrac) const;

        void update_eb_flags(const Geometry &geom,
                       EBCellFlagFab& flags) const;
#endif

        virtual void update_boundary_cells(const Box& box,
                                           const Geometry &geom,
                                           FArrayBox &prim,
                                           EB_OPTIONAL(const FArrayBox &vfrac,)
                                           const Real time) const;

        virtual void face_bc(const int dir,
                             Box const& box,
                             const FArrayBox &src,
                             FArrayBox& dest,
                             const int numcomp,
                             Geometry const& geom,
                             EB_OPTIONAL(const EBCellFlagFab &flag,)
                             const Real time,
                             const bool do_all=false) const;

        virtual void update_face_prim(const Box& box, const Geometry &geom,
                                      Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                                      Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                                      EB_OPTIONAL(const EBCellFlagFab &flag,)
                                      const Real time,
                                      const bool do_all=false) const;

        virtual void load_state_for_flux(Vector<Array4<const Real> > &face,
                                                 int i, int j, int k, Vector<Real> &S) const;

        virtual void calc_fluxes(const Box& box,
                                 Vector<Array<FArrayBox, AMREX_SPACEDIM>> &r_lo,
                                 Vector<Array<FArrayBox, AMREX_SPACEDIM>> &r_hi,
                                 Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                 EB_OPTIONAL(const EBCellFlagFab& flag,)
                                 const Real *dx, const Real dt,
                                 FArrayBox *shock = nullptr) const;

        virtual void correct_face_prim(const Box& box,
                                       Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                                       Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                                       const Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                       EB_OPTIONAL(const EBCellFlagFab &flag,)
                                       const Real *dx,
                                       const Real dt) const;

        virtual void calc_viscous_fluxes(const Box& box,
                                         Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                         const Box& pbox,
                                         const Vector<FArrayBox> &prim,
                                         EB_OPTIONAL(const EBCellFlagFab& flag,)
                                         const Real* dx) const {return;}


        virtual void calc_divergence(const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     FArrayBox& du,
                                     const Real *dx,
                                     const Real dt) const;

        virtual void calc_divergence(const Box& box,
                                     FArrayBox &U,
                                     FArrayBox& du,
                                     EB_OPTIONAL(const EBCellFlagFab &flag,)
                                     const Real *dx) const;

#ifdef AMREX_USE_EB

        bool eb_all_regular;
        int eb_index;
        const EB2::IndexSpace* eb2_index;
        Vector<std::unique_ptr<BoundaryEB>> eb_bcs;

        virtual void set_eb_bc(const sol::table &bc_def);

        virtual void set_eb_divergence();

        static bool check_covered_stencil(Array4<const EBCellFlag> const& flag, int i, int j, int k, int d, int stencil_length);

        virtual Vector<Real> calc_wall_normal_slopes(const Vector<int> &slope_idx,
                                                     const Vector<Real> &wall_value,
                                                     Array<Real,3> &wall_normal,
                                                     Array<Real, AMREX_SPACEDIM> &wall_centre,
                                                     Array4<const Real> const &prim4,
                                                     const int i, const int j, const int k) const;

        virtual void calc_wall_fluxes(const Box& box,
                                      const Vector<FArrayBox>& prim,
                                      Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                      const EBCellFlagFab& flag,
                                      const CutFab &bc_idx,
                                      const FArrayBox& bcent,
                                      const FArrayBox& bnorm,
                                      const Array<const FArrayBox*, AMREX_SPACEDIM> &afrac,
                                      const Real *dx,
                                      const Real dt) const;

        virtual void get_wall_value(const Box& box,
                                    std::map<BoundaryEB::EBType,FArrayBox*> bcs_data,
                                    const EBCellFlagFab& flag,
                                    const CutFab &bc_idx,
                                    const FArrayBox &bcent,
                                    const FArrayBox &bnorm,
                                    const Real t,
                                    const Real *dx,
                                    const Real *prob_lo) const;

        std::unique_ptr<DivergenceEB> eb_div;

#endif

        virtual void write_info(nlohmann::json& js) const;

        virtual std::string str() const;
};

template <typename D>
std::unique_ptr<State> StateBuilder(const sol::table& def)
{
    std::string type = def["type"];
    if (type == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

PhysicsFactory<State>& GetStateFactory();

#endif // MFP_GENERIC_H
