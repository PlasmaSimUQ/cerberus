#ifndef MHDSTATE_H
#define MHDSTATE_H
#include "forward.hpp"

class MHDState;

#include "MFP_eulerian.H"
#include "MFP_bc.H"
#include "MFP.H"
#include "MFP_mhd_defs.H"
#include "MFP_reconstruction.H"
#include "MFP_mhd_riemann.H"
#include "MFP_mhd_bc.H"
#include "MFP_mhd_shockdetector.H"


class MHDState : public EulerianState
{
public:
    MHDState();
    MHDState(const sol::table& def);
    ~MHDState();

    static bool registered;
    static std::string tag;

    static Vector<std::string> cons_names;
    static Vector<std::string> prim_names;
    static Array<int,2> cons_vector_idx;
    static Array<int,2> prim_vector_idx;
    static std::map<std::string, int> bc_names;
    static Vector<set_bc> bc_set;

    Real gamma;

    Real div_transport;

    std::unique_ptr<MHDRiemannSolver> flux_solver;
    std::unique_ptr<MHDShockDetector> shock_detector;

    virtual StateType get_type() const override {return StateType::MHD;}

    size_t n_cons() const {return +MHDDef::ConsIdx::NUM;}
    size_t n_prim() const {return +MHDDef::PrimIdx::NUM;}

#ifdef AMREX_USE_EB
    Vector<std::unique_ptr<MHDBoundaryEB>> eb_bcs;
    virtual void set_eb_bc(const sol::table &bc_def) override;
    virtual size_t get_eb_bc_size() const override {return eb_bcs.size();}
#endif

    virtual int get_num_grow() const override {return num_grow;}
    virtual void set_num_grow(int n) override {num_grow = std::max(n, num_grow);}

    virtual int get_num_cons() const override {return +MHDDef::ConsIdx::NUM;}
    virtual const Vector<std::string>& get_cons_names() const override {return cons_names;}

    virtual Vector<std::string> get_plot_output_names() const override;
    virtual void get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string,FArrayBox>& out,
                                 Vector<std::string>& updated
                             #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                             #endif
                                 ) const override;

    void set_udf();

    void set_flux();
    void set_shock_detector();
    void set_refinement();

    virtual void variable_setup(Vector<int> periodic) override;
    virtual void init_data(MFP* mfp, const Real time) override;
    virtual void init_from_lua() override;

    bool cons2prim(Array<Real, +MHDDef::ConsIdx::NUM>& U, Array<Real, +MHDDef::PrimIdx::NUM> &Q) const;
    void prim2cons(Array<Real, +MHDDef::PrimIdx::NUM> &Q, Array<Real, +MHDDef::ConsIdx::NUM> &U) const;

    bool prim_valid(const Array<Real, +MHDDef::PrimIdx::NUM>& Q) const;
    bool cons_valid(const Array<Real, +MHDDef::ConsIdx::NUM>& U) const;

    RealArray get_speed_from_cons(const Array<Real, +MHDDef::ConsIdx::NUM> &U) const;
    RealArray get_speed_from_prim(const Array<Real,+MHDDef::PrimIdx::NUM>& Q) const;

    virtual Real get_allowed_time_step(MFP* mfp) const override;

    void calc_velocity(const Box& box,
                       FArrayBox& cons,
                       FArrayBox &prim
                   #ifdef AMREX_USE_EB
                       ,const FArrayBox& vfrac
                   #endif
                       ) const;

    virtual void calc_primitives(const Box& box,
                                 FArrayBox& cons,
                                 FArrayBox &prim,
                                 const Real* dx,
                                 const Real t,
                                 const Real* prob_lo
                             #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                             #endif
                                 ) const override;

    virtual void update_boundary_cells(const Box& box,
                                       const Geometry &geom,
                                       FArrayBox &prim,
                                   #ifdef AMREX_USE_EB
                                       const FArrayBox& vfrac,
                                   #endif
                                       const Real time) const override;

    virtual void calc_reconstruction(const Box& box,
                                     FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                 #ifdef AMREX_USE_EB
                                     ,const EBCellFlagFab &flag
                                     ,const FArrayBox &vfrac
                                 #endif
                                     ) const override;

    virtual void calc_time_averaged_faces(const Box& box,
                                          const FArrayBox &prim,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                                      #ifdef AMREX_USE_EB
                                          const EBCellFlagFab& flag,
                                      #endif
                                          const Real* dx,
                                          Real dt) const override;

    virtual void face_bc(const int dir,
                         Box const& box,
                         const FArrayBox& src,
                         FArrayBox& dest,
                         const Geometry &geom,
                     #ifdef AMREX_USE_EB
                         const EBCellFlagFab& flag,
                     #endif
                         const Real time,
                         const bool do_all) const override;

    void load_state_for_flux(const Array4<const Real> &face,
                             int i, int j, int k,
                             Array<Real,+MHDDef::PrimIdx::NUM> &S) const;

    virtual void calc_fluxes(const Box& box,
                             FArrayBox &cons,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                             Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                         #ifdef AMREX_USE_EB
                             const EBCellFlagFab& flag,
                         #endif
                             const Real *dx,
                             const Real dt) const override;

    virtual void correct_face_prim(const Box& box,
                                   Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                                   Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                                   const Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                               #ifdef AMREX_USE_EB
                                   const EBCellFlagFab &flag,
                               #endif
                                   const Real *dx,
                                   const Real dt) const override;



#ifdef AMREX_USE_EB


    virtual void calc_wall_fluxes(const Box& box,
                                  const FArrayBox &prim,
                                  Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                  const EBCellFlagFab& flag,
                                  const CutFab &bc_idx,
                                  const FArrayBox& bcent,
                                  const FArrayBox &bnorm,
                                  const Array<const FArrayBox*, AMREX_SPACEDIM> &afrac,
                                  const Real *dx,
                                  const Real dt) const override;
#endif


    virtual void write_info(nlohmann::json& js) const override;

    static MHDState& get_state_global(const int idx) {return static_cast<MHDState&>(MFP::get_state(idx));}
    static MHDState& get_state(const std::string& name) {return static_cast<MHDState&>(MFP::get_state(name));}

};

#endif // MHDSTATE_H
