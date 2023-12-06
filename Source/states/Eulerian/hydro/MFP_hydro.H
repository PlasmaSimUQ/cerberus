#ifndef FLUIDSTATE_H
#define FLUIDSTATE_H
#include "forward.hpp"

class HydroState;

#include "MFP.H"
#include "MFP_bc.H"
#include "MFP_eulerian.H"
#include "MFP_hydro_bc.H"
#include "MFP_hydro_defs.H"
#include "MFP_hydro_gas.H"
#include "MFP_hydro_riemann.H"
#include "MFP_hydro_shockdetector.H"
#include "MFP_hydro_viscous.H"
#include "MFP_reconstruction.H"

class HydroState : public EulerianState
{
  public:
    HydroState();
    HydroState(const sol::table& def);
    ~HydroState();

    static bool registered;
    static std::string tag;

    static Vector<std::string> cons_names;
    static Vector<std::string> prim_names;
    static Array<int, 1> flux_vector_idx;
    static Array<int, 1> cons_vector_idx;
    static Array<int, 1> prim_vector_idx;
    static std::map<std::string, int> bc_names;
    static Vector<set_bc> bc_set;

    int n_tracers;
    static std::string multicomp_prim_name;
    static std::string multicomp_cons_name;
    std::string get_multicomp_name(const std::string& name, const int i) const
    {
        return name + "_" + num2str(i);
    }

    std::unique_ptr<HydroGas> gas;
    std::unique_ptr<HydroViscous> viscous;

#ifdef AMREX_USE_EB
    virtual void set_eb_bc(const sol::table& bc_def) override;
#endif

    virtual StateType get_type() const override { return StateType::Hydro; }

    virtual size_t n_cons() const override { return +HydroDef::ConsIdx::NUM + n_tracers; }
    virtual size_t n_prim() const override { return +HydroDef::PrimIdx::NUM + n_tracers; }

    virtual const Vector<int> get_cons_vector_idx() const override
    {
        return arr2vec(cons_vector_idx);
    }
    virtual const Vector<int> get_prim_vector_idx() const override
    {
        return arr2vec(prim_vector_idx);
    }

    virtual bool is_viscous() const { return viscous != nullptr; }

    virtual int get_num_grow() const override { return num_grow; }
    virtual void set_num_grow(int n) override { num_grow = std::max(n, num_grow); }

    virtual const Vector<std::string>& get_cons_names() const override { return cons_names; }
    virtual const Vector<std::string>& get_prim_names() const override { return prim_names; }

    virtual Vector<std::string> get_plot_output_names() const override;
    virtual void get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string, FArrayBox>& out,
                                 Vector<std::string>& updated
#ifdef AMREX_USE_EB
                                 ,
                                 const FArrayBox& vfrac
#endif
    ) const override;

    void set_gas();
    void set_viscosity();

    Real init_from_number_density(std::map<std::string, Real> data);
    void set_udf();

    void set_flux();
    void set_shock_detector();
    void set_refinement();

    virtual void variable_setup(Vector<int> periodic) override;
    virtual void init_data(MFP* mfp, const Real time) override;
    virtual void init_from_lua() override;

    virtual Real get_allowed_time_step(MFP* mfp) const override;

    virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const { return gas->cons2prim(U, Q); }
    virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const { return gas->prim2cons(Q, U); }

    virtual RealArray get_speed_from_cons(const Vector<Real>& U) const override
    {
        return gas->get_speed_from_cons(U);
    }
    virtual RealArray get_speed_from_prim(const Vector<Real>& Q) const override
    {
        return gas->get_speed_from_prim(Q);
    }

    void calc_velocity(const Box& box,
                       FArrayBox& cons,
                       FArrayBox& prim
#ifdef AMREX_USE_EB
                       ,
                       const FArrayBox& vfrac
#endif
    ) const;

    virtual void calc_reconstruction(const Box& box,
                                     FArrayBox& prim,
                                     Array<FArrayBox, AMREX_SPACEDIM>& rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM>& rhi
#ifdef AMREX_USE_EB
                                     ,
                                     const EBCellFlagFab& flag,
                                     const FArrayBox& vfrac
#endif
    ) const override;

    void calc_diffusion_terms(const FArrayBox& prim,
                              FArrayBox& diff
#ifdef AMREX_USE_EB
                              ,
                              const EBCellFlagFab& flag
#endif
    ) const;

    virtual void calc_viscous_fluxes(const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM>& fluxes,
                                     const FArrayBox& prim,
#ifdef AMREX_USE_EB
                                     const EBCellFlagFab& flag,
#endif
                                     const Real* dx) const override;

#ifdef AMREX_USE_EB

    virtual void calc_viscous_fluxes_eb(const Box& box,
                                        Array<FArrayBox, AMREX_SPACEDIM>& fluxes,
                                        const FArrayBox& prim,
    #ifdef AMREX_USE_EB
                                        const EBCellFlagFab& flag,
    #endif
                                        const Real* dx) const;

#endif

    void calc_current_and_charge(const Box& box,
                                 const FArrayBox& cons,
                                 FArrayBox* cd,
                                 FArrayBox* J
#ifdef AMREX_USE_EB
                                 ,
                                 const FArrayBox& vfrac
#endif
    ) const;

    virtual void write_info(nlohmann::json& js) const override;

    static HydroState& get_state_global(const int idx)
    {
        return static_cast<HydroState&>(MFP::get_state(idx));
    }
    static HydroState& get_state(const std::string& name)
    {
        return static_cast<HydroState&>(MFP::get_state(name));
    }
};

#endif  // FLUIDSTATE_H
