#ifndef MFP_HYDRO_H
#define MFP_HYDRO_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include "MFP_bc.H"
#include "MFP_utility.H"
#include "MFP_state.H"

using namespace amrex;

struct HydroState : public State
{
public:
    HydroState();
    HydroState(const sol::table& def);
    ~HydroState();

    enum class ConsIdx : int {
        Density=0,
        Xmom,
        Ymom,
        Zmom,
        Eden,
        Tracer,
        NUM,
    };

    enum class PrimIdx : int {
        Density=0,
        Xvel,
        Yvel,
        Zvel,
        Prs,
        Temp,
        Alpha,
        NUM,
    };

    const std::map<std::string, int> bc_names = {{"interior",  PhysBCType::interior},
                                                 {"inflow",    PhysBCType::inflow},
                                                 {"outflow",   PhysBCType::outflow},
                                                 {"symmetry",  PhysBCType::symmetry},
                                                 {"slipwall",  PhysBCType::slipwall},
                                                 {"noslipwall",PhysBCType::noslipwall}};

    static bool registered;
    static std::string tag;

    static Vector<std::string> cons_names;
    static Vector<std::string> prim_names;
    static Vector<set_bc> bc_set;

    bool enforce_positivity = false;
    const Vector<int> enforce_positivity_idx = {+ConsIdx::Density, +ConsIdx::Eden};

    bool extra_slope_limits = true;

    virtual void init_from_lua() override;
#ifdef AMREX_USE_EB
    virtual void set_eb_bc(const sol::table &bc_def) override;
#endif

    Real init_from_number_density(std::map<std::string, Real> data);
    virtual void set_udf() override;

    virtual const std::string& get_tag() const override {return tag;}
    virtual int get_type() const override {return +StateType::isHydro;}
    virtual AssociatedType get_association_type() const override;
    virtual int n_prim() const override {return +PrimIdx::NUM;}
    virtual int n_cons() const override {return +ConsIdx::NUM;}

    virtual bool get_enforce_positivity() const override {return enforce_positivity;}
    virtual Vector<int> get_enforce_positivity_idx() const override {return enforce_positivity_idx;}

    virtual const Vector<std::string>& get_cons_names() const override;
    virtual const Vector<std::string>& get_prim_names() const override;
    virtual const Vector<set_bc>& get_bc_set() const override;

    virtual Real get_mass(Real alpha) const override;
    virtual Real get_mass(const Vector<Real> &U) const override;

    virtual Real get_charge(Real alpha) const override;
    virtual Real get_charge(const Vector<Real> &U) const override;

    virtual Real get_gamma(Real alpha) const override;
    virtual Real get_gamma(const Vector<Real> &U) const override;

    virtual Real get_cp(Real alpha) const override;
    virtual Real get_cp(const Vector<Real> &U) const override;

    virtual bool cons2prim(Vector<Real>& U) const override;
    virtual bool cons2prim(Vector<autodiff::dual>& U) const override;

    virtual void prim2cons(Vector<Real>& U) const override;
    virtual bool prim_valid(Vector<Real>& Q) const override;
    virtual bool cons_valid(Vector<Real>& U) const override;

    virtual int get_cons_density_idx() const override {
        return +ConsIdx::Density;
    }
    virtual const Vector<int> get_prim_vector_idx() const override {
        return {+PrimIdx::Xvel};
    }
    virtual const Vector<int> get_cons_vector_idx() const override {
        return {+ConsIdx::Xmom};
    }

    virtual const Vector<int> get_prs_idx() const override {
        return {+PrimIdx::Prs};
    }
    virtual const Vector<int> get_nrg_idx() const override {
        return {+ConsIdx::Eden};
    }

    virtual Real get_density_from_cons(const Vector<Real> &U) const override
    {
        return U[+ConsIdx::Density];
    }
    virtual Real get_density_from_prim(const Vector<Real> &Q) const override
    {
        return Q[+PrimIdx::Density];
    }
    virtual dual get_density_from_prim(const Vector<dual> &Q) const override
    {
        return Q[+PrimIdx::Density];
    }

    virtual int get_cons_tracer_idx() const override {return +ConsIdx::Tracer;}
    virtual int get_prim_alpha_idx() const override {return +PrimIdx::Alpha;}

    virtual Real get_alpha_from_cons(const Vector<Real> &U) const override
    {
        return U[+ConsIdx::Tracer]/U[+ConsIdx::Density];
    }
    virtual Real get_alpha_from_prim(const Vector<Real> &Q) const override
    {
        return Q[+PrimIdx::Alpha];
    }
    virtual dual get_alpha_from_prim(const Vector<dual> &Q) const override
    {
        return Q[+PrimIdx::Alpha];
    }

    virtual Real get_energy_from_cons(const Vector<Real> &U) const override;

    virtual Real get_temperature_from_cons(const Vector<Real> &U) const override;
    virtual Real get_temperature_from_prim(const Vector<Real> &Q) const override;
    virtual dual get_temperature_from_prim(const Vector<dual> &Q) const override;


    virtual RealArray get_speed_from_cons(const Vector<Real> &U) const override;
    virtual RealArray get_speed_from_prim(const Vector<Real> &Q) const override;

    virtual RealArray get_current_from_cons(const Vector<Real> &U) const override;

    void calc_reconstruction(const Box& box,
                             FArrayBox &prim,
                             Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                             Array<FArrayBox, AMREX_SPACEDIM> &rhi
                             EB_OPTIONAL(,const EBCellFlagFab& flag)
                             EB_OPTIONAL(,const FArrayBox &vfrac)
                             ) const override;

    virtual Real local_shock_detector(const Vector<Real> &L,
                                      const Vector<Real> &R) const override;

    void get_state_values(const Box& box,
                                  const FArrayBox& src,
                                  std::map<std::string,FArrayBox>& out,
                                  Vector<std::string>& updated
                                  EB_OPTIONAL(,const FArrayBox& vfrac)
                                  ) const override;


    void calc_velocity(const Box& box,
                       FArrayBox& src,
                       FArrayBox& vel
                       EB_OPTIONAL(,const EBCellFlagFab& flag)
                       ) const override;

    virtual void calc_viscous_fluxes(const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const Box& pbox,
                                     const Vector<FArrayBox> &prim,
                                     EB_OPTIONAL(const EBCellFlagFab &flag,)
                                     const Real* dx) const override;

    void calc_neutral_diffusion_terms(const Box& box,
                                      const FArrayBox &prim,
                                      FArrayBox& diff
                                      EB_OPTIONAL(,const EBCellFlagFab& flag)
                                      ) const;

    void calc_neutral_viscous_fluxes(const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const Box &pbox,
                                     const FArrayBox &prim,
                                     EB_OPTIONAL(const EBCellFlagFab &flag,)
                                     const Real* dx) const;

    void calc_charged_viscous_fluxes(int passed_idx,
                                     int ion_idx,
                                     int electron_idx,
                                     int em_idx,
                                     const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const Box& pbox,
                                     const Vector<FArrayBox>& prim,
                                     EB_OPTIONAL(const EBCellFlagFab &flag,)
                                     const Real* dx,
                                     FArrayBox& diff) const;

    void calc_ion_diffusion_terms(const Box& box,const Vector<FArrayBox>& prim,
                                  State& EMstate,Array4<const Real> const& prim_EM4,
                                  State& ELEstate,Array4<const Real> const& prim_ELE4,
                                  FArrayBox& diff
                                  EB_OPTIONAL(,const EBCellFlagFab& flag)
                                  ) const ;

    void calc_ion_viscous_fluxes(const Box& box,
                                 Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                 const Box &pbox,
                                 const Vector<FArrayBox>& prim,
                                 EB_OPTIONAL(const EBCellFlagFab& flag,)
                                 const Real* dx) const;

    void calc_electron_diffusion_terms(const Box& box,const Vector<FArrayBox>& prim,
                                       State& EMstate,Array4<const Real> const& prim_EM4,
                                       State& ELEstate,Array4<const Real> const& prim_ELE4,
                                       FArrayBox& diff
                                       EB_OPTIONAL(,const EBCellFlagFab& flag)
                                       ) const ;
    
    void calc_electron_viscous_fluxes(const Box& box, Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                      const Box &pbox, const Vector<FArrayBox>& prim,
                                      EB_OPTIONAL(const EBCellFlagFab& flag,)
                                      const Real* dx) const;

#ifdef AMREX_USE_EB
    void calc_neutral_viscous_fluxes_eb(const Box& box, Array<FArrayBox,
                                        AMREX_SPACEDIM> &fluxes,
                                        const Box& pbox,
                                        const FArrayBox &prim,
                                        const EBCellFlagFab& flag,
                                        const Real* dx) const;
#endif

    virtual void write_info(nlohmann::json &js) const override;

    virtual std::string str() const override;
};

#endif // MFP_HYDRO_H
