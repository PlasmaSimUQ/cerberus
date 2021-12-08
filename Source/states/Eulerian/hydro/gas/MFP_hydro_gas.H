#ifndef MFP_HYDRO_GAS_H
#define MFP_HYDRO_GAS_H

#include "MFP.H"
#include "MFP_hydro_defs.H"
#include "MFP_factory.H"

//-------------------------------------------------------------------------


class HydroGas
{
public:
    HydroGas();
    ~HydroGas();

    size_t n_species() const {return mass.size();}
    size_t n_tracers() const {return n_species() - 1;}

    virtual void get_alpha_fractions_from_prim(const Vector<Real> &Q, Vector<Real> &alpha, const int tracer_idx=+HydroDef::PrimIdx::NUM) const {}
    virtual void get_alpha_fractions_from_cons(const Vector<Real> &U, Vector<Real> &alpha, const int density_idx=+HydroDef::ConsIdx::Density, const int tracer_idx=+HydroDef::ConsIdx::NUM) const {}

    virtual Real get_mass_from_cons(const Vector<Real> &U, const int density_idx=+HydroDef::ConsIdx::Density, const int tracer_idx=+HydroDef::ConsIdx::NUM) const {return 0;}
    virtual Real get_charge_from_cons(const Vector<Real> &U, const int density_idx=+HydroDef::ConsIdx::Density, const int tracer_idx=+HydroDef::ConsIdx::NUM) const {return 0;}
    virtual Real get_gamma_from_cons(const Vector<Real> &U, const int density_idx=+HydroDef::ConsIdx::Density, const int tracer_idx=+HydroDef::ConsIdx::NUM) const {return 0;}
    virtual Real get_cp_from_cons(const Vector<Real> &U, const int density_idx=+HydroDef::ConsIdx::Density, const int tracer_idx=+HydroDef::ConsIdx::NUM) const {return 0;}

    virtual Real get_mass_from_prim(const Vector<Real> &Q, const int idx=+HydroDef::PrimIdx::NUM) const {return 0;}
    virtual Real get_charge_from_prim(const Vector<Real> &Q, const int tracer_idx=+HydroDef::PrimIdx::NUM) const {return 0;}
    virtual Real get_gamma_from_prim(const Vector<Real> &Q, const int idx=+HydroDef::PrimIdx::NUM) const {return 0;}
    virtual Real get_cp_from_prim(const Vector<Real> &Q, const int tracer_idx=+HydroDef::PrimIdx::NUM) const {return 0;}

    virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const {return true;}
    virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const {}

    virtual void define_rho_p_T(Vector<Real>& Q) const {}
    Real get_density_from_number_density(const Real nd, Vector<Real>& Q) const;

    virtual bool prim_valid(const Vector<Real>& Q) const {return true;}
    virtual bool cons_valid(const Vector<Real>& U) const {return true;}
    virtual Real get_energy_from_cons(const Vector<Real>& U) const {return 0;}
    virtual Real get_internal_energy_density_from_cons(const Vector<Real>& U) const {return 0;}
    virtual Real get_temperature_from_cons(const Vector<Real>& U) const {return 0;}
    virtual Real get_temperature_from_prim(const Vector<Real>& Q) const {return 0;}
    virtual RealArray get_speed_from_cons(const Vector<Real>& U) const {return {AMREX_D_DECL(0,0,0)};}
    virtual RealArray get_speed_from_prim(const Vector<Real>& Q) const {return {AMREX_D_DECL(0,0,0)};}


    virtual std::string get_tag() const {return "";}
    virtual void write_info(nlohmann::json &js) const;

    int idx;
    bool mass_const, charge_const;
    Vector<Real> mass, charge;
    Vector<std::string> comp_names;
};

template <typename D>
std::unique_ptr<HydroGas> HydroGasBuilder(const sol::table& def)
{
    const int global_idx = def["global_idx"];
    const sol::table& gas = def["gas"].get_or(sol::table());

    if (gas.valid()) {

        // load the diffusion coefficients

        if (gas["type"] == D::tag) {
            return std::unique_ptr<D>(new D(global_idx, gas));
        } else {
            return nullptr;
        }
    } else {
        return nullptr;
    }
}

ClassFactory<HydroGas>& GetHydroGasFactory();


#endif // MFP_HYDRO_GAS_H
