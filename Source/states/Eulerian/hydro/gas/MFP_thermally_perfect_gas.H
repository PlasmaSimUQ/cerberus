#ifndef THERMALLYPERFECTGAS_H
#define THERMALLYPERFECTGAS_H

#include "MFP_hydro_gas.H"


class ThermallyPerfectGas  : public HydroGas
{
public:
    ThermallyPerfectGas();
    ThermallyPerfectGas(const int global_idx, const sol::table &def);
    ~ThermallyPerfectGas();

    virtual Real get_gamma_from_cons(const Vector<Real> &U, const int density_idx=+HydroDef::ConsIdx::Density, const int tracer_idx=+HydroDef::ConsIdx::NUM) const override;
    virtual Real get_cp_from_cons(const Vector<Real> &U, const int density_idx=+HydroDef::ConsIdx::Density, const int tracer_idx=+HydroDef::ConsIdx::NUM) const override;

    virtual Real get_gamma_from_prim(const Vector<Real> &Q, const int idx=+HydroDef::PrimIdx::NUM) const override;
    virtual Real get_cp_from_prim(const Vector<Real> &Q, const int tracer_idx=+HydroDef::PrimIdx::NUM) const override;

    virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const override;
    virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const override;
    virtual void define_rho_p_T(Vector<Real>& Q) const override;

    virtual Real get_temperature_from_cons(const Vector<Real>& U) const override;
    virtual RealArray get_speed_from_cons(const Vector<Real>& U) const override;
    virtual RealArray get_speed_from_prim(const Vector<Real>& Q) const override;

    virtual std::string get_tag() const override {return tag;}
    virtual void write_info(nlohmann::json &js) const override;

    static bool registered;
    static std::string tag;

    bool gamma_const;
    Vector<Real> gamma;
};

#endif // THERMALLYPERFECTGAS_H
