#ifndef MFP_HYDRO_VISCOUS_H
#define MFP_HYDRO_VISCOUS_H

#include "MFP.H"
#include "MFP_hydro_defs.H"
#include "MFP_factory.H"

//-------------------------------------------------------------------------


class HydroViscous
{
public:
    HydroViscous();
    virtual ~HydroViscous() = default;

    enum class CoeffIdx {
        Temp=0,
        Kappa,
        Mu,
        NUM
    };

    virtual void get_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa) const {}

    virtual std::string get_tag() const {return "";}
    virtual const std::map<std::string, Real> get_refs() const {return {};}
    virtual std::string str() const;

    Real get_min_dt(MFP *mfp) const;
    virtual Real calc_stability_condition(Vector<Real>& U) const {return 0.0;}

    int idx;
    bool active;
};

template <typename D>
std::unique_ptr<HydroViscous> HydroViscousBuilder(const sol::table& def)
{
    std::unique_ptr<HydroViscous> viscous;

    const int global_idx = def["global_idx"];
    const sol::table& visc = def["viscosity"].get_or(sol::table());

    if (visc.valid()) {

        // load the diffusion coefficients

        if (visc["type"] == D::tag) {
            return std::unique_ptr<D>(new D(global_idx, visc));
        } else {
            return nullptr;
        }
    } else {
        return nullptr;
    }
}

ClassFactory<HydroViscous>& GetHydroViscousFactory();

// ====================================================================================

class Sutherland : public HydroViscous
{
public:
    Sutherland();
    Sutherland(const int global_idx, const sol::table &def);
    ~Sutherland();

    Real mu_0, T0, S, Prandtl;

    virtual std::string get_tag() const override {return tag;}
    virtual const std::map<std::string, Real> get_refs() const override {
        return {{"mu_0", mu_0},{"T0", T0}, {"S",S}, {"Pr",Prandtl}};
    }

    virtual Real calc_stability_condition(Vector<Real> &U) const override;
    virtual void get_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa) const override;

    static bool registered;
    static std::string tag;
};

// ====================================================================================

class PowerLaw : public HydroViscous
{
public:
    PowerLaw();
    PowerLaw(const int global_idx, const sol::table& def);
    ~PowerLaw();

    enum NeutralDiffusiveVariable {
        Temp = 0,
        Mu,
        Kappa,
        NUM_NEUTRAL_DIFF
    };

    Real mu_0, T0, n, Prandtl;

    virtual std::string get_tag() const override {return tag;}
    virtual const std::map<std::string, Real> get_refs() const override {
        return {{"mu_0", mu_0},{"T0", T0}, {"n",n}, {"Pr",Prandtl}};
    }

    virtual Real calc_stability_condition(Vector<Real> &U) const override;
    virtual void get_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa) const override;

    static bool registered;
    static std::string tag;
};

// ====================================================================================

class UserDefinedViscosity : public HydroViscous
{
public:
    UserDefinedViscosity();
    UserDefinedViscosity(const int global_idx, const sol::table& def);
    ~UserDefinedViscosity();

    enum NeutralDiffusiveVariable {
        Temp = 0,
        Mu,
        Kappa,
        NUM_NEUTRAL_DIFF
    };

    Real mu_0, T0, n, Prandtl;

    virtual std::string get_tag() const override {return tag;}
    virtual const std::map<std::string, Real> get_refs() const override {
        return {{"mu_0", mu_0}, {"Pr",Prandtl}};
    }

    virtual Real calc_stability_condition(Vector<Real> &U) const override;
    virtual void get_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa) const override;

    static bool registered;
    static std::string tag;
};

// ====================================================================================

#endif // MFP_HYDRO_VISCOUS_H
