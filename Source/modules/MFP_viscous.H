#ifndef MFP_VISCOUS_H
#define MFP_VISCOUS_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include "MFP_factory.H"
#include "sol.hpp"

using amrex::Real;
using amrex::Vector;
using amrex::Array;
// ====================================================================================

class Viscous
{
public:
  Viscous();
  ~Viscous();

  enum DiffusionType {
    Neutral = 0,
    Ion,
    Electron
  };

  enum NeutralDiffusionCoeffs {
    NeutralTemp = 0,
    NeutralMu,
    NeutralKappa,
    NUM_NEUTRAL_DIFF_COEFFS
  };

  enum IonDiffusionCoeffs {
    NUM_ION_DIFF_COEFFS
  };

  enum ElectronDiffusionCoeffs {
    NUM_ELE_DIFF_COEFFS
  };

  virtual int get_type();
  virtual int get_num();
  virtual void get_neutral_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa);
  virtual void get_ion_coeffs(){return;}
  virtual void get_electron_coeffs(){return;}
  virtual void update_linked_states();
  virtual Vector<int> get_linked_states() {return {idx};}
  virtual Real get_max_speed(const Vector<Vector<Real>> &U);
  
  virtual const std::string& get_tag() const = 0;
  virtual const std::map<std::string, Real> get_refs() const {return {};}
  virtual std::string str() const;

  int idx;
  Real cfl;
};

// ====================================================================================

class Sutherland : public Viscous
{
public:
  Sutherland();
  Sutherland(const int global_idx, const sol::table &def);
  ~Sutherland();

  Real mu_0, T0, S, Prandtl;

  virtual int get_type() override;
  virtual int get_num() override;
  virtual void get_neutral_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa) override;
  virtual Vector<int> get_linked_states() override {return {idx};}
  virtual Real get_max_speed(const Vector<Vector<Real>> &U) override;
  static bool valid_state(const int idx);
  virtual const std::string& get_tag() const override {return tag;}
  virtual const std::map<std::string, Real> get_refs() const override {
      return {{"mu_0", mu_0},{"T0", T0}, {"S",S}, {"Pr",Prandtl}};
  }
  static bool registered;
  static std::string tag;
};

// ====================================================================================

class PowerLaw : public Viscous
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

  virtual int get_type() override;
  virtual int get_num() override;
  virtual void get_neutral_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa) override;
  virtual Vector<int> get_linked_states() override {return {idx};}
  virtual Real get_max_speed(const Vector<Vector<Real>> &U) override;
  static bool valid_state(const int idx);
  virtual const std::string& get_tag() const override {return tag;}
  virtual const std::map<std::string, Real> get_refs() const override {
      return {{"mu_0", mu_0},{"T0", T0}, {"n",n}, {"Pr",Prandtl}};
  }
  static bool registered;
  static std::string tag;
};

// ====================================================================================

class UserDefinedViscosity : public Viscous
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

  virtual int get_type() override;
  virtual int get_num() override;
  virtual void get_neutral_coeffs(const Vector<Real> &Q, Real &T, Real &mu, Real &kappa) override;
  virtual Real get_max_speed(const Vector<Vector<Real>>& U) override;
  static bool valid_state(const int idx);
  virtual const std::string& get_tag() const override {return tag;}
  virtual const std::map<std::string, Real> get_refs() const override {
      return {{"mu_0", mu_0}, {"Pr",Prandtl}};
  }
  static bool registered;
  static std::string tag;
};

// ====================================================================================


template <typename D>
std::unique_ptr<Viscous> ViscousBuilder(const sol::table& def)
{
    std::unique_ptr<Viscous> viscous;

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



PhysicsFactory<Viscous>& GetViscousFactory();

#endif // MFP_DIFFUSION_H
