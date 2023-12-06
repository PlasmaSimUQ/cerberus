#ifndef HYDROSHOCKDETECTOR_H
#define HYDROSHOCKDETECTOR_H

#include "MFP_eulerian_riemann.H"
#include "MFP_factory.H"
#include "MFP_hydro_defs.H"
#include "json.hpp"
#include "sol.hpp"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

using namespace amrex;

class HydroShockDetector : public ShockDetector
{
  public:
    HydroShockDetector();
    ~HydroShockDetector();
};

template <typename D>
std::unique_ptr<HydroShockDetector> HydroShockDetectorBuilder(const sol::table& def)
{
    //    int global_idx = def["global_idx"];

    if (def["name"] == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

ClassFactory<HydroShockDetector>& GetHydroShockDetectorFactory();

//-----------------------------------------------------------------------------

class PressureJumpShockDetector : public HydroShockDetector
{
  public:
    PressureJumpShockDetector();
    PressureJumpShockDetector(const sol::table& def);

    virtual Real solve(Vector<Real>& L, Vector<Real>& R) const override;

    virtual const std::string get_tag() const override { return tag; }
    virtual void write_info(nlohmann::json& js) const override;

    static std::string tag;
    static bool registered;

  private:
    Real shock_threshold;
};

#endif  // SHOCKDETECTOR_H
