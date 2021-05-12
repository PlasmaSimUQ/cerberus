#ifndef SHOCKDETECTOR_H
#define SHOCKDETECTOR_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "sol.hpp"
#include "json.hpp"
#include "MFP_factory.H"

using namespace amrex;

class State;


class ShockDetector
{
public:
    ShockDetector();
    ~ShockDetector();
    virtual Real solve(Vector<Real> &L,Vector<Real> &R) const;

    int idx;
    State* istate;

    virtual const std::string& get_tag() const = 0;
    virtual void write_info(nlohmann::json& js) const {return;}
};

template <typename D>
std::unique_ptr<ShockDetector> ShockDetectorBuilder(const sol::table& def)
{
    int global_idx = def["global_idx"];

    if ((def["name"] == D::tag) && (D::valid_state(global_idx))) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

PhysicsFactory<ShockDetector>& GetShockDetectorFactory();

//-----------------------------------------------------------------------------

class HydroShockDetector : public ShockDetector
{
public:
    HydroShockDetector();
    HydroShockDetector(const sol::table &def);

    virtual Real solve(Vector<Real> &L,
                       Vector<Real> &R) const override;


    static bool valid_state(const int idx);
    virtual const std::string& get_tag() const override {return tag;}
    virtual void write_info(nlohmann::json& js) const override;

    static std::string tag;
    static bool registered;

private:
    Real shock_threshold;

};

#endif // SHOCKDETECTOR_H
