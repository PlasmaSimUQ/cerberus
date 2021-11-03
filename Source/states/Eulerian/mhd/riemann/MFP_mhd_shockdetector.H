#ifndef MHDSHOCKDETECTOR_H
#define MHDSHOCKDETECTOR_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "sol.hpp"
#include "json.hpp"
#include "MFP_factory.H"
#include "MFP_mhd_defs.H"

using namespace amrex;

class State;


class MHDShockDetector
{
public:
    MHDShockDetector();
    ~MHDShockDetector();
    virtual Real solve(Array<Real,+MHDDef::PrimIdx::NUM> &L,
                       Array<Real,+MHDDef::PrimIdx::NUM> &R) const;

    int idx;
    State* istate;

    virtual std::string get_tag() const {return "";}
    virtual void write_info(nlohmann::json& js) const {return;}
};

template <typename D>
std::unique_ptr<MHDShockDetector> MHDShockDetectorBuilder(const sol::table& def)
{
    int global_idx = def["global_idx"];

    if (def["name"] == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

ClassFactory<MHDShockDetector>& GetMHDShockDetectorFactory();

//-----------------------------------------------------------------------------

class PressureJumpShockDetectorMHD : public MHDShockDetector
{
public:
    PressureJumpShockDetectorMHD();
    PressureJumpShockDetectorMHD(const sol::table &def);

    virtual Real solve(Array<Real,+MHDDef::PrimIdx::NUM> &L,
                       Array<Real,+MHDDef::PrimIdx::NUM> &R) const override;


    virtual std::string get_tag() const override {return tag;}
    virtual void write_info(nlohmann::json& js) const override;

    static std::string tag;
    static bool registered;

private:
    Real shock_threshold;

};

#endif // SHOCKDETECTOR_H
