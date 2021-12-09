#ifndef ACTION_H
#define ACTION_H

#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include "sol.hpp"
#include "MFP_factory.H"
#include "MFP.H"

class Action;
class MFP;

using namespace amrex;

struct UpdateData {

    enum class Status {
        Inactive=0,
        Local,
        Filled,
        Expired,
        Changed,
        Zero,
    };

    Status dU_status = Status::Inactive;
    Status U_status = Status::Inactive;
    MultiFab U, dU;
};

class Action
{
public:
    Action();
    ~Action();

    enum class ActionType {
        Acceleration = 0,
        CTU,
        Flux,
        Elliptic,
        HydroTracer,
        Lorentz,
        Plasma5,
#ifdef SYMPLECTIC
        Symplectic,
#endif
        Current,
        UserDefined,
        DampDivergenceMHD,
        Collisions,
        GasKinetics
    };

    int action_idx;
    std::string name;
    Vector<int> state_indexes;

    // copy state data
    void get_data(MFP* mfp, Vector<Array<int,2>>& options, Vector<UpdateData>& update, const Real time) const;
    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const {}

    // calculate time and space derivatives for situations where we can accumulate the effects
    virtual void calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt) {}
    virtual void calc_spatial_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt, const Real flux_register_scale) {}

    // directly update quantities where an accumulation is not possible, these will overwrite stuff!!
    virtual void apply_time_derivative(MFP* mfp, const Real time, const Real dt) {}
    virtual void apply_spatial_derivative(MFP* mfp, const Real time, const Real dt) {}

    // just do something, typically at the end of the update main loop
    virtual void apply_change(MFP* mfp, const Real time, const Real dt) {}

    virtual Real get_allowed_time_step(MFP* mfp) {return std::numeric_limits<Real>::max();}

    virtual std::string get_tag() const = 0;
    virtual ActionType get_type() const = 0;
};

template <typename D>
std::unique_ptr<Action> ActionBuilder(const sol::table& def)
{
    int idx = def["src_idx"];

    if (def["type"] == D::tag) {
        return std::unique_ptr<D>(new D(idx, def));
    } else {
        return nullptr;
    }
}

ClassFactory<Action> &GetActionFactory();

#endif // ACTION_H
