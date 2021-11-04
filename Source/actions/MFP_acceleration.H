#ifndef ACCELERATION_H
#define ACCELERATION_H

#include "MFP_action.H"
#include "MFP_hydro.H"

class Acceleration : public Action
{
public:
    Acceleration();
    Acceleration(const int idx, const sol::table &def);
    ~Acceleration();

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    virtual void calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt) override;

    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::Acceleration;}

    Vector<HydroState*> species;

    Array<Real,3> acc;

    static bool registered;
    static std::string tag;
};

#endif // ACCELERATION_H
