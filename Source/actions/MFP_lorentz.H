#ifndef LORENTZ_H
#define LORENTZ_H

#include "MFP_action.H"
#include "MFP_field.H"
#include "MFP_hydro.H"

class Lorentz : public Action
{
public:
    Lorentz();
    Lorentz(const int idx, const sol::table &def);
    ~Lorentz();

    virtual void calc_time_derivative(MFP* mfp, Vector<std::pair<int,MultiFab>>& dU,  const Real time, const Real dt) override;

    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::Lorentz;}

    FieldState* field;
    Vector<HydroState*> species;

    static bool registered;
    static std::string tag;
};

#endif // LORENTZ_H
