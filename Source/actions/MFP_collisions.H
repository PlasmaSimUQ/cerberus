#ifndef COLLISIONS_H
#define COLLISIONS_H

#include "MFP_action.H"
#include "MFP_hydro.H"

class Collisions : public Action
{
public:
    Collisions();
    Collisions(const int idx, const sol::table &def);
    ~Collisions();

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    virtual void calc_time_derivative(MFP* mfp, Vector<UpdateData>& update,  const Real time, const Real dt) override;

    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::Lorentz;}

    Vector<HydroState*> species;
    Vector<Vector<Real> > ccs;

    static bool registered;
    static std::string tag;
};

#endif // LORENTZ_H
