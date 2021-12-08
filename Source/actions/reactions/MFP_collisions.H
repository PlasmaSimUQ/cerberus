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

    struct SpeciesInfo {
        std::string name;
        int idx;
        int state_idx;
        int alpha_idx;
        Real rho, m, q, T, q2, n;
        Array<Real,3> vel;
        Array<Real,+HydroDef::ConsIdx::NUM> delta;
    };

    void init_species_info(const sol::table &def);

    virtual void calc_update(){}

    Vector<HydroState*> states;
    Vector<SpeciesInfo> species_info;

};

#endif // LORENTZ_H
