#ifndef PLASMA5_H
#define PLASMA5_H

#include "MFP_action.H"
#include "MFP_field.H"
#include "MFP_hydro.H"

class Plasma5 : public Action
{
  public:
    Plasma5();
    Plasma5(const int idx, const sol::table& def);
    ~Plasma5();

    enum class TimeIntegrator {
        ForwardsEuler = 0,
        BackwardsEuler,
    };

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    virtual void calc_time_derivative(MFP* mfp,
                                      Vector<UpdateData>& update,
                                      const Real time,
                                      const Real dt) override;

    void explicit_solve(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt);
    void implicit_solve(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt);

    virtual Real get_allowed_time_step(MFP* mfp) override;

    virtual std::string get_tag() const override { return tag; }
    virtual ActionType get_type() const override { return ActionType::Plasma5; }

    FieldState* field;
    Vector<HydroState*> species;

    TimeIntegrator solver;

    static bool registered;
    static std::string tag;
};

#endif  // PLASMA5_H
