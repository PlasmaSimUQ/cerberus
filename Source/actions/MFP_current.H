#ifndef CURRENT_H
#define CURRENT_H

#include "MFP_action.H"
#include "MFP_field.H"

class Current : public Action
{
  public:
    Current();
    Current(const int idx, const sol::table& def);
    ~Current();

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    virtual void calc_time_derivative(MFP* mfp,
                                      Vector<UpdateData>& update,
                                      const Real time,
                                      const Real dt) override;

    virtual std::string get_tag() const override { return tag; }
    virtual ActionType get_type() const override { return ActionType::Current; }

    FieldState* field;
    Array<Optional3D1VFunction, 3> current;

    static bool registered;
    static std::string tag;

    Real scale_factor;
};

#endif  // CURRENT_H
