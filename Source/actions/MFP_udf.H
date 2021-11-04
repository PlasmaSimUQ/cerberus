#ifndef USERDEFINED_H
#define USERDEFINED_H

#include "MFP_action.H"
#include "MFP_optional_func.H"
#include "MFP_eulerian.H"


class UserDefined : public Action
{
public:
    UserDefined();
    UserDefined(const int idx, const sol::table &def);
    ~UserDefined();

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    virtual void calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt) override;

    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::UserDefined;}

    static bool registered;
    static std::string tag;

    EulerianState* state;

    Vector<Optional3D1VFunction> terms;
};

#endif // USERDEFINED_H
