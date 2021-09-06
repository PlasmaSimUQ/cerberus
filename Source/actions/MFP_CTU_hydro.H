#ifndef CTU_H
#define CTU_H

#include "MFP_action.H"
#include "MFP_hydro.H"

class HydroCTU : public Action
{
public:
    HydroCTU();
    HydroCTU(const int idx, const sol::table &def);
    ~HydroCTU();

    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::CTU;}

    virtual void calc_spatial_derivative(MFP* mfp, Vector<std::pair<int,MultiFab>>& dU, const Real time, const Real dt) override;

    Vector<EulerianState*> states;

    bool do_CTU;


    static bool registered;
    static std::string tag;
};

#endif // CTU_H
