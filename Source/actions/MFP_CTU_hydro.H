#ifndef CTU_H
#define CTU_H

#include "MFP_action.H"
#include "MFP_hydro.H"

class HydroCTU : public Action
{
  public:
    HydroCTU();
    HydroCTU(const int idx, const sol::table& def);
    ~HydroCTU();

    virtual std::string get_tag() const override { return tag; }
    virtual ActionType get_type() const override { return ActionType::CTU; }

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    virtual void calc_spatial_derivative(MFP* mfp,
                                         Vector<UpdateData>& dU,
                                         const Real time,
                                         const Real dt,
                                         const Real flux_register_scale) override;

    Vector<EulerianState*> states;

    bool do_CTU;

    static bool registered;
    static std::string tag;
};

#endif  // CTU_H
