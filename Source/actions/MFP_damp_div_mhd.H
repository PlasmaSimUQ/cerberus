#ifndef DAMPDIVERGENCEMHD_H
#define DAMPDIVERGENCEMHD_H

#include "MFP_action.H"
#include "MFP_mhd.H"

class DampDivergenceMHD : public Action
{
  public:
    DampDivergenceMHD();
    ~DampDivergenceMHD();
    DampDivergenceMHD(const int idx, const sol::table& def);

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    virtual void calc_time_derivative(MFP* mfp,
                                      Vector<UpdateData>& update,
                                      const Real time,
                                      const Real dt) override;

    virtual std::string get_tag() const override { return tag; }
    virtual ActionType get_type() const override { return ActionType::DampDivergenceMHD; }

    static bool registered;
    static std::string tag;

    MHDState* mhd;

    Real div_damp, div_transport;
};

#endif  // DAMPDIVERGENCEMHD_H
