#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include "MFP_action.H"
#include "MFP_hydro.H"
#include "MFP_field.H"
#include "MFP_mhd.H"

class Elliptic : public Action
{
public:
    Elliptic();
    Elliptic(const int idx, const sol::table &def);
    ~Elliptic();

    virtual void apply_change(MFP* mfp, const Real time, const Real dt) override;

    void solve_static_fields(MFP* mfp, const Real time);
    void project_divergence(MFP* mfp, const Real time);
    void project_divergence_field(MFP* mfp, const Real time);
    void project_divergence_mhd(MFP* mfp, const Real time);
    void solve_divergence(MFP* mfp,
                          const int state_idx,
                          const int vector_idx,
                          const int phi_idx,
                          const BCRec &bc,
                          Vector<MultiFab*> S_cc_ptr={});


    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::Elliptic;}

    State::StateType select;
    FieldState* field;
    MHDState* mhd;
    Vector<HydroState*> species;

    // projection method?
    int projection;

    // relative tolerance
    Real reltol;

    // absolute tolerance
    Real abstol;

    static bool registered;
    static std::string tag;
};

#endif // ELLIPTIC_H
