#ifndef HYDROAUSMDV_H
#define HYDROAUSMDV_H

#include "MFP_Riemann.H"

class HydroAUSMDV : public RiemannSolver
{
public:
    HydroAUSMDV();
    HydroAUSMDV(const int i);

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) const override;

    virtual int get_flux_size() const override {return 6;} // keep aligned with HydroState::ConsIdx::NUM

    static bool valid_state(const int idx);
    virtual const std::string& get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;
};

#endif // HYDROAUSMDV_H
