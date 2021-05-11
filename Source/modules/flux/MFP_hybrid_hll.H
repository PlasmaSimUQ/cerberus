#ifndef HYDROHYBRIDHLL_H
#define HYDROHYBRIDHLL_H

#include "MFP_Riemann.H"
#include "MFP_hydro_hlle.H"
#include "MFP_hllc.H"
#include "MFP_hydro.H"


class HydroHybridHLL : public RiemannSolver
{
public:
    HydroHybridHLL();
    HydroHybridHLL(const int i);

    HydroHLLC hllc;
    HydroHLLE hlle;

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) override;

    virtual int get_flux_size() const override {return 6;} // keep aligned with HydroState::ConsIdx::NUM


    static bool valid_state(const int idx);
    virtual const std::string& get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

private:
    Vector<Real> F_hlle, F_hllc;
};

#endif // HYDROHYBRIDHLL_H
