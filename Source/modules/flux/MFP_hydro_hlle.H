#ifndef HYDROHLLE_H
#define HYDROHLLE_H

#include "MFP_Riemann.H"
#include "MFP_hydro.H"

class HydroHLLE : public RiemannSolver
{
public:
    HydroHLLE();
    HydroHLLE(const int i);

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) override;

    virtual int get_flux_size() const override {return +HydroState::ConsIdx::NUM;}

    static bool valid_state(const int idx);
    virtual const std::string& get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;

};

#endif // HYDROHLLE_H
