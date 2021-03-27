#ifndef MHDHLLE_H
#define MHDHLLE_H

#include "MFP_Riemann.H"

class MhdHLLE : public RiemannSolver
{
public:

    MhdHLLE();
    MhdHLLE(const int i);

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) const override;

    virtual int get_flux_size() const override {return 10;} // keep aligned with MhdState::ConsIdx::NUM

    static bool valid_state(const int idx);
    virtual const std::string& get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;
};

#endif // MHDHLLE_H
