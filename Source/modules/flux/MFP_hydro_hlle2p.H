#ifndef HLLE2P_H
#define HLLE2P_H

#include "MFP_Riemann.H"

class HLLE2P : public RiemannSolver
{
public:
    HLLE2P();
    HLLE2P(const int i);

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) override;

    virtual int get_flux_size() const override {return 7;} // keep aligned with Hydro2PState::ConsIdx::NUM

    static bool valid_state(const int idx);
    virtual const std::string& get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;
};

#endif // HLLE2P_H
