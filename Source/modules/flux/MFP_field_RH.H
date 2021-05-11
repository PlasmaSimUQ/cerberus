#ifndef FIELDRH_H
#define FIELDRH_H

#include "MFP_Riemann.H"

class FieldRH : public RiemannSolver
{
public:

    FieldRH();
    FieldRH( const int i);

    virtual void solve(Vector<Real> &L,
                       Vector<Real> &R,
                       Vector<Real> &F,
                       Real* shk) override;

    virtual int get_flux_size() const override {return 10;} // keep aligned with FieldState::ConsIdx::NUM

    static bool valid_state(const int idx);
    virtual const std::string& get_tag() const override {return tag;}
    static std::string tag;
    static bool registered;
};

#endif // FIELDRH_H
