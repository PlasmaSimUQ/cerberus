#ifndef DAMPDIVERGENCECORRECTION_H
#define DAMPDIVERGENCECORRECTION_H

#include "MFP_source.H"

//-----------------------------------------------------------------------------

class DampDivergenceCorrection : public SourceTerm
{
public:
    DampDivergenceCorrection();
    DampDivergenceCorrection(const sol::table& def);
    ~DampDivergenceCorrection();
    virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const override;
    virtual int fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const override;

    static bool valid_state(const int global_idx);
    static bool valid_solver(const int solve_idx);

    virtual SourceType get_type() const override {return SourceType::DampDivergence;}
    virtual const std::string& get_tag() const override {return tag;}

    static std::string tag;
    static bool registered;

    Real damping; // this is set at init

    virtual void write_info(nlohmann::json& js) const override;

};

#endif // DAMPDIVERGENCECORRECTION_H
