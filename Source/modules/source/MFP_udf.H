#ifndef USERDEFINEDSOURCE_H
#define USERDEFINEDSOURCE_H

#include "MFP_source.H"

//-----------------------------------------------------------------------------

class UserDefinedSource : public SourceTerm
{
public:
    UserDefinedSource();
    UserDefinedSource(const sol::table& def);
    ~UserDefinedSource();
    virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const override;

    static bool valid_state(const int global_idx);
    static bool valid_solver(const int solve_idx);

    virtual SourceType get_type() const override {return SourceType::User;}
    virtual const std::string& get_tag() const override {return tag;}
    virtual std::string print() const override;

    static std::string tag;
    static bool registered;

    Vector<Vector<Optional3D1VFunction>> terms; // [local_id][variable_id]

};

#endif // USERDEFINEDSOURCE_H
