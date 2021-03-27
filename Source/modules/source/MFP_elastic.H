#ifndef ELASTIC_H
#define ELASTIC_H

#include "MFP_source.H"
#include "MFP_interp.H"


class Elastic : public SourceTerm
{
    public:
        Elastic();
        Elastic(const sol::table& def);
        ~Elastic();

        static Vector<Real> elastic(const Vector<Real> &y0, const Vector<OffsetIndex> &offsets, const Interp2D &mom_interp, const Interp2D &nrg_interp);

        virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const override;
        virtual int fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const override;

        static bool valid_state(const int global_idx);
        static bool valid_solver(const int solve_idx);

        virtual SourceType get_type() const override {return SourceType::Elastic;}
        virtual const std::string& get_tag() const override {return tag;}

        static std::string tag;
        static bool registered;

        Interp2D mom_interp;
        Interp2D nrg_interp;
};


#endif
