#ifndef CURRENT_H
#define CURRENT_H

#include "MFP_source.H"

#include "sol.hpp"

class CurrentSource : public SourceTerm
{
    public:
        CurrentSource();
        CurrentSource(const sol::table& def);
        ~CurrentSource();

        Array<Real,3> get_current(Real x, Real y, Real z, Real t) const;

        static Vector<dual> apply_current(const Vector<dual> &y0,
                                          const Vector<OffsetIndex> &offsets,
                                          const Array<Real, 3> &j);

        virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const override;

        virtual bool has_face_src() const override {return true;}
        virtual int face_src(Real x, Real y, Real z, Real t, Vector<Real> &y0, Array<Vector<Real>, AMREX_SPACEDIM> &ydot_lo, Array<Vector<Real>, AMREX_SPACEDIM> &ydot_hi) const override;

        static bool valid_state(const int global_idx);
        static bool valid_solver(const int solve_idx);

        virtual bool has_charge_density() const override {return true;}
        virtual void calc_charge_density(const Box& box,
                                         const Real* prob_lo,
                                         const Real* dx,
                                         Real time,
                                         const Vector<FArrayBox*>& src,
                                         FArrayBox& cd,
                                         FArrayBox& J
                                         EB_OPTIONAL(,const Vector<const EBCellFlagFab*>& flag)
                                         ) const;

        virtual SourceType get_type() const override {return SourceType::Current;}
        virtual const std::string& get_tag() const override {return tag;}

        static std::string tag;
        static bool registered;

        Array<Optional3D1VFunction,3> current;
};

#endif // CURRENT_H
