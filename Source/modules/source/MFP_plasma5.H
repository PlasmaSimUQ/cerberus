#ifndef PLASMA5_H
#define PLASMA5_H

#include "MFP_source.H"

//-----------------------------------------------------------------------------

class Plasma5 : public SourceTerm
{
    public:
        Plasma5();
        Plasma5(const sol::table& def);
        ~Plasma5();
        static Vector<Real> source(const Vector<Real> &y, const Vector<OffsetIndex> &apply);
        virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const override;
        virtual int fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const override;

        virtual bool has_charge_density() const override {return true;}
        virtual void calc_charge_density(const Box& box,
                                         const Real* prob_lo,
                                         const Real* dx,
                                         Real time,
                                         const Vector<FArrayBox*>& src,
                                         FArrayBox& cd,
                                         FArrayBox& J
                                         EB_OPTIONAL(,const Vector<const EBCellFlagFab*>& flag)
                                         ) const override;

        virtual bool has_face_src() const override {return true;}
        virtual int face_src(Real x, Real y, Real z, Real t, Vector<Real> &y0, Array<Vector<Real>, AMREX_SPACEDIM> &ydot_lo, Array<Vector<Real>, AMREX_SPACEDIM> &ydot_hi) const override;

        virtual Real get_max_freq(Vector<Real> &y) const override;
        virtual bool has_freq() const override {return true;}

        static bool valid_state(const int global_idx);
        static bool valid_solver(const int solve_idx);

        virtual SourceType get_type() const override {return SourceType::Plasma5;}
        virtual const std::string& get_tag() const override {return tag;}
        virtual std::string print() const override;

        static std::string tag;
        static bool registered;

        Vector<int> linked_electron_idx, linked_ion_idx;
        int linked_EM_idx;

};

#endif // PLASMA5_H
