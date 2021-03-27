#ifndef HYDROTWOPRESSURE_H
#define HYDROTWOPRESSURE_H

#include "MFP_source.H"

//-----------------------------------------------------------------------------

class HydroTwoPressure : public SourceTerm
{
    public:
        HydroTwoPressure();
        HydroTwoPressure(const sol::table& def);
        ~HydroTwoPressure();
        static Vector<Real> source(const Vector<Real> &y, const Vector<Real> &slopes, const Vector<OffsetIndex> &apply);
        virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const override;
        virtual int fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const override;

        virtual int num_slopes() const override;
        virtual void calc_slopes(const Box& box,
                                 Vector<FArrayBox*>& src_dat,
                                 Vector<FArrayBox>& slopes,
                                 EB_OPTIONAL(Vector<const EBCellFlagFab*> &flags,)
                                 const Real *dx) const override;

        virtual void retrieve_slopes(
                Vector<FArrayBox>& slopes,
                const int i,
                const int j,
                const int k) override;

        static bool valid_state(const int global_idx);
        static bool valid_solver(const int solve_idx);

        int linked_em;
        Vector<int> linked_hydro;

        virtual SourceType get_type() const override {return SourceType::Hydro2P;}
        virtual const std::string& get_tag() const override {return tag;}

        static std::string tag;
        static bool registered;
};
#endif // HYDROTWOPRESSURE_H
