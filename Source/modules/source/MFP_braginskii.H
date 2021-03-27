#ifndef MFP_BRAGINSKII_H
#define MFP_BRAGINSKII_H

 #include "MFP_source.H" 

class BraginskiiSource : public SourceTerm
{
public:
    BraginskiiSource();
    BraginskiiSource(const sol::table& def);
    ~BraginskiiSource();

    Vector<Real> get_slopes() const;

    Vector<Real> source(const Vector<Real> &y, const Vector<OffsetIndex> &apply) const;
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

    virtual Real get_max_freq(Vector<Real> &y) const override;
    virtual bool has_freq() const override {return true;}

    static bool valid_state(const int global_idx);
    static bool valid_solver(const int solve_idx);

    virtual SourceType get_type() const override {return SourceType::Braginskii;}
    virtual const std::string& get_tag() const override {return tag;}
    virtual std::string print() const override;

    static std::string tag;
    static bool registered;


    int linked_em, linked_hydro;
};

#endif // MFP_BRAGINSKII_H
