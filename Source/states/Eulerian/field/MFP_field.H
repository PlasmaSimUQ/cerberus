#ifndef FIELDSTATE_H
#define FIELDSTATE_H

#include "MFP.H"
#include "MFP_bc.H"
#include "MFP_eulerian.H"
#include "MFP_field_bc.H"
#include "MFP_field_defs.H"
#include "MFP_field_riemann.H"
#include "MFP_state.H"

class FieldState : public EulerianState
{
  public:
    FieldState();
    FieldState(const sol::table& def);
    ~FieldState();

    static bool registered;
    static std::string tag;

    static Vector<std::string> cons_names;
    static Array<int, +FieldDef::VectorIdx::Cons> vector_idx;

    static std::map<std::string, int> bc_names;
    static Vector<set_bc> bc_set;

    int is_static;
    Real fastest_speed;

    Real div_speed, div_damping;
    Real relative_div_speed;

    virtual StateType get_type() const override { return StateType::Field; }
    virtual size_t n_cons() const override { return +FieldDef::ConsIdx::NUM; }
    virtual size_t n_prim() const override { return +FieldDef::ConsIdx::NUM; }
    virtual int get_num_grow() const override { return num_grow; }
    virtual void set_num_grow(int n) override { num_grow = std::max(n, num_grow); }
    virtual bool is_transported() const override { return !is_static; }

    virtual const Vector<int> get_cons_vector_idx() const override { return arr2vec(vector_idx); }
    virtual const Vector<int> get_prim_vector_idx() const override { return arr2vec(vector_idx); }

    virtual const Vector<std::string>& get_cons_names() const override { return cons_names; }
    virtual const Vector<std::string>& get_prim_names() const override { return cons_names; }
    virtual Vector<std::string> get_plot_output_names() const override;

    virtual void get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string, FArrayBox>& out,
                                 Vector<std::string>& updated
#ifdef AMREX_USE_EB
                                 ,
                                 const FArrayBox& vfrac
#endif
    ) const override;

    void set_udf();
    void set_flux();
    void set_refinement();

    virtual void variable_setup(Vector<int> periodic) override;
    virtual void init_from_lua() override;

#ifdef AMREX_USE_EB
    virtual void set_eb_bc(const sol::table& bc_def) override;
#endif

    virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const override;
    virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const override;

    virtual RealArray get_speed_from_cons(const Vector<Real>& U) const
    {
        return {AMREX_D_DECL(fastest_speed, fastest_speed, fastest_speed)};
    }
    virtual RealArray get_speed_from_prim(const Vector<Real>& Q) const
    {
        return {AMREX_D_DECL(fastest_speed, fastest_speed, fastest_speed)};
    }

    virtual Real get_allowed_time_step(MFP* mfp) const override;

#ifdef AMREX_USE_EB

    virtual void get_wall_value(const Box& box,
                                Vector<FArrayBox*> bcs_data,
                                const EBCellFlagFab& flag,
                                const CutFab& bc_idx,
                                const FArrayBox& bcent,
                                const FArrayBox& bnorm,
                                const Real t,
                                const Real* dx,
                                const Real* prob_lo) const override;
#endif

    void write_info(nlohmann::json& js) const override;

    static FieldState& get_state_global(const int idx)
    {
        return static_cast<FieldState&>(MFP::get_state(idx));
    }
    static FieldState& get_state(const std::string& name)
    {
        return static_cast<FieldState&>(MFP::get_state(name));
    }
};

#endif  // FIELDSTATE_H
