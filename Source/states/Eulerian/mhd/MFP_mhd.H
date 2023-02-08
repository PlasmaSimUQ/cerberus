#ifndef MHDSTATE_H
#define MHDSTATE_H
#include "forward.hpp"

class MHDState;

#include "MFP_eulerian.H"
#include "MFP_bc.H"
#include "MFP.H"
#include "MFP_mhd_defs.H"
#include "MFP_reconstruction.H"
#include "MFP_mhd_riemann.H"
#include "MFP_mhd_bc.H"
#include "MFP_mhd_shockdetector.H"


class MHDState : public EulerianState
{
public:
    MHDState();
    MHDState(const sol::table& def);
    ~MHDState();

    static bool registered;
    static std::string tag;

    static Vector<std::string> cons_names;
    static Vector<std::string> prim_names;
    static Array<int,2> cons_vector_idx;
    static Array<int,2> prim_vector_idx;
    static std::map<std::string, int> bc_names;
    static Vector<set_bc> bc_set;

    Real gamma;

    Real div_transport;

    virtual StateType get_type() const override {return StateType::MHD;}

    virtual size_t n_cons() const override {return +MHDDef::ConsIdx::NUM;}
    virtual size_t n_prim() const override {return +MHDDef::PrimIdx::NUM;}

    virtual const Vector<int> get_cons_vector_idx() const override {return arr2vec(cons_vector_idx);}
    virtual const Vector<int> get_prim_vector_idx() const override {return arr2vec(prim_vector_idx);}

#ifdef AMREX_USE_EB
    virtual void set_eb_bc(const sol::table &bc_def) override;
#endif

    virtual int get_num_grow() const override {return num_grow;}
    virtual void set_num_grow(int n) override {num_grow = std::max(n, num_grow);}

    virtual const Vector<std::string>& get_cons_names() const override {return cons_names;}
    virtual const Vector<std::string>& get_prim_names() const override {return prim_names;}

    void set_udf();

    void set_flux();
    void set_shock_detector();
    void set_refinement();

    virtual void variable_setup(Vector<int> periodic) override;
    virtual void init_from_lua() override;

    virtual bool cons2prim(Vector<Real>& U, Vector<Real> &Q) const override;
    virtual void prim2cons(Vector<Real> &Q, Vector<Real> &U) const override;

    bool prim_valid(const Vector<Real>& Q) const;
    bool cons_valid(const Vector<Real>& U) const;

    RealArray get_speed_from_cons(const Vector<Real> &U) const;
    RealArray get_speed_from_prim(const Vector<Real>& Q) const;

    void calc_velocity(const Box& box,
                       FArrayBox& cons,
                       FArrayBox &prim
                   #ifdef AMREX_USE_EB
                       ,const FArrayBox& vfrac
                   #endif
                       ) const;


    virtual void write_info(nlohmann::json& js) const override;

    static MHDState& get_state_global(const int idx) {return static_cast<MHDState&>(MFP::get_state(idx));}
    static MHDState& get_state(const std::string& name) {return static_cast<MHDState&>(MFP::get_state(name));}

};

#endif // MHDSTATE_H
