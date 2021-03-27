#ifndef MFP_SOURCE_H
#define MFP_SOURCE_H

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include <functional>
#include <map>

#include "MFP_global.H"
#include "MFP_state.H"
#include "MFP_ode_system.h"
#include "MFP_reconstruction.H"

using amrex::Array;
using amrex::Vector;
using amrex::Real;

// check that the supplied index gives a state that is valid for the thing that will be
// applied to it
typedef bool (*state_valid)(const int idx);

enum class SourceType : int {
    Null = -1,
    Hydro2P,
    DampDivergence,
    DampDiv,
    Plasma5,
    Collisions,
    Braginskii,
    User,
    Ionisation,
    Recombination,
    Elastic,
    Acceleration,
    Current,
    NUM
};

struct OffsetIndex {
        OffsetIndex() : local(-1), global(-1), solver(-1), valid(true) {}
        int local;
        int global;
        int solver;
        bool valid;
};

class SourceTerm
{
    public:
        SourceTerm();
        SourceTerm(Vector<int> apply, const std::string label);
        ~SourceTerm();

        void get_includes(const sol::table& def, state_valid valid, Vector<std::string> &includes, Vector<int> &index);
        void set_reconstruction(const sol::table& def);

        ODESystem* parent;
        Vector<OffsetIndex> offsets; // {global index, solver index}
        std::string name;

        void set_local_idx_validity(const int local_idx, const bool valid);
        void set_global_idx_validity(const int global_idx, const bool valid);

        virtual int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0) const;
        virtual int fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const;
        virtual int num_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Vector<Real> &J) const;

        virtual bool has_face_src() const {return false;}
        virtual int face_src(Real x, Real y, Real z, Real t, Vector<Real> &y0, Array<Vector<Real>, AMREX_SPACEDIM> &ydot_lo, Array<Vector<Real>, AMREX_SPACEDIM> &ydot_hi) const {return 0;}

        virtual bool has_charge_density() const {return false;}
        virtual void calc_charge_density(const Box& box,
                                         const Real* prob_lo,
                                         const Real* dx,
                                         Real time,
                                         const Vector<FArrayBox*>& src,
                                         FArrayBox& cd,
                                         FArrayBox& J
                                         EB_OPTIONAL(,const Vector<const EBCellFlagFab*>& flag)
                                         ) const {return;}

        bool get_hydro_prim(Vector<Real> &y, Vector<Real> &hydro_prim, const int offset=0);

        virtual bool has_freq() const {return false;}
        virtual Real get_max_freq(Vector<Real> &y) const;

        // slopes
        virtual int num_slopes() const;
        virtual void calc_slopes(const Box& box,
                                 Vector<FArrayBox*>& src_dat,
                                 Vector<FArrayBox>& slopes,
                                 EB_OPTIONAL(Vector<const EBCellFlagFab*> &flags,)
                                 const Real *dx) const {return;}

        virtual void retrieve_slopes(
                Vector<FArrayBox>& slopes,
                const int i,
                const int j,
                const int k) {return;}

        Vector<Real> slope;

        virtual const std::string& get_name() const {return name;}
        virtual SourceType get_type() const {return SourceType::Null;}
        virtual const std::string& get_tag() const = 0;

        std::unique_ptr<Reconstruction> reconstruction;

        virtual std::string print() const {return get_tag() + " : " + get_name();}
        virtual void write_info(nlohmann::json& js) const {return;}

};



template <typename D>
std::unique_ptr<SourceTerm> SourceTermBuilder(const sol::table& def)
{
    if (def["type"] == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

PhysicsFactory<SourceTerm>& GetSourceTermFactory();

#endif // MFP_SOURCE_H
