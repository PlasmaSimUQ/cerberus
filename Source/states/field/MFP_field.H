#ifndef MFP_FIELD_H
#define MFP_FIELD_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include "MFP_bc.H"
#include "MFP_utility.H"
#include "MFP_state.H"
using namespace amrex;

struct FieldState : public State
{
    public:
        FieldState();
        FieldState(const sol::table &def);
        ~FieldState();

        enum class ConsIdx : int {
            Dx=0,
            Dy,
            Dz,
            Bx,
            By,
            Bz,
            phi,
            psi,
            mu,
            ep,
            NUM
        };

        enum class PrimIdx : int {
            Dx=0,
            Dy,
            Dz,
            Bx,
            By,
            Bz,
            phi,
            psi,
            mu,
            ep,
            NUM
        };

        enum class FluxIdx : int {
            Dx=0,
            Dy,
            Dz,
            Bx,
            By,
            Bz,
            phi,
            psi,
            mu,
            ep,
            NUM
        };

        static Vector<int> flux_vector_idx;
        static Vector<int> cons_vector_idx;
        static Vector<int> prim_vector_idx;

        const std::map<std::string, int> bc_names = {{"interior",  PhysBCType::interior},
                                                     {"inflow",    PhysBCType::inflow},
                                                     {"outflow",   PhysBCType::outflow},
                                                     {"symmetry",  PhysBCType::symmetry},
                                                     {"asymmetry",  4}};

        Real fastest_speed = 0.0;
        bool is_static = false;

        static bool registered;
        static std::string tag;

        static Vector<std::string> cons_names;
        static Vector<std::string> input_names;
        static Vector<set_bc> bc_set;

        virtual void init_from_lua() override;

        virtual bool is_transported() const override {return !is_static;}
        virtual bool has_field() const override {return true;}

#ifdef AMREX_USE_EB
        virtual void set_eb_bc(const sol::table &bc_def) override;
#endif

        virtual const std::string& get_tag() const override {return tag;}
        virtual int get_type() const override {return +StateType::isField;}
        virtual AssociatedType get_association_type() const override {return AssociatedType::Field;}
        virtual int n_prim() const override {return +PrimIdx::NUM;}
        virtual int n_cons() const override {return +ConsIdx::NUM;}
        virtual int n_flux() const override {return +FluxIdx::NUM;}
        virtual bool prim_valid(Vector<Real>& Q) const override {return true;}
        virtual bool cons_valid(Vector<Real>& U) const override {return true;}

        virtual const Vector<std::string>& get_cons_names() const override;
        virtual const Vector<std::string>& get_prim_names() const override;

        virtual const Vector<std::pair<int,Real>> init_with_value() const override {
            return {
                {+PrimIdx::phi, 0.0},
                {+PrimIdx::psi, 0.0},
                {+PrimIdx::mu, 1.0},
                {+PrimIdx::ep, 1.0},
            };
        }

        virtual const Vector<set_bc>& get_bc_set() const override;

        virtual bool cons2prim(Vector<Real>& U) const override;
        virtual void prim2cons(Vector<Real>& U) const override;

        virtual int get_cons_density_idx() const override {return -1;}

        virtual const Vector<int> get_prim_vector_idx() const override {
            return prim_vector_idx;
        }

        virtual const Vector<int> get_flux_vector_idx() const override {
            return flux_vector_idx;
        }

        virtual const Vector<int> get_cons_vector_idx() const override {
            return cons_vector_idx;
        }

        virtual const int get_cons_D_idx() const override {return +ConsIdx::Dx;}
        virtual const int get_cons_B_idx() const override {return +ConsIdx::Bx;}

        virtual const int get_cons_phi_idx() const override {return +ConsIdx::phi;}
        virtual const int get_cons_psi_idx() const override {return +ConsIdx::psi;}

        virtual Real get_energy_from_cons(const Vector<Real> &U) const override;

        virtual RealArray get_speed_from_cons(const Vector<Real> &U) const override;
        virtual RealArray get_speed_from_prim(const Vector<Real> &Q) const override;

        virtual void get_state_values(const Box& box,
                                      const FArrayBox& src,
                                      std::map<std::string,FArrayBox>& out,
                                      Vector<std::string>& updated
                                      EB_OPTIONAL(,const FArrayBox& vfrac)
                                      ) const override;

        void calc_primitives(const Box& box,
                             FArrayBox& cons,
                             FArrayBox& prim,
                             const Real* dx,
                             const Real t,
                             const Real* prob_lo
                             EB_OPTIONAL(,const FArrayBox& vfrac)
                             ) const override;

        void calc_reconstruction(const Box& box,
                                 FArrayBox &prim,
                                 Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                 Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                 EB_OPTIONAL(,const EBCellFlagFab &flag)
                                 EB_OPTIONAL(,const FArrayBox &vfrac)
                                 ) const override;

        virtual void load_state_for_flux(Vector<Array4<const Real> > &face,
                                                 int i, int j, int k, Vector<Real> &S) const override;


        virtual void write_info(nlohmann::json& js) const override;
};

#endif // MFP_FIELD_H
