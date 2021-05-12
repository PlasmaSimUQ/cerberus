#ifndef MFP_MHD_H
#define MFP_MHD_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include "MFP_bc.H"
#include "MFP_utility.H"
#include "MFP_state.H"
using namespace amrex;

class MhdState : public State
{
    public:
        MhdState();
        MhdState(const sol::table &def);
        virtual ~MhdState();

        enum class ConsIdx : int {
            Density=0,
            Xmom,
            Ymom,
            Zmom,
            Eden,
            Tracer,
            Bx,
            By,
            Bz,
            psi,
            NUM
        };

        enum class PrimIdx : int {
            Density=0,
            Xvel,
            Yvel,
            Zvel,
            Prs,
            Alpha,
            Bx,
            By,
            Bz,
            psi,
            NUM,
        };

        enum class FluxIdx : int {
            Density=0,
            Xvel,
            Yvel,
            Zvel,
            Prs,
            Alpha,
            Gamma,
            Mass,
            Bx,
            By,
            Bz,
            psi,
            NUM
        };

        static Vector<int> flux_vector_idx;
        static Vector<int> cons_vector_idx;
        static Vector<int> prim_vector_idx;

        const std::map<std::string, int> bc_names = {{"interior",  PhysBCType::interior},
                                                     {"inflow",    PhysBCType::inflow},
                                                     {"outflow",   PhysBCType::outflow},
                                                     {"symmetry",  PhysBCType::symmetry},
                                                     {"slipwall",  PhysBCType::slipwall},
                                                     {"noslipwall",PhysBCType::noslipwall}};

        static bool registered;
        static std::string tag;

        static Vector<std::string> cons_names;
        static Vector<std::string> prim_names;
        static Vector<set_bc> bc_set;

        bool enforce_positivity = false;
        const Vector<int> enforce_positivity_idx = {+ConsIdx::Density, +ConsIdx::Eden};

        virtual void init_from_lua() override;

        virtual bool has_field() const override {return true;}

        virtual const std::string& get_tag() const override {return tag;}
        virtual int get_type() const override {return +StateType::isMHD;}
        virtual AssociatedType get_association_type() const override {return AssociatedType::MHD;}
        virtual int n_prim() const override {return +PrimIdx::NUM;}
        virtual int n_cons() const override {return +ConsIdx::NUM;}
        virtual int n_flux() const override {return +FluxIdx::NUM;}

        virtual bool get_enforce_positivity() const override {return enforce_positivity;}
        virtual Vector<int> get_enforce_positivity_idx() const override {return enforce_positivity_idx;}

        virtual const Vector<std::string>& get_cons_names() const override;
        virtual const Vector<std::string>& get_prim_names() const override;

        virtual const Vector<std::pair<int,Real>> init_with_value() const override {
            return {{+PrimIdx::psi, 0.0}};
        }

        virtual const Vector<set_bc>& get_bc_set() const override;

        virtual Real get_mass(Real alpha) const override;
        virtual Real get_mass(const Vector<Real> &U) const override;

        virtual Real get_gamma(Real alpha) const override;
        virtual Real get_gamma(const Vector<Real> &U) const override;

        virtual Real get_cp(Real alpha) const override;
        virtual Real get_cp(const Vector<Real> &U) const override;

        virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const override;
        virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const override;
        virtual bool prim_valid(Vector<Real>& Q) const override;
        virtual bool cons_valid(Vector<Real>& U) const override;


        virtual int get_cons_density_idx() const override {return +ConsIdx::Density;}

        virtual const Vector<int> get_prim_vector_idx() const override {
            return prim_vector_idx;
        }
        virtual const Vector<int> get_flux_vector_idx() const override {
            return flux_vector_idx;
        }
        virtual const Vector<int> get_cons_vector_idx() const override {
            return cons_vector_idx;
        }

        virtual const Vector<int> get_prs_idx() const override {
            return {+PrimIdx::Prs};
        }
        virtual const Vector<int> get_nrg_idx() const override {
            return {+ConsIdx::Eden};
        }

        virtual const int get_cons_B_idx() const override {return +ConsIdx::Bx;}

        virtual const int get_cons_psi_idx() const override {
            return +ConsIdx::psi;
        }

        virtual Real get_density_from_cons(const Vector<Real> &U) const override
        {
            return U[+ConsIdx::Density];
        }
        virtual Real get_density_from_prim(const Vector<Real> &Q) const override
        {
            return Q[+PrimIdx::Density];
        }

        virtual int get_cons_tracer_idx() const override {return +ConsIdx::Tracer;}
        virtual int get_prim_alpha_idx() const override {return +PrimIdx::Alpha;}

        virtual Real get_alpha_from_cons(const Vector<Real> &U) const override
        {
            return U[+ConsIdx::Tracer]/U[+ConsIdx::Density];
        }
        virtual Real get_alpha_from_prim(const Vector<Real> &Q) const override
        {
            return Q[+PrimIdx::Alpha]/Q[+PrimIdx::Density];
        }

        virtual Real get_energy_from_cons(const Vector<Real> &U) const override;

        virtual Real get_temperature_from_cons(const Vector<Real> &U) const override;
        virtual Real get_temperature_from_prim(const Vector<Real> &Q) const override;

        virtual RealArray get_speed_from_cons(const Vector<Real> &U) const override;
        virtual RealArray get_speed_from_prim(const Vector<Real> &Q) const override;

        virtual void update_div_clean(const Real* dx) override;

        virtual void get_state_values(const Box& box,
                                      const FArrayBox& src,
                                      std::map<std::string,FArrayBox>& out,
                                      Vector<std::string>& updated
                                      EB_OPTIONAL(,const FArrayBox& vfrac)
                                      ) const override;

        virtual void calc_velocity(const Box& box,
                                   FArrayBox& src,
                                   FArrayBox& vel
                                   EB_OPTIONAL(,const EBCellFlagFab& flag)
                                   ) const override;

        virtual void load_state_for_flux(Vector<Array4<const Real> > &face,
                                                 int i, int j, int k, Vector<Real> &S) const override;

        virtual void write_info(nlohmann::json& js) const override;

};

#endif // MFP_MHD_H
