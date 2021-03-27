#ifndef BOUNDARYEB_H
#define BOUNDARYEB_H
#ifdef AMREX_USE_EB

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

#include "sol.hpp"
#include "MFP_Riemann.H"
#include "MFP_viscous.H"


using namespace amrex;


class BoundaryEB
{
public:
    BoundaryEB();
    ~BoundaryEB();

    enum EBType {
        Dirichlet = 0,
        Slip,
        NoSlip,
        PerfectConductor,
        ScalarPotential,
        VectorPotential,
        SurfaceCharge,
        SurfaceCurrent,
        Collection,
    };

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Vector<Real> &cell_state,
                       Vector<Real> &normal_slope,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F,
                       const Real* dx) const {return;}

    virtual std::map<EBType, Vector<Real>> get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t ) const {return std::map<EBType, Vector<Real>> {};}

    virtual const bool need_slopes() const;
    virtual const Vector<int> get_slopes() const;
    virtual const bool has_charge() const {return false;}
    virtual Vector<Real> local2global(const Array<Array<Real,3>,3> &wall_coord) const;
    virtual const EBType& get_type() const {return type;}
    virtual const bool time_varying() const {return false;}

    RiemannSolver* flux_solver = nullptr;

    EBType type;
};


//-----------------------------------------------------------------------------

// custom boundary state where the user is able to optionally define any component
// of the primitive state vector
// Note that x points into the wall

class DirichletWall : public BoundaryEB
{
public:
    DirichletWall();
    DirichletWall(RiemannSolver* flux,
                  const Vector<std::string> &names,
                  const Vector<int> &vec_idx,
                  const sol::table &bc_def);
    ~DirichletWall();

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Vector<Real> &state,
                       Vector<Real> &normal_slope,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F,
                       const Real* dx) const override;

    virtual const EBType& get_type() const override {return type;}

    Vector<std::pair<int, Real>> wall_value;
    Vector<int> vector_idx;
    EBType type = EBType::Dirichlet;

    static std::string tag;
};

#endif
#endif // BOUNDARYEB_H
