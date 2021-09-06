#ifndef MFP_HYDRO_BC_H
#define MFP_HYDRO_BC_H
#ifdef AMREX_USE_EB

#include "sol.hpp"
#include "map"
#include "MFP_hydro_viscous.H"
#include "MFP_hydro_riemann.H"

class HydroBoundaryEB
{
public:
    HydroBoundaryEB();
    ~HydroBoundaryEB();

    enum EBType {
        Dirichlet = 0,
        Slip,
        NoSlip
    };


    EBType type;

    virtual const EBType get_type() const = 0;

    virtual bool get_slope() const { return false;}

    HydroRiemannSolver* flux_solver = nullptr;

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Array<Real,AMREX_SPACEDIM> wall_centre,
                       Vector<Real> &cell_state,
                       Array4<const Real> const &prim4,
                       const int i, const int j, const int k, const Real *dx,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F) const = 0;


};

// custom boundary state where the user is able to optionally define any component
// of the primitive state vector
// Note that x points into the wall

class DirichletWall : public HydroBoundaryEB
{
public:
    DirichletWall();
    DirichletWall(HydroRiemannSolver* flux,
                  const sol::table &bc_def);
    ~DirichletWall();

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Array<Real,AMREX_SPACEDIM> wall_centre,
                       Vector<Real> &cell_state,
                       Array4<const Real> const &prim4,
                       const int i, const int j, const int k, const Real *dx,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F) const override;


    virtual const EBType get_type() const override {return EBType::Dirichlet;}

    Vector<std::pair<int, Real>> wall_value;

    static std::string tag;
};


//-----------------------------------------------------------------------------

// standard slip wall

class HydroSlipWall : public HydroBoundaryEB
{
public:
    HydroSlipWall();
    HydroSlipWall(HydroRiemannSolver* flux);
    ~HydroSlipWall();

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Array<Real,AMREX_SPACEDIM> wall_centre,
                       Vector<Real> &cell_state,
                       Array4<const Real> const &prim4,
                       const int i, const int j, const int k, const Real *dx,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F) const override;

    virtual const EBType get_type() const override {return type;}

    EBType type = EBType::Slip;
    static std::string tag;

};

//-----------------------------------------------------------------------------

// no slip wall using the diffusive fluxes defined by the diffusion model used
// by the hydro state that it is interacting with

class HydroNoSlipWall : public HydroBoundaryEB
{
public:
    HydroNoSlipWall();
    HydroNoSlipWall(HydroRiemannSolver* flux, HydroViscous* visc, const sol::table &bc_def);
    ~HydroNoSlipWall();

    virtual bool get_slope() const override { return true;}

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Array<Real,AMREX_SPACEDIM> wall_centre,
                       Vector<Real> &cell_state,
                       Array4<const Real> const &prim4,
                       const int i, const int j, const int k, const Real *dx,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F) const override;

    virtual const EBType get_type() const override {return EBType::NoSlip;}

    HydroViscous* viscous = nullptr;
    Array<Real,3> wall_velocity;
    Real wall_temp;

    static std::string tag;
};

#endif
#endif // MFP_HYDRO_BC_H
