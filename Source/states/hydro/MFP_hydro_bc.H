#ifndef MFP_HYDRO_BC_H
#define MFP_HYDRO_BC_H
#ifdef AMREX_USE_EB

#include "MFP_eb_bc.H"
#include "sol.hpp"
#include "map"

//-----------------------------------------------------------------------------

// standard slip wall

class HydroSlipWall : public BoundaryEB
{
public:
    HydroSlipWall();
    HydroSlipWall(RiemannSolver* flux);
    ~HydroSlipWall();

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Vector<Real> &state,
                       Vector<Real> &normal_slope,
                       Array<Vector<Real>, AMREX_SPACEDIM> &F,
                       const Real *dx) const override;

    virtual const EBType& get_type() const override {return type;}

    EBType type = EBType::Slip;
    static std::string tag;

};

//-----------------------------------------------------------------------------

// no slip wall using the diffusive fluxes defined by the diffusion model used
// by the hydro state that it is interacting with

class HydroNoSlipWall : public BoundaryEB
{
public:
    HydroNoSlipWall();
    HydroNoSlipWall(RiemannSolver* flux, Viscous* visc, const sol::table &bc_def);
    ~HydroNoSlipWall();



    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Vector<Real> &state,
                       Vector<Real> &normal_slope,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F,
                       const Real* dx) const override;

    const bool need_slopes() const override;
    const Vector<int> get_slopes() const override;
    Vector<Real> local2global(const Array<Array<Real,3>,3> &wall_coord) const override;
    virtual const EBType& get_type() const override {return type;}

    Viscous* viscous = nullptr;
    Vector<int> slopes_idx;
    Vector<Real> wall_velocity;
    Real wall_temp;

    EBType type = EBType::NoSlip;
    static std::string tag;
};

#endif
#endif // MFP_HYDRO_BC_H
