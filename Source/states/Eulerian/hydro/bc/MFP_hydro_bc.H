#ifndef MFP_HYDRO_BC_H
#define MFP_HYDRO_BC_H
#ifdef AMREX_USE_EB

    #include "MFP_eulerian_bc.H"
    #include "MFP_hydro_riemann.H"
    #include "MFP_hydro_viscous.H"
    #include "map"
    #include "sol.hpp"

class HydroBoundaryEB : public EulerianBoundaryEB
{
  public:
    HydroBoundaryEB();
    ~HydroBoundaryEB();

    enum EBType { Dirichlet = 0, Slip, NoSlip, MultiState };

    EBType type;

    RiemannSolver* flux_solver = nullptr;
};

// custom boundary state where the user is able to optionally define any component
// of the primitive state vector
// Note that x points into the wall

class DirichletWall : public HydroBoundaryEB
{
  public:
    DirichletWall();
    DirichletWall(int idx, RiemannSolver* flux, const sol::table& bc_def);
    ~DirichletWall();

    virtual void solve(Array<Array<Real, 3>, 3>& wall_coord,
                       Array<Real, AMREX_SPACEDIM> wall_centre,
                       const Vector<Array4<const Real>>& all_prim,
                       const int i,
                       const int j,
                       const int k,
                       const Real* dx,
                       Array<Vector<Real>, AMREX_SPACEDIM>& F) override;

    virtual const int get_type() const override { return EBType::Dirichlet; }

    Vector<std::pair<int, Real>> wall_value;
    Vector<Real> normal_flux;
    Vector<Real> cell_state;
    Vector<Real> wall_state;

    static std::string tag;
};

//-----------------------------------------------------------------------------

// standard slip wall

class HydroSlipWall : public HydroBoundaryEB
{
  public:
    HydroSlipWall();
    HydroSlipWall(int idx, RiemannSolver* flux);
    ~HydroSlipWall();

    virtual void solve(Array<Array<Real, 3>, 3>& wall_coord,
                       Array<Real, AMREX_SPACEDIM> wall_centre,
                       const Vector<Array4<const Real>>& all_prim,
                       const int i,
                       const int j,
                       const int k,
                       const Real* dx,
                       Array<Vector<Real>, AMREX_SPACEDIM>& F) override;

    virtual const int get_type() const override { return type; }

    EBType type = EBType::Slip;
    static std::string tag;
    Vector<Real> normal_flux;
    Vector<Real> cell_state;
    Vector<Real> wall_state;
};

//-----------------------------------------------------------------------------

// no slip wall using the diffusive fluxes defined by the diffusion model used
// by the hydro state that it is interacting with

class HydroNoSlipWall : public HydroBoundaryEB
{
  public:
    HydroNoSlipWall();
    HydroNoSlipWall(int idx, RiemannSolver* flux, HydroViscous* visc, const sol::table& bc_def);
    ~HydroNoSlipWall();

    virtual bool get_slope() const override { return true; }

    virtual void solve(Array<Array<Real, 3>, 3>& wall_coord,
                       Array<Real, AMREX_SPACEDIM> wall_centre,
                       const Vector<Array4<const Real>>& all_prim,
                       const int i,
                       const int j,
                       const int k,
                       const Real* dx,
                       Array<Vector<Real>, AMREX_SPACEDIM>& F) override;

    virtual const int get_type() const override { return EBType::NoSlip; }

    HydroViscous* viscous = nullptr;
    Array<Real, 3> wall_velocity;
    Real wall_temp;
    Vector<Real> normal_flux;
    Vector<Real> cell_state;
    Vector<Real> wall_state;

    static std::string tag;
};

//-----------------------------------------------------------------------------

// a wall that needs other states info
// NOTE: When defining this BC in the input script, the state that is being
// remotely accessed must be one of the states included by the action within which
// the EB BC functions are called (i.e 'CTU')

class MultiStateWall : public HydroBoundaryEB
{
  public:
    MultiStateWall();
    MultiStateWall(const State& this_state,
                   const State& em_state,
                   RiemannSolver* flux,
                   const sol::table& bc_def);
    ~MultiStateWall();

    virtual bool get_slope() const override { return true; }

    virtual void update_data_idx(Vector<size_t>& glob2dat) override
    {
        state_idx = glob2dat[state_idx];
        em_state_idx = glob2dat[em_state_idx];
    }

    virtual void solve(Array<Array<Real, 3>, 3>& wall_coord,
                       Array<Real, AMREX_SPACEDIM> wall_centre,
                       const Vector<Array4<const Real>>& all_prim,
                       const int i,
                       const int j,
                       const int k,
                       const Real* dx,
                       Array<Vector<Real>, AMREX_SPACEDIM>& F) override;

    virtual const int get_type() const override { return EBType::MultiState; }

    static std::string tag;
    Vector<Real> normal_flux;
    Vector<Real> cell_hydro_state;
    Vector<Real> cell_em_state;
    Vector<Real> wall_state;
    int em_state_idx;
};

#endif
#endif  // MFP_HYDRO_BC_H
