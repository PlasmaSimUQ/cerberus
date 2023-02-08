#ifndef MFP_MHD_BC_H
#define MFP_MHD_BC_H
#ifdef AMREX_USE_EB

#include "sol.hpp"
#include "map"
#include "MFP_mhd_riemann.H"
#include "MFP_eulerian_bc.H"

class MHDBoundaryEB : public EulerianBoundaryEB
{
public:
    MHDBoundaryEB();
    ~MHDBoundaryEB();

    enum EBType {
        Dirichlet = 0,
        Slip,
        NoSlip
    };


    EBType type;

    MHDRiemannSolver* flux_solver = nullptr;

};

// custom boundary state where the user is able to optionally define any component
// of the primitive state vector
// Note that x points into the wall

class DirichletWallMHD : public MHDBoundaryEB
{
public:
    DirichletWallMHD();
    DirichletWallMHD(int idx, MHDRiemannSolver* flux,
                  const sol::table &bc_def);
    ~DirichletWallMHD();

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Array<Real,AMREX_SPACEDIM> wall_centre,
                       const Vector<Array4<const Real>>& all_prim,
                       const int i, const int j, const int k, const Real *dx,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F) override;


    virtual const int get_type() const override {return (int) EBType::Dirichlet;}

    Vector<std::pair<int, Real>> wall_value;
    static std::string tag;

    Vector<Real> cell_state, wall_state, normal_flux;
};


//-----------------------------------------------------------------------------

// standard slip wall

class MHDSlipWall : public MHDBoundaryEB
{
public:
    MHDSlipWall();
    MHDSlipWall(int idx, MHDRiemannSolver* flux);
    ~MHDSlipWall();

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Array<Real,AMREX_SPACEDIM> wall_centre,
                       const Vector<Array4<const Real>>& all_prim,
                       const int i, const int j, const int k, const Real *dx,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F) override;

    virtual const int get_type() const override {return (int) EBType::Slip;}

    EBType type = EBType::Slip;
    static std::string tag;
    Vector<Real> cell_state, wall_state, normal_flux;

};

#endif
#endif // MFP_MHD_BC_H
