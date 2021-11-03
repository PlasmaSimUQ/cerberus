#ifndef MFP_MHD_BC_H
#define MFP_MHD_BC_H
#ifdef AMREX_USE_EB

#include "sol.hpp"
#include "map"
#include "MFP_mhd_riemann.H"

class MHDBoundaryEB
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

    virtual const EBType get_type() const = 0;

    virtual bool get_slope() const { return false;}

    MHDRiemannSolver* flux_solver = nullptr;

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Array<Real,AMREX_SPACEDIM> wall_centre,
                       Array<Real,+MHDDef::PrimIdx::NUM> &cell_state,
                       Array4<const Real> const &prim4,
                       const int i, const int j, const int k, const Real *dx,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F) const = 0;


};

// custom boundary state where the user is able to optionally define any component
// of the primitive state vector
// Note that x points into the wall

class DirichletWallMHD : public MHDBoundaryEB
{
public:
    DirichletWallMHD();
    DirichletWallMHD(MHDRiemannSolver* flux,
                  const sol::table &bc_def);
    ~DirichletWallMHD();

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Array<Real,AMREX_SPACEDIM> wall_centre,
                       Array<Real,+MHDDef::PrimIdx::NUM> &cell_state,
                       Array4<const Real> const &prim4,
                       const int i, const int j, const int k, const Real *dx,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F) const override;


    virtual const EBType get_type() const override {return EBType::Dirichlet;}

    Vector<std::pair<int, Real>> wall_value;

    static std::string tag;
};


//-----------------------------------------------------------------------------

// standard slip wall

class MHDSlipWall : public MHDBoundaryEB
{
public:
    MHDSlipWall();
    MHDSlipWall(MHDRiemannSolver* flux);
    ~MHDSlipWall();

    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Array<Real,AMREX_SPACEDIM> wall_centre,
                       Array<Real,+MHDDef::PrimIdx::NUM> &cell_state,
                       Array4<const Real> const &prim4,
                       const int i, const int j, const int k, const Real *dx,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F) const override;

    virtual const EBType get_type() const override {return type;}

    EBType type = EBType::Slip;
    static std::string tag;

};

#endif
#endif // MFP_MHD_BC_H
