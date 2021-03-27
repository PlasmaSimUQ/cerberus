#ifndef MFP_FIELD_BC_H
#define MFP_FIELD_BC_H
#ifdef AMREX_USE_EB

#include "MFP_eb_bc.H"
#include "MFP_state.H"
#include "MFP_field.H"
#include "MFP_optional_func.H"
#include "sol.hpp"

//-----------------------------------------------------------------------------

class ConductingWall : public BoundaryEB
{
public:
    ConductingWall();
    ConductingWall(const State &istate, const sol::table &bc_def);
    ~ConductingWall();



    virtual void solve(Array<Array<Real,3>,3> &wall_coord,
                       Vector<Real> &state,
                       Vector<Real> &normal_slope,
                       Array<Vector<Real>,AMREX_SPACEDIM> &F,
                       const Real* dx) const override;

    virtual const bool has_charge() const {return true;}
    virtual const EBType& get_type() const override {return type;}

    Vector<int> vector_idx;

    bool B1_defined, B2_defined;
    Real wall_B1, wall_B2;

    bool D_defined;
    Real wall_D;

    EBType type = EBType::PerfectConductor;
    static std::string tag;

};

//-----------------------------------------------------------------------------

class ScalarPotentialWall : public BoundaryEB
{
public:
    ScalarPotentialWall();
    ScalarPotentialWall(const sol::table &bc_def);
    ~ScalarPotentialWall();

    virtual std::map<EBType, Vector<Real>> get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t) const override;

    virtual const bool has_charge() const {return true;}
    virtual const EBType& get_type() const override {return type;}
    virtual const bool time_varying() const override {return phi.has_func();}

    Optional3D1VFunction phi;

    EBType type = EBType::ScalarPotential;
    static std::string tag;

};

//-----------------------------------------------------------------------------

class SurfaceChargeWall : public BoundaryEB
{
public:
    SurfaceChargeWall();
    SurfaceChargeWall(const sol::table &bc_def);
    ~SurfaceChargeWall();

    virtual std::map<BoundaryEB::EBType, Vector<Real> > get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t) const override;

    virtual const EBType& get_type() const override {return type;}
    virtual const bool time_varying() const override {return charge_density.has_func();}

    Optional3D1VFunction charge_density;

    EBType type = EBType::SurfaceCharge;
    static std::string tag;

};

//-----------------------------------------------------------------------------

class SurfaceCurrentWall : public BoundaryEB
{
public:
    SurfaceCurrentWall();
    SurfaceCurrentWall(const sol::table &bc_def);
    ~SurfaceCurrentWall();

    virtual std::map<BoundaryEB::EBType, Vector<Real> > get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t) const override;

    virtual const EBType& get_type() const override {return type;}
    virtual const bool time_varying() const override {return current_1.has_func() || current_2.has_func();}

    Optional3D1VFunction current_1, current_2;

    EBType type = EBType::SurfaceCurrent;
    static std::string tag;

};

//-----------------------------------------------------------------------------

class VectorPotentialWall : public BoundaryEB
{
public:
    VectorPotentialWall();
    VectorPotentialWall(const sol::table &bc_def);
    ~VectorPotentialWall();

    virtual std::map<BoundaryEB::EBType, Vector<Real> > get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t) const override;

    virtual const EBType& get_type() const override {return type;}
    virtual const bool time_varying() const override {return has_function;}

    bool align_with_boundary;
    bool has_function;
    Optional3D1VFunction A_0, A_1, A_2;

    EBType type = EBType::VectorPotential;
    static std::string tag;

};

//-----------------------------------------------------------------------------

class CollectionWall : public BoundaryEB
{
public:
    CollectionWall();
    CollectionWall(const sol::table &bc_def);
    ~CollectionWall();

    virtual std::map<BoundaryEB::EBType, Vector<Real> > get_wall_state(const Array<Real,AMREX_SPACEDIM> wall_centre, const Array<Array<Real,3>,3> &wall_coord, const Real t) const override;

    virtual const EBType& get_type() const override {return type;}
    virtual const bool time_varying() const override;

    Vector<std::unique_ptr<BoundaryEB>> bcs;

    EBType type = EBType::Collection;
    static std::string tag;

};

//-----------------------------------------------------------------------------



#endif
#endif // MFP_FIELD_BC_H
