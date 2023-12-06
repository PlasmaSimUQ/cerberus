#ifndef MFP_FIELD_BC_H
#define MFP_FIELD_BC_H
#ifdef AMREX_USE_EB

    #include "MFP_eulerian_bc.H"
    #include "MFP_field_riemann.H"
    #include "MFP_optional_func.H"
    #include "sol.hpp"

class FieldBoundaryEB : public EulerianBoundaryEB
{
  public:
    FieldBoundaryEB();
    ~FieldBoundaryEB();

    enum class EBType {
        Dirichlet = 0,
        PerfectConductor,
        ScalarPotential,
        VectorPotential,
        SurfaceCharge,
        SurfaceCurrent,
        Collection,
        NUM
    };

    virtual const bool has_charge() const { return false; }
    virtual const bool time_varying() const { return false; }

    RiemannSolver* flux_solver = nullptr;

    EBType type;
};

//-----------------------------------------------------------------------------

class DefinedWall : public FieldBoundaryEB
{
  public:
    DefinedWall();
    DefinedWall(int idx, RiemannSolver* flux, const sol::table& bc_def);
    ~DefinedWall();

    void solve(Array<Array<Real, 3>, 3>& wall_coord,
               Array<Real, AMREX_SPACEDIM> wall_centre,
               const Vector<Array4<const Real>>& all_prim,
               const int i,
               const int j,
               const int k,
               const Real* dx,
               Array<Vector<Real>, AMREX_SPACEDIM>& F) override;

    virtual const int get_type() const override { return (int)EBType::Dirichlet; }

    Vector<std::pair<int, Real>> wall_value;

    EBType type = EBType::Dirichlet;
    static std::string tag;

    Vector<Real> cell_state, wall_state, normal_flux;
};

//-----------------------------------------------------------------------------

class ConductingWall : public FieldBoundaryEB
{
  public:
    ConductingWall();
    ConductingWall(int idx, RiemannSolver* flux, const sol::table& bc_def);
    ~ConductingWall();

    void solve(Array<Array<Real, 3>, 3>& wall_coord,
               Array<Real, AMREX_SPACEDIM> wall_centre,
               const Vector<Array4<const Real>>& all_prim,
               const int i,
               const int j,
               const int k,
               const Real* dx,
               Array<Vector<Real>, AMREX_SPACEDIM>& F) override;

    virtual const bool has_charge() const override { return true; }
    virtual const int get_type() const override { return (int)EBType::PerfectConductor; }

    bool B1_defined, B2_defined;
    Real wall_B1, wall_B2;

    bool D_defined;
    Real wall_D;

    EBType type = EBType::PerfectConductor;
    static std::string tag;

    Vector<Real> cell_state, wall_state, normal_flux;
};

//-----------------------------------------------------------------------------

class ScalarPotentialWall : public FieldBoundaryEB
{
  public:
    ScalarPotentialWall();
    ScalarPotentialWall(int idx, const sol::table& bc_def);
    ~ScalarPotentialWall();

    virtual Vector<Vector<Real>> get_wall_state(const Array<Real, AMREX_SPACEDIM> wall_centre,
                                                const Array<Array<Real, 3>, 3>& wall_coord,
                                                const Real t) const override;

    virtual const bool has_charge() const override { return true; }
    virtual const int get_type() const override { return (int)EBType::ScalarPotential; }
    virtual const bool time_varying() const override { return phi.has_func(); }

    Optional3D1VFunction phi;

    EBType type = EBType::ScalarPotential;
    static std::string tag;
};

//-----------------------------------------------------------------------------

class SurfaceChargeWall : public FieldBoundaryEB
{
  public:
    SurfaceChargeWall();
    SurfaceChargeWall(int idx, const sol::table& bc_def);
    ~SurfaceChargeWall();

    virtual Vector<Vector<Real>> get_wall_state(const Array<Real, AMREX_SPACEDIM> wall_centre,
                                                const Array<Array<Real, 3>, 3>& wall_coord,
                                                const Real t) const override;

    virtual const int get_type() const override { return (int)EBType::SurfaceCharge; }
    virtual const bool time_varying() const override { return charge_density.has_func(); }

    Optional3D1VFunction charge_density;

    EBType type = EBType::SurfaceCharge;
    static std::string tag;
};

//-----------------------------------------------------------------------------

class SurfaceCurrentWall : public FieldBoundaryEB
{
  public:
    SurfaceCurrentWall();
    SurfaceCurrentWall(int idx, const sol::table& bc_def);
    ~SurfaceCurrentWall();

    virtual Vector<Vector<Real>> get_wall_state(const Array<Real, AMREX_SPACEDIM> wall_centre,
                                                const Array<Array<Real, 3>, 3>& wall_coord,
                                                const Real t) const override;

    virtual const int get_type() const override { return (int)EBType::SurfaceCurrent; }
    virtual const bool time_varying() const override
    {
        return current_1.has_func() || current_2.has_func();
    }

    Optional3D1VFunction current_1, current_2;

    EBType type = EBType::SurfaceCurrent;
    static std::string tag;
};

//-----------------------------------------------------------------------------

class VectorPotentialWall : public FieldBoundaryEB
{
  public:
    VectorPotentialWall();
    VectorPotentialWall(int idx, const sol::table& bc_def);
    ~VectorPotentialWall();

    virtual Vector<Vector<Real>> get_wall_state(const Array<Real, AMREX_SPACEDIM> wall_centre,
                                                const Array<Array<Real, 3>, 3>& wall_coord,
                                                const Real t) const override;

    virtual const int get_type() const override { return (int)EBType::VectorPotential; }
    virtual const bool time_varying() const override { return has_function; }

    bool align_with_boundary;
    bool has_function;
    Optional3D1VFunction A_0, A_1, A_2;

    EBType type = EBType::VectorPotential;
    static std::string tag;
};

//-----------------------------------------------------------------------------

class CollectionWall : public FieldBoundaryEB
{
  public:
    CollectionWall();
    CollectionWall(int idx, const sol::table& bc_def);
    ~CollectionWall();

    virtual Vector<Vector<Real>> get_wall_state(const Array<Real, AMREX_SPACEDIM> wall_centre,
                                                const Array<Array<Real, 3>, 3>& wall_coord,
                                                const Real t) const override;

    virtual const int get_type() const override { return (int)EBType::Collection; }
    virtual const bool time_varying() const override;

    Vector<std::unique_ptr<FieldBoundaryEB>> bcs;

    EBType type = EBType::Collection;
    static std::string tag;
};

    //-----------------------------------------------------------------------------

#endif
#endif  // MFP_FIELD_BC_H
