#include "AABB.h"
#include "EBGeometry.hpp"
#include "MFP_polyspline.H"
#include "sol.hpp"

#include <AMReX_Array.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
using namespace amrex;

class ReadSTL
{
  public:
    ReadSTL();
    ReadSTL(const std::string& stl_file);

    void read_file(const std::string& stl_file);

    bool valid_tri(const Array<Real, 9>& facet);

    void get_limits(const Array<Real, 9>& facet, Vector<Real>& upper, Vector<Real>& lower);

    void grow_tree();

    Real query(AMREX_D_DECL(Real x, Real y, Real z));

    const std::string str() const;

    static void register_with_lua(sol::state& lua);

    Vector<Array<Real, 9>> facets;

    aabb::Tree tree;
};

//=============================================================================
// EBGeometry-backed STL signed-distance reader
//=============================================================================

class ReadEBGeometrySTL
{
  public:
    // EBGeometry precision. Keep this as Real to match AMReX/Cerberus builds.
    using T = Real;

    // Metadata attached to facets. We do not currently use metadata in Cerberus.
    using Meta = int;

    // BVH branching factor. PaintEB uses K = 4.
    static constexpr int K = 4;

    using Vec3 = EBGeometry::Vec3T<T>;
    using BV = EBGeometry::BoundingVolumes::AABBT<T>;
    using SDF = EBGeometry::FastCompactMeshSDF<T, Meta, BV, K>;

    ReadEBGeometrySTL();

    explicit ReadEBGeometrySTL(const std::string& stl_file);

    ReadEBGeometrySTL(const std::string& stl_file, bool flip_sign);

    void read_file(const std::string& stl_file);

    Real query(AMREX_D_DECL(Real x, Real y, Real z)) const;

    const std::string str() const;

    void set_flip_sign(bool flip_sign);

    bool get_flip_sign() const;

    static void register_with_lua(sol::state& lua);

  private:
    std::shared_ptr<SDF> m_sdf;
    std::string m_filename;
    bool m_flip_sign = false;
};

//=====
