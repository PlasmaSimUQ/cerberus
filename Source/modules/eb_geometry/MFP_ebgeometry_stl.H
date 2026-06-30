#ifndef MFP_EBGEOMETRY_STL_H
#define MFP_EBGEOMETRY_STL_H

// EBGeometry-backed STL signed-distance reader (standalone module).
//
// This class is intentionally independent of the Cerberus / AMReX embedded
// boundary (EB) infrastructure so that STL geometries can be read and queried
// even in fluid-only builds (USE_EB=FALSE).  The EBGeometry include path
// ($EBGEOMETRY_HOME) is added unconditionally by the top-level Makefile, so
// nothing here depends on AMREX_USE_EB.
//
// EBGeometry is intrinsically 3D, so this module is only meaningful for
// AMREX_SPACEDIM > 1 (the surrounding Make.package is included unconditionally,
// which means the 1D build must still compile this translation unit to empty).
//
// This first pass is a faithful relocation of the previous implementation that
// lived in Source/modules/eb/MFP_read_geom.{h,cpp}; the Tier 0 / Tier 1
// memory-footprint reductions will be applied here in subsequent passes.

#if AMREX_SPACEDIM > 1

    #include "sol.hpp"

    #include <AMReX_Array.H>
    #include <AMReX_REAL.H>
    #include <AMReX_SPACE.H>

    // EBGeometry umbrella header (header-only library).
    #include "EBGeometry.hpp"

    #include <memory>
    #include <string>

using namespace amrex;

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

//=============================================================================
// ReadEBGeometrySTL_TriMesh  --  Tier 0 memory-reduced variant.
//
// Identical public API to ReadEBGeometrySTL, but with two changes that shrink
// the resident SDF (and speed up queries) without altering behaviour:
//
//   (1) SDF type FastCompactMeshSDF -> FastTriMeshSDF. The triangle SDF stores
//       self-contained triangles (positions + baked vertex/edge pseudonormals)
//       and does NOT retain the DCEL half-edge/vertex connectivity graph that
//       FastCompactMeshSDF keeps alive.
//   (2) Geometry precision T  Real(double) -> float, halving per-coordinate
//       storage in the triangles and BVH nodes.
//
// Trade-off: loses query-time mesh topology (getClosestTriangles has no
// adjacency) and uses single precision -- neither is exercised by Cerberus,
// which only needs a signed-distance value. See the Tier 0 design notes.
//=============================================================================

class ReadEBGeometrySTL_TriMesh
{
  public:
    // Tier 0 (2): single (float) precision for the geometry / BVH.
    using T = float;

    // Metadata attached to facets. We do not currently use metadata in Cerberus.
    using Meta = int;

    // BVH branching factor. PaintEB uses K = 4.
    static constexpr int K = 4;

    using Vec3 = EBGeometry::Vec3T<T>;
    using BV = EBGeometry::BoundingVolumes::AABBT<T>;

    // Tier 0 (1): self-contained triangle SDF (drops the retained DCEL graph).
    using SDF = EBGeometry::FastTriMeshSDF<T, Meta, BV, K>;

    ReadEBGeometrySTL_TriMesh();

    explicit ReadEBGeometrySTL_TriMesh(const std::string& stl_file);

    ReadEBGeometrySTL_TriMesh(const std::string& stl_file, bool flip_sign);

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

#endif  // AMREX_SPACEDIM > 1
#endif  // MFP_EBGEOMETRY_STL_H
