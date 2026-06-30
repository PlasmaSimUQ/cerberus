#ifndef MFP_EBGEOMETRY_NODESHARED_H
#define MFP_EBGEOMETRY_NODESHARED_H

// Node-shared EBGeometry STL/SDF support.
//
// STEP 1 (this file, no MPI): a flat, pointer-free, cache-friendly BVH over a
// triangle soup, queried with the same ordered nearest-first descent + pruning
// that EBGeometry's LinearBVH uses. The flat layout (contiguous POD arrays,
// integer child indices) is what later (Tier 1, Steps 2-3) lets the data live in
// an MPI-3 shared-memory window so a node holds ONE copy instead of one per rank.
//
// The free functions below (build_triangle_soup / build_flat_bvh / flat_query)
// are deliberately split out and operate on raw pointers so the exact same build
// and query code is reused by the shared-memory backend in later steps.

#if AMREX_SPACEDIM > 1

    #include "sol.hpp"

    #include <AMReX_REAL.H>
    #include <AMReX_SPACE.H>

    // EBGeometry umbrella header (header-only; include path is unconditional).
    #include "EBGeometry.hpp"

    #include <cstdint>
    #include <string>
    #include <vector>

using namespace amrex;

namespace mfp_ebgeom
{

// Geometry precision (Tier 0: float) and primitive type. Kept identical to
// ReadEBGeometrySTL_TriMesh so results match.
using SdfT = float;
using SdfMeta = int;
using SdfVec3 = EBGeometry::Vec3T<SdfT>;
using SdfTriangle = EBGeometry::Triangle<SdfT, SdfMeta>;

// 4-ary tree (matches EBGeometry's K=4) with a small leaf cutoff.
static constexpr int FLATBVH_K = 4;
static constexpr int FLATBVH_LEAF_MAX = 4;

// Pointer-free BVH node. child[] holds node INDICES (depth-first array), never
// pointers, so the array is relocatable into a shared-memory window. A leaf has
// prim_count > 0 (and prim_offset into the triangle array); an internal node has
// prim_count == 0 and up to FLATBVH_K child indices (-1 = absent).
struct FlatBVHNode {
    float lo[3];
    float hi[3];
    std::int32_t child[FLATBVH_K];
    std::int32_t prim_offset;
    std::int32_t prim_count;
};

// Build a triangle soup IDENTICAL to EBGeometry::FastTriMeshSDF: read the STL into
// a DCEL mesh and copy face normal + vertex positions + vertex pseudonormals +
// EDGE pseudonormals (from the edges, not the vertices) into each triangle. This
// is why we do not use Parser::readIntoTriangles, which sets edge normals to the
// vertex normals and would not match ReadEBGeometrySTL_TriMesh near edges.
std::vector<SdfTriangle> build_triangle_soup(const std::string& stl_file);

// Build a contiguous, depth-first, K=4 BVH over the soup. nodes_out is the flat
// node array; tris_out is the triangles reordered so each leaf's triangles are
// contiguous (referenced by prim_offset/prim_count).
void build_flat_bvh(const std::vector<SdfTriangle>& soup,
                    std::vector<FlatBVHNode>& nodes_out,
                    std::vector<SdfTriangle>& tris_out);

// Signed distance at p. Stack-based ordered nearest-first descent with AABB
// pruning, reusing EBGeometry::Triangle::signedDistance per primitive. Operates
// on raw pointers so it serves both the heap arrays (Step 1) and a shared window
// (Steps 2-3). Returns the canonical signed distance (no flip applied).
SdfT flat_query(const FlatBVHNode* nodes,
                std::int64_t n_nodes,
                const SdfTriangle* tris,
                const SdfVec3& p);

}  // namespace mfp_ebgeom

// ---------------------------------------------------------------------------
// Step 1 heap-backed reader + validation/benchmark scaffold.
// Registered to Lua as usertype "FlatTriMeshSDF" (same API as the other readers)
// plus a free function "flat_trimesh_self_test(file, n)".
// ---------------------------------------------------------------------------
class FlatTriMeshSDF
{
  public:
    FlatTriMeshSDF();
    explicit FlatTriMeshSDF(const std::string& stl_file);
    FlatTriMeshSDF(const std::string& stl_file, bool flip_sign);

    void read_file(const std::string& stl_file);

    Real query(AMREX_D_DECL(Real x, Real y, Real z)) const;

    const std::string str() const;

    void set_flip_sign(bool flip_sign);
    bool get_flip_sign() const;

    static void register_with_lua(sol::state& lua);

    // Correctness + speed check vs ReadEBGeometrySTL_TriMesh on an n^3 grid over
    // the (padded) mesh bounding box. Prints max/mean |error| and query timings.
    static void self_test(const std::string& stl_file, int n);

  private:
    std::vector<mfp_ebgeom::FlatBVHNode> m_nodes;
    std::vector<mfp_ebgeom::SdfTriangle> m_tris;
    std::string m_filename;
    bool m_flip_sign = false;
};

#endif  // AMREX_SPACEDIM > 1
#endif  // MFP_EBGEOMETRY_NODESHARED_H
