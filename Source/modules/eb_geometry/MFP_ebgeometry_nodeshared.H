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

    // MPI types (MPI_Comm/MPI_Win) for the shared-memory window. We gate the
    // window path on MPI_VERSION >= 3 (the MPI library's own capability), not on
    // AMReX's optional USE_MPI3/BL_USE_MPI3 build flag: our code calls MPI-3
    // directly, so it only needs an MPI-3 library, regardless of AMReX's flags.
    #if defined(AMREX_USE_MPI)
        #include <mpi.h>
    #endif

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

// On-window layout header, stored at byte offset 0 of the shared buffer. Every
// field is a count or a byte OFFSET (never a pointer), so it is position
// independent and valid in each rank's distinct mapping of the same memory.
struct SharedSDFHeader {
    std::int64_t magic;
    std::int64_t n_nodes;
    std::int64_t n_tris;
    std::int64_t nodes_off;
    std::int64_t tris_off;
    std::int64_t total_bytes;
};

static constexpr std::int64_t SHARED_SDF_MAGIC = 0x4d46505342564831LL;  // "MFPSBVH1"

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

// ---------------------------------------------------------------------------
// Tier 1 step 2: node-shared SDF.
//
// The flat arrays live in an MPI-3 shared-memory window, so all ranks on a node
// map ONE physical copy that the node-lead builds once. Falls back to a private
// heap copy when MPI-3 is unavailable, OR when there is a single rank per node
// (m_node_size == 1) where a window would share with nobody. Query reuses
// mfp_ebgeom::flat_query unchanged.
//
// COLLECTIVE: build_shared() and free_shared() must be called by ALL ranks of
// the node together (guaranteed via SPMD read_config / the collective self-test).
// Stores the canonical signed distance (no flip; flip is a per-handle concern in
// step 3).
// ---------------------------------------------------------------------------
class NodeSharedTriMeshSDF
{
  public:
    NodeSharedTriMeshSDF() = default;
    ~NodeSharedTriMeshSDF();

    NodeSharedTriMeshSDF(const NodeSharedTriMeshSDF&) = delete;
    NodeSharedTriMeshSDF& operator=(const NodeSharedTriMeshSDF&) = delete;

    // Collective. Build (lead) + map (all) the shared SDF for stl_file.
    void build_shared(const std::string& stl_file);

    // Collective. Release the window/communicator. Safe to call more than once.
    void free_shared();

    // Canonical signed distance (no flip). Read-only; concurrent-safe.
    Real query(AMREX_D_DECL(Real x, Real y, Real z)) const;

    static void register_with_lua(sol::state& lua);

    // Collective validation: correctness vs FlatTriMeshSDF (IO rank), cross-rank
    // agreement (Allreduce), and a per-node allocation report.
    static void self_test(const std::string& stl_file, int n);

  private:
    void build_heap(const std::string& stl_file);  // private-copy path

    const mfp_ebgeom::FlatBVHNode* m_nodes = nullptr;
    const mfp_ebgeom::SdfTriangle* m_tris = nullptr;
    std::int64_t m_n_nodes = 0;
    std::int64_t m_n_tris = 0;

    int m_node_rank = 0;
    int m_node_size = 1;
    bool m_using_window = false;
    std::int64_t m_alloc_bytes = 0;
    std::string m_filename;

    // fallback / single-rank-shortcut private storage
    std::vector<mfp_ebgeom::FlatBVHNode> m_nodes_heap;
    std::vector<mfp_ebgeom::SdfTriangle> m_tris_heap;

#if defined(AMREX_USE_MPI) && (MPI_VERSION >= 3)
    char* m_base = nullptr;
    MPI_Comm m_node_comm = MPI_COMM_NULL;
    MPI_Win m_win = MPI_WIN_NULL;
#endif
};

#endif  // AMREX_SPACEDIM > 1
#endif  // MFP_EBGEOMETRY_NODESHARED_H
