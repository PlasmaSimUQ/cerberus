#include "MFP_ebgeometry_nodeshared.H"

#if AMREX_SPACEDIM > 1

    #include "MFP_ebgeometry_stl.H"  // ReadEBGeometrySTL_TriMesh (golden reference)

    #include <AMReX.H>
    #include <AMReX_BLProfiler.H>
    #include <AMReX_ParallelDescriptor.H>
    #include <AMReX_Print.H>
    #include <algorithm>
    #include <cmath>
    #include <cstdint>
    #include <cstring>
    #include <limits>
    #include <sstream>

namespace mfp_ebgeom
{

namespace
{

constexpr float INF = std::numeric_limits<float>::max();

// Distance from point (px,py,pz) to an AABB; 0 if the point is inside the box.
inline float aabb_min_dist(const FlatBVHNode& n, float px, float py, float pz)
{
    const float dx = std::fmax(std::fmax(n.lo[0] - px, px - n.hi[0]), 0.0f);
    const float dy = std::fmax(std::fmax(n.lo[1] - py, py - n.hi[1]), 0.0f);
    const float dz = std::fmax(std::fmax(n.lo[2] - pz, pz - n.hi[2]), 0.0f);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Per-triangle build record (AABB + centroid + index into the soup).
struct BuildItem {
    float lo[3];
    float hi[3];
    float c[3];
    int idx;
};

// Recursive 4-ary top-down BVH builder. Emits nodes depth-first and copies each
// leaf's triangles contiguously into tris.
struct FlatBVHBuilder {
    const std::vector<SdfTriangle>& soup;
    std::vector<BuildItem>& items;
    std::vector<FlatBVHNode>& nodes;
    std::vector<SdfTriangle>& tris;

    // Count-median (object median) split of items[begin,end) along the longest
    // centroid axis. Count-based => each side is strictly smaller, guaranteeing
    // termination even with coincident centroids.
    int median_split(int begin, int end)
    {
        float clo[3] = {INF, INF, INF};
        float chi[3] = {-INF, -INF, -INF};
        for (int i = begin; i < end; ++i) {
            for (int d = 0; d < 3; ++d) {
                clo[d] = std::fmin(clo[d], items[i].c[d]);
                chi[d] = std::fmax(chi[d], items[i].c[d]);
            }
        }
        int axis = 0;
        float ext = chi[0] - clo[0];
        for (int d = 1; d < 3; ++d) {
            const float e = chi[d] - clo[d];
            if (e > ext) {
                ext = e;
                axis = d;
            }
        }
        const int mid = (begin + end) / 2;
        std::nth_element(
          items.begin() + begin,
          items.begin() + mid,
          items.begin() + end,
          [axis](const BuildItem& a, const BuildItem& b) { return a.c[axis] < b.c[axis]; });
        return mid;
    }

    int recurse(int begin, int end)
    {
        const int my = static_cast<int>(nodes.size());
        nodes.push_back(FlatBVHNode {});  // reserve depth-first slot

        // node AABB = union of item AABBs
        FlatBVHNode node {};
        for (int d = 0; d < 3; ++d) {
            node.lo[d] = INF;
            node.hi[d] = -INF;
        }
        for (int i = begin; i < end; ++i) {
            for (int d = 0; d < 3; ++d) {
                node.lo[d] = std::fmin(node.lo[d], items[i].lo[d]);
                node.hi[d] = std::fmax(node.hi[d], items[i].hi[d]);
            }
        }

        if (end - begin <= FLATBVH_LEAF_MAX) {
            node.prim_offset = static_cast<std::int32_t>(tris.size());
            node.prim_count = static_cast<std::int32_t>(end - begin);
            for (int k = 0; k < FLATBVH_K; ++k) node.child[k] = -1;
            for (int i = begin; i < end; ++i) tris.push_back(soup[items[i].idx]);
        } else {
            node.prim_offset = -1;
            node.prim_count = 0;
            // up to 4 buckets via 3 median splits
            const int b2 = median_split(begin, end);
            const int b1 = median_split(begin, b2);
            const int b3 = median_split(b2, end);
            const int bnd[FLATBVH_K + 1] = {begin, b1, b2, b3, end};
            for (int k = 0; k < FLATBVH_K; ++k) {
                node.child[k] = (bnd[k] < bnd[k + 1]) ? recurse(bnd[k], bnd[k + 1]) : -1;
            }
        }

        nodes[my] = node;  // write back after children exist (indexed by int -> safe)
        return my;
    }
};

}  // namespace

std::vector<SdfTriangle> build_triangle_soup(const std::string& stl_file)
{
    BL_PROFILE("mfp_ebgeom::build_triangle_soup");

    // Read into a DCEL mesh (computes vertex/edge pseudonormals), then flatten to
    // self-contained triangles exactly as EBGeometry::FastTriMeshSDF does. The
    // DCEL mesh is released when this function returns.
    auto mesh = EBGeometry::Parser::readIntoDCEL<SdfT, SdfMeta>(stl_file);

    std::vector<SdfTriangle> soup;
    const auto& faces = mesh->getFaces();
    soup.reserve(faces.size());

    for (const auto& f : faces) {
        const auto vertices = f->gatherVertices();
        const auto edges = f->gatherEdges();
        if (vertices.size() != 3 || edges.size() != 3) continue;  // STL is all triangles

        SdfTriangle tri;
        tri.setNormal(f->getNormal());
        tri.setVertexPositions(
          {vertices[0]->getPosition(), vertices[1]->getPosition(), vertices[2]->getPosition()});
        tri.setVertexNormals(
          {vertices[0]->getNormal(), vertices[1]->getNormal(), vertices[2]->getNormal()});
        tri.setEdgeNormals({edges[0]->getNormal(), edges[1]->getNormal(), edges[2]->getNormal()});
        tri.setMetaData(SdfMeta(0));

        soup.push_back(tri);
    }

    return soup;
}

void build_flat_bvh(const std::vector<SdfTriangle>& soup,
                    std::vector<FlatBVHNode>& nodes_out,
                    std::vector<SdfTriangle>& tris_out)
{
    BL_PROFILE("mfp_ebgeom::build_flat_bvh");

    nodes_out.clear();
    tris_out.clear();
    if (soup.empty()) return;

    std::vector<BuildItem> items(soup.size());
    for (std::size_t i = 0; i < soup.size(); ++i) {
        const auto& V = soup[i].getVertexPositions();
        for (int d = 0; d < 3; ++d) {
            const float v0 = V[0][d];
            const float v1 = V[1][d];
            const float v2 = V[2][d];
            items[i].lo[d] = std::fmin(v0, std::fmin(v1, v2));
            items[i].hi[d] = std::fmax(v0, std::fmax(v1, v2));
            items[i].c[d] = (v0 + v1 + v2) / 3.0f;
        }
        items[i].idx = static_cast<int>(i);
    }

    nodes_out.reserve(2 * soup.size() / FLATBVH_LEAF_MAX + 8);
    tris_out.reserve(soup.size());

    FlatBVHBuilder builder {soup, items, nodes_out, tris_out};
    builder.recurse(0, static_cast<int>(soup.size()));
}

SdfT flat_query(const FlatBVHNode* nodes,
                std::int64_t n_nodes,
                const SdfTriangle* tris,
                const SdfVec3& p)
{
    if (nodes == nullptr || n_nodes <= 0) return INF;

    const float px = p[0];
    const float py = p[1];
    const float pz = p[2];

    float best = INF;      // signed distance of nearest triangle so far
    float best_mag = INF;  // its magnitude (drives pruning)

    // Fixed-size local stack (no heap). Depth ~ log4(N); 128 is ample.
    int stack[128];
    int sp = 0;
    stack[sp++] = 0;  // root

    while (sp > 0) {
        const FlatBVHNode& n = nodes[stack[--sp]];

        // Re-check prune at pop: best_mag may have tightened since we pushed it.
        if (aabb_min_dist(n, px, py, pz) >= best_mag) continue;

        if (n.prim_count > 0) {  // leaf
            const int e = n.prim_offset + n.prim_count;
            for (int t = n.prim_offset; t < e; ++t) {
                const float d = tris[t].signedDistance(p);
                const float m = std::fabs(d);
                if (m < best_mag) {
                    best_mag = m;
                    best = d;
                }
            }
        } else {  // internal: ordered nearest-first descent
            struct Child {
                float dist;
                int idx;
            } c[FLATBVH_K];
            int nc = 0;
            for (int k = 0; k < FLATBVH_K; ++k) {
                const int ci = n.child[k];
                if (ci < 0) continue;
                c[nc].idx = ci;
                c[nc].dist = aabb_min_dist(nodes[ci], px, py, pz);
                ++nc;
            }
            // Sort by BV distance DESCENDING so the nearest child is pushed last
            // and therefore popped first (LIFO) -> tightens best_mag early.
            std::sort(c, c + nc, [](const Child& a, const Child& b) { return a.dist > b.dist; });
            for (int k = 0; k < nc; ++k) {
                if (c[k].dist < best_mag && sp < 128) stack[sp++] = c[k].idx;
            }
        }
    }

    return best;
}

}  // namespace mfp_ebgeom

// ===========================================================================
// FlatTriMeshSDF
// ===========================================================================

FlatTriMeshSDF::FlatTriMeshSDF() {}

FlatTriMeshSDF::FlatTriMeshSDF(const std::string& stl_file) { read_file(stl_file); }

FlatTriMeshSDF::FlatTriMeshSDF(const std::string& stl_file, bool flip_sign) : m_flip_sign(flip_sign)
{
    read_file(stl_file);
}

void FlatTriMeshSDF::read_file(const std::string& stl_file)
{
    BL_PROFILE("FlatTriMeshSDF::read_file");

    m_filename = stl_file;

    auto soup = mfp_ebgeom::build_triangle_soup(stl_file);
    if (soup.empty()) {
        amrex::Abort("FlatTriMeshSDF: no triangles read from STL file '" + stl_file + "'");
    }
    mfp_ebgeom::build_flat_bvh(soup, m_nodes, m_tris);
}

Real FlatTriMeshSDF::query(AMREX_D_DECL(Real x, Real y, Real z)) const
{
    BL_PROFILE("FlatTriMeshSDF::query");

    if (m_nodes.empty()) {
        amrex::Abort("FlatTriMeshSDF::query called before a valid STL file was loaded");
    }

    const mfp_ebgeom::SdfT px = static_cast<mfp_ebgeom::SdfT>(x);
    const mfp_ebgeom::SdfT py = static_cast<mfp_ebgeom::SdfT>(y);
    #if AMREX_SPACEDIM == 3
    const mfp_ebgeom::SdfT pz = static_cast<mfp_ebgeom::SdfT>(z);
    #else
    const mfp_ebgeom::SdfT pz = static_cast<mfp_ebgeom::SdfT>(0);
    #endif

    const mfp_ebgeom::SdfVec3 p(px, py, pz);

    const Real d =
      static_cast<Real>(mfp_ebgeom::flat_query(m_nodes.data(),
                                               static_cast<std::int64_t>(m_nodes.size()),
                                               m_tris.data(),
                                               p));

    return m_flip_sign ? -d : d;
}

const std::string FlatTriMeshSDF::str() const
{
    BL_PROFILE("FlatTriMeshSDF::str");

    std::stringstream ss;
    ss << "FlatTriMeshSDF\n";
    ss << "  filename  : " << m_filename << "\n";
    ss << "  n_nodes   : " << m_nodes.size() << "\n";
    ss << "  n_tris    : " << m_tris.size() << "\n";
    ss << "  flip_sign : " << m_flip_sign << "\n";
    return ss.str();
}

void FlatTriMeshSDF::set_flip_sign(bool flip_sign) { m_flip_sign = flip_sign; }

bool FlatTriMeshSDF::get_flip_sign() const { return m_flip_sign; }

void FlatTriMeshSDF::register_with_lua(sol::state& lua)
{
    BL_PROFILE("FlatTriMeshSDF::register_with_lua");

    lua.new_usertype<FlatTriMeshSDF>("FlatTriMeshSDF",
                                     sol::constructors<FlatTriMeshSDF(const std::string&),
                                                       FlatTriMeshSDF(const std::string&, bool)>(),
                                     "query",
                                     &FlatTriMeshSDF::query,
                                     "set_flip_sign",
                                     &FlatTriMeshSDF::set_flip_sign,
                                     "get_flip_sign",
                                     &FlatTriMeshSDF::get_flip_sign,
                                     "str",
                                     &FlatTriMeshSDF::str);

#ifdef AMREX_DEBUG
    lua.set_function("flat_trimesh_self_test", &FlatTriMeshSDF::self_test);
#endif
}

#ifdef AMREX_DEBUG
void FlatTriMeshSDF::self_test(const std::string& stl_file, int n)
{
    BL_PROFILE("FlatTriMeshSDF::self_test");

    const int N = std::max(2, n);

    amrex::Print() << "[FlatTriMeshSDF::self_test] file='" << stl_file << "' grid=" << N << "^"
                   << AMREX_SPACEDIM << "\n";

    FlatTriMeshSDF flat(stl_file);
    ReadEBGeometrySTL_TriMesh golden(stl_file);

    // Sample box = padded mesh bounding box (from the flat BVH root node).
    const auto& root = flat.m_nodes[0];
    double lo[3], hi[3];
    for (int d = 0; d < 3; ++d) {
        const double L = root.lo[d];
        const double H = root.hi[d];
        const double pad = 0.25 * (H - L) + 1e-6;
        lo[d] = L - pad;
        hi[d] = H + pad;
    }

    const int NK = (AMREX_SPACEDIM == 3) ? N : 1;

    auto coord = [&](int i, int d) { return lo[d] + (i + 0.5) * (hi[d] - lo[d]) / N; };

    // --- correctness ---
    double max_err = 0.0, sum_err = 0.0;
    long cnt = 0;
    for (int k = 0; k < NK; ++k) {
        const double z = (AMREX_SPACEDIM == 3) ? coord(k, 2) : 0.0;
        amrex::ignore_unused(z);
        for (int j = 0; j < N; ++j) {
            const double y = coord(j, 1);
            for (int i = 0; i < N; ++i) {
                const double x = coord(i, 0);
                const double a = golden.query(AMREX_D_DECL(x, y, z));
                const double b = flat.query(AMREX_D_DECL(x, y, z));
                const double e = std::abs(a - b);
                max_err = std::max(max_err, e);
                sum_err += e;
                ++cnt;
            }
        }
    }

    // --- timing (separate passes so neither caches the other's work) ---
    volatile double sink = 0.0;
    double t0 = amrex::ParallelDescriptor::second();
    for (int k = 0; k < NK; ++k) {
        const double z = (AMREX_SPACEDIM == 3) ? coord(k, 2) : 0.0;
        amrex::ignore_unused(z);
        for (int j = 0; j < N; ++j) {
            const double y = coord(j, 1);
            for (int i = 0; i < N; ++i) sink += golden.query(AMREX_D_DECL(coord(i, 0), y, z));
        }
    }
    const double t_gold = amrex::ParallelDescriptor::second() - t0;

    t0 = amrex::ParallelDescriptor::second();
    for (int k = 0; k < NK; ++k) {
        const double z = (AMREX_SPACEDIM == 3) ? coord(k, 2) : 0.0;
        amrex::ignore_unused(z);
        for (int j = 0; j < N; ++j) {
            const double y = coord(j, 1);
            for (int i = 0; i < N; ++i) sink += flat.query(AMREX_D_DECL(coord(i, 0), y, z));
        }
    }
    const double t_flat = amrex::ParallelDescriptor::second() - t0;

    double extent = 0.0;
    for (int d = 0; d < 3; ++d) extent = std::max(extent, hi[d] - lo[d]);

    amrex::Print() << "  triangles    : " << flat.m_tris.size()
                   << ", nodes : " << flat.m_nodes.size() << "\n"
                   << "  queries      : " << cnt << "\n"
                   << "  max |err|    : " << max_err << "  (" << (max_err / extent)
                   << " of extent)\n"
                   << "  mean |err|   : " << (cnt ? sum_err / cnt : 0.0) << "\n"
                   << "  golden time  : " << t_gold << " s  (" << (cnt / t_gold) << " q/s)\n"
                   << "  flat   time  : " << t_flat << " s  (" << (cnt / t_flat) << " q/s)\n"
                   << "  speed ratio  : " << (t_flat > 0 ? t_gold / t_flat : 0.0)
                   << "x (flat vs golden)\n";
}
#endif  // AMREX_DEBUG

// ===========================================================================
// NodeSharedTriMeshSDF (Tier 1 step 2: MPI-3 shared-memory window)
// ===========================================================================

namespace {

#if defined(AMREX_USE_MPI) && (MPI_VERSION >= 3)

inline std::int64_t align16(std::int64_t x)
{
    return (x + 15) & ~static_cast<std::int64_t>(15);
}

// Compute the contiguous window layout: [Header | pad | Nodes | pad | Tris].
mfp_ebgeom::SharedSDFHeader make_layout(std::int64_t n_nodes, std::int64_t n_tris)
{
    mfp_ebgeom::SharedSDFHeader h{};
    h.magic = mfp_ebgeom::SHARED_SDF_MAGIC;
    h.n_nodes = n_nodes;
    h.n_tris = n_tris;

    std::int64_t off = align16(static_cast<std::int64_t>(sizeof(mfp_ebgeom::SharedSDFHeader)));
    h.nodes_off = off;
    off += n_nodes * static_cast<std::int64_t>(sizeof(mfp_ebgeom::FlatBVHNode));
    off = align16(off);
    h.tris_off = off;
    off += n_tris * static_cast<std::int64_t>(sizeof(mfp_ebgeom::SdfTriangle));
    h.total_bytes = off;
    return h;
}

// Single-rank-per-node shortcut warning (emitted once, on the IO rank).
void warn_single_rank_once()
{
    static bool warned = false;
    if (!warned) {
        warned = true;
        amrex::Print() << "[NodeSharedTriMeshSDF] Single-rank-per-node: shortcut to the heap "
                          "path (no shared-memory window to build).\n";
    }
}

#endif  // AMREX_USE_MPI && MPI_VERSION >= 3

}  // namespace

void NodeSharedTriMeshSDF::build_heap(const std::string& stl_file)
{
    auto soup = mfp_ebgeom::build_triangle_soup(stl_file);
    if (soup.empty()) {
        amrex::Abort("NodeSharedTriMeshSDF: no triangles read from STL file '" + stl_file + "'");
    }
    mfp_ebgeom::build_flat_bvh(soup, m_nodes_heap, m_tris_heap);

    m_nodes = m_nodes_heap.data();
    m_tris = m_tris_heap.data();
    m_n_nodes = static_cast<std::int64_t>(m_nodes_heap.size());
    m_n_tris = static_cast<std::int64_t>(m_tris_heap.size());
    m_using_window = false;
    m_alloc_bytes = m_n_nodes * static_cast<std::int64_t>(sizeof(mfp_ebgeom::FlatBVHNode)) +
                    m_n_tris * static_cast<std::int64_t>(sizeof(mfp_ebgeom::SdfTriangle));
}

void NodeSharedTriMeshSDF::build_shared(const std::string& stl_file)
{
    BL_PROFILE("NodeSharedTriMeshSDF::build_shared");
    m_filename = stl_file;

#if defined(AMREX_USE_MPI) && (MPI_VERSION >= 3)
    // node-local communicator: the ranks that can share physical memory
    BL_MPI_REQUIRE(MPI_Comm_split_type(amrex::ParallelDescriptor::Communicator(),
                                       MPI_COMM_TYPE_SHARED,
                                       0,
                                       MPI_INFO_NULL,
                                       &m_node_comm));
    BL_MPI_REQUIRE(MPI_Comm_rank(m_node_comm, &m_node_rank));
    BL_MPI_REQUIRE(MPI_Comm_size(m_node_comm, &m_node_size));

    if (m_node_size <= 1) {
        // Single rank on this node: a window would share with nobody. Shortcut to
        // a private heap copy (see Decision 1) and warn once.
        warn_single_rank_once();
        BL_MPI_REQUIRE(MPI_Comm_free(&m_node_comm));
        m_node_comm = MPI_COMM_NULL;
        build_heap(stl_file);
        return;
    }

    // (1) lead builds the flat arrays; peers idle (build is paid once per node)
    std::vector<mfp_ebgeom::FlatBVHNode> nodes;
    std::vector<mfp_ebgeom::SdfTriangle> tris;
    mfp_ebgeom::SharedSDFHeader hdr{};
    if (m_node_rank == 0) {
        auto soup = mfp_ebgeom::build_triangle_soup(stl_file);
        if (soup.empty()) {
            amrex::Abort("NodeSharedTriMeshSDF: no triangles read from STL file '" + stl_file +
                         "'");
        }
        mfp_ebgeom::build_flat_bvh(soup, nodes, tris);
        hdr = make_layout(static_cast<std::int64_t>(nodes.size()),
                          static_cast<std::int64_t>(tris.size()));
    }

    // (2) allocate the shared window: lead requests all bytes, peers request 0
    const MPI_Aint bytes = (m_node_rank == 0) ? static_cast<MPI_Aint>(hdr.total_bytes) : 0;
    BL_MPI_REQUIRE(MPI_Win_allocate_shared(bytes, 1, MPI_INFO_NULL, m_node_comm, &m_base, &m_win));
    BL_MPI_REQUIRE(
      MPI_Win_lock_all(MPI_MODE_NOCHECK, m_win));  // passive epoch for the window's life

    // (3) everyone maps the lead's segment (same physical bytes, per-rank address)
    MPI_Aint qsize = 0;
    int qdisp = 0;
    BL_MPI_REQUIRE(MPI_Win_shared_query(m_win, 0, &qsize, &qdisp, &m_base));

    // (4) lead writes header + arrays via direct store into the mapped buffer
    if (m_node_rank == 0) {
        auto* H = reinterpret_cast<mfp_ebgeom::SharedSDFHeader*>(m_base);
        *H = hdr;
        std::memcpy(m_base + hdr.nodes_off,
                    nodes.data(),
                    static_cast<std::size_t>(hdr.n_nodes) * sizeof(mfp_ebgeom::FlatBVHNode));
        std::memcpy(m_base + hdr.tris_off,
                    tris.data(),
                    static_cast<std::size_t>(hdr.n_tris) * sizeof(mfp_ebgeom::SdfTriangle));
    }
    // lead's heap arrays (nodes/tris) freed here; only the window copy remains

    // (5) publish: make the lead's stores visible to all readers
    BL_MPI_REQUIRE(MPI_Win_sync(m_win));
    BL_MPI_REQUIRE(MPI_Barrier(m_node_comm));
    BL_MPI_REQUIRE(MPI_Win_sync(m_win));

    // (6) all ranks set query pointers from the now-visible header
    const auto* H = reinterpret_cast<const mfp_ebgeom::SharedSDFHeader*>(m_base);
    if (H->magic != mfp_ebgeom::SHARED_SDF_MAGIC) {
        amrex::Abort(
          "NodeSharedTriMeshSDF: shared header magic mismatch (memory-model / sync error)");
    }
    m_n_nodes = H->n_nodes;
    m_n_tris = H->n_tris;
    m_nodes = reinterpret_cast<const mfp_ebgeom::FlatBVHNode*>(m_base + H->nodes_off);
    m_tris = reinterpret_cast<const mfp_ebgeom::SdfTriangle*>(m_base + H->tris_off);
    m_using_window = true;
    m_alloc_bytes = (m_node_rank == 0) ? hdr.total_bytes : 0;
#else
    // No MPI-3 available: private heap copy per process.
    build_heap(stl_file);
#endif
}

void NodeSharedTriMeshSDF::free_shared()
{
#if defined(AMREX_USE_MPI) && (MPI_VERSION >= 3)
    if (m_win != MPI_WIN_NULL) {
        BL_MPI_REQUIRE(MPI_Barrier(m_node_comm));
        BL_MPI_REQUIRE(MPI_Win_unlock_all(m_win));
        BL_MPI_REQUIRE(MPI_Win_free(&m_win));  // sets m_win = MPI_WIN_NULL
        m_base = nullptr;
    }
    if (m_node_comm != MPI_COMM_NULL) {
        BL_MPI_REQUIRE(MPI_Comm_free(&m_node_comm));  // sets m_node_comm = MPI_COMM_NULL
    }
#endif
    m_nodes = nullptr;
    m_tris = nullptr;
    m_n_nodes = 0;
    m_n_tris = 0;
    m_using_window = false;
    m_nodes_heap.clear();
    m_nodes_heap.shrink_to_fit();
    m_tris_heap.clear();
    m_tris_heap.shrink_to_fit();
}

NodeSharedTriMeshSDF::~NodeSharedTriMeshSDF() { free_shared(); }

Real NodeSharedTriMeshSDF::query(AMREX_D_DECL(Real x, Real y, Real z)) const
{
    BL_PROFILE("NodeSharedTriMeshSDF::query");

    if (m_nodes == nullptr) {
        amrex::Abort("NodeSharedTriMeshSDF::query called before build_shared()");
    }

    const mfp_ebgeom::SdfT px = static_cast<mfp_ebgeom::SdfT>(x);
    const mfp_ebgeom::SdfT py = static_cast<mfp_ebgeom::SdfT>(y);
    #if AMREX_SPACEDIM == 3
    const mfp_ebgeom::SdfT pz = static_cast<mfp_ebgeom::SdfT>(z);
    #else
    const mfp_ebgeom::SdfT pz = static_cast<mfp_ebgeom::SdfT>(0);
    #endif
    const mfp_ebgeom::SdfVec3 p(px, py, pz);

    return static_cast<Real>(mfp_ebgeom::flat_query(m_nodes, m_n_nodes, m_tris, p));
}

// --- filename-keyed registry (step 3) --------------------------------------

std::map<std::string, std::shared_ptr<NodeSharedTriMeshSDF>> NodeSharedTriMeshSDF::s_registry;

std::shared_ptr<NodeSharedTriMeshSDF>
NodeSharedTriMeshSDF::get_or_create(const std::string& stl_file)
{
    BL_PROFILE("NodeSharedTriMeshSDF::get_or_create");

    auto it = s_registry.find(stl_file);
    if (it != s_registry.end()) return it->second;  // dedup; no MPI call on a hit

    auto res = std::make_shared<NodeSharedTriMeshSDF>();
    res->build_shared(stl_file);  // COLLECTIVE on first use
    s_registry.emplace(stl_file, res);
    return res;
}

void NodeSharedTriMeshSDF::clear_all()
{
    BL_PROFILE("NodeSharedTriMeshSDF::clear_all");

    // std::map iterates in sorted (filename) order, identical on every rank, so
    // the collective free_shared() calls line up across ranks.
    for (auto& kv : s_registry) {
        if (kv.second) kv.second->free_shared();
    }
    s_registry.clear();
}

void NodeSharedTriMeshSDF::register_with_lua(sol::state& lua)
{
    BL_PROFILE("NodeSharedTriMeshSDF::register_with_lua");
#ifdef AMREX_DEBUG
    lua.set_function("node_shared_self_test", &NodeSharedTriMeshSDF::self_test);
#else
    amrex::ignore_unused(lua);
#endif
}

#ifdef AMREX_DEBUG
void NodeSharedTriMeshSDF::self_test(const std::string& stl_file, int n)
{
    BL_PROFILE("NodeSharedTriMeshSDF::self_test");

    const int N = std::max(2, n);

    amrex::Print() << "[NodeSharedTriMeshSDF::self_test] file='" << stl_file << "' grid=" << N << "^"
                   << AMREX_SPACEDIM << "\n";

    NodeSharedTriMeshSDF shared;
    shared.build_shared(stl_file);  // collective

    // Sample box = padded mesh bounding box (root node AABB; identical on all ranks).
    const mfp_ebgeom::FlatBVHNode& root = shared.m_nodes[0];
    double lo[3], hi[3];
    for (int d = 0; d < 3; ++d) {
        const double L = root.lo[d];
        const double H = root.hi[d];
        const double pad = 0.25 * (H - L) + 1e-6;
        lo[d] = L - pad;
        hi[d] = H + pad;
    }

    const int NK = (AMREX_SPACEDIM == 3) ? N : 1;
    auto coord = [&](int i, int d) { return lo[d] + (i + 0.5) * (hi[d] - lo[d]) / N; };

    // (1) value correctness on the global IO rank vs a Step-1 FlatTriMeshSDF.
    // Expect EXACTLY 0: the lead built the same arrays via the same build_flat_bvh.
    double max_err = 0.0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        FlatTriMeshSDF flat(stl_file);
        for (int k = 0; k < NK; ++k) {
            const double z = (AMREX_SPACEDIM == 3) ? coord(k, 2) : 0.0;
            amrex::ignore_unused(z);
            for (int j = 0; j < N; ++j) {
                const double y = coord(j, 1);
                for (int i = 0; i < N; ++i) {
                    const double x = coord(i, 0);
                    const double a = flat.query(AMREX_D_DECL(x, y, z));
                    const double b = shared.query(AMREX_D_DECL(x, y, z));
                    max_err = std::max(max_err, std::abs(a - b));
                }
            }
        }
    }

    // (2) cross-rank agreement: every rank's query checksum must match.
    double checksum = 0.0;
    for (int k = 0; k < NK; ++k) {
        const double z = (AMREX_SPACEDIM == 3) ? coord(k, 2) : 0.0;
        amrex::ignore_unused(z);
        for (int j = 0; j < N; ++j) {
            const double y = coord(j, 1);
            for (int i = 0; i < N; ++i) checksum += shared.query(AMREX_D_DECL(coord(i, 0), y, z));
        }
    }
    double cmin = checksum, cmax = checksum;
    amrex::ParallelDescriptor::ReduceRealMin(cmin);
    amrex::ParallelDescriptor::ReduceRealMax(cmax);

    // (3) memory: per-node-lead allocation and world sum (peers contribute 0).
    double alloc_sum = static_cast<double>(shared.m_alloc_bytes);
    double alloc_max = static_cast<double>(shared.m_alloc_bytes);
    amrex::ParallelDescriptor::ReduceRealSum(alloc_sum);
    amrex::ParallelDescriptor::ReduceRealMax(alloc_max);

    amrex::Print() << "  using_window : " << shared.m_using_window
                   << ", node_size : " << shared.m_node_size << "\n"
                   << "  triangles    : " << shared.m_n_tris << ", nodes : " << shared.m_n_nodes
                   << "\n"
                   << "  max |err| vs FlatTriMeshSDF (IO rank) : " << max_err << "\n"
                   << "  cross-rank checksum spread            : " << (cmax - cmin) << "\n"
                   << "  alloc bytes (per node-lead)           : " << alloc_max << "\n"
                   << "  alloc bytes (world sum)               : " << alloc_sum
                   << "  [= per-lead * num_nodes; peers = 0]\n";

    shared.free_shared();  // collective free at a collective point
}
#endif  // AMREX_DEBUG

// ===========================================================================
// ReadEBGeometrySTL_TriMesh_NodeShared (Tier 1 step 3: user-facing Lua handle)
// ===========================================================================

ReadEBGeometrySTL_TriMesh_NodeShared::ReadEBGeometrySTL_TriMesh_NodeShared() {}

ReadEBGeometrySTL_TriMesh_NodeShared::ReadEBGeometrySTL_TriMesh_NodeShared(
  const std::string& stl_file) :
    m_filename(stl_file)
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh_NodeShared::ctor");
    m_res = NodeSharedTriMeshSDF::get_or_create(stl_file);  // COLLECTIVE on first use
}

ReadEBGeometrySTL_TriMesh_NodeShared::ReadEBGeometrySTL_TriMesh_NodeShared(
  const std::string& stl_file, bool flip_sign) :
    m_filename(stl_file), m_flip_sign(flip_sign)
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh_NodeShared::ctor");
    m_res = NodeSharedTriMeshSDF::get_or_create(stl_file);  // COLLECTIVE on first use
}

Real ReadEBGeometrySTL_TriMesh_NodeShared::query(AMREX_D_DECL(Real x, Real y, Real z)) const
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh_NodeShared::query");

    if (!m_res) {
        amrex::Abort("ReadEBGeometrySTL_TriMesh_NodeShared::query before a file was loaded");
    }
    // Backend returns the canonical signed distance; flip is applied per handle.
    const Real d = m_res->query(AMREX_D_DECL(x, y, z));
    return m_flip_sign ? -d : d;
}

const std::string ReadEBGeometrySTL_TriMesh_NodeShared::str() const
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh_NodeShared::str");

    std::stringstream ss;
    ss << "ReadEBGeometrySTL_TriMesh_NodeShared\n";
    ss << "  filename  : " << m_filename << "\n";
    ss << "  has_sdf   : " << static_cast<bool>(m_res) << "\n";
    ss << "  flip_sign : " << m_flip_sign << "\n";
    return ss.str();
}

void ReadEBGeometrySTL_TriMesh_NodeShared::set_flip_sign(bool flip_sign)
{
    m_flip_sign = flip_sign;
}

bool ReadEBGeometrySTL_TriMesh_NodeShared::get_flip_sign() const { return m_flip_sign; }

void ReadEBGeometrySTL_TriMesh_NodeShared::register_with_lua(sol::state& lua)
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh_NodeShared::register_with_lua");

    lua.new_usertype<ReadEBGeometrySTL_TriMesh_NodeShared>(
      "ReadEBGeometrySTL_TriMesh_NodeShared",
      sol::constructors<ReadEBGeometrySTL_TriMesh_NodeShared(const std::string&),
                        ReadEBGeometrySTL_TriMesh_NodeShared(const std::string&, bool)>(),
      "query",
      &ReadEBGeometrySTL_TriMesh_NodeShared::query,
      "set_flip_sign",
      &ReadEBGeometrySTL_TriMesh_NodeShared::set_flip_sign,
      "get_flip_sign",
      &ReadEBGeometrySTL_TriMesh_NodeShared::get_flip_sign,
      "str",
      &ReadEBGeometrySTL_TriMesh_NodeShared::str);
}

#endif  // AMREX_SPACEDIM > 1
