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

    lua.set_function("flat_trimesh_self_test", &FlatTriMeshSDF::self_test);
}

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

#endif  // AMREX_SPACEDIM > 1
