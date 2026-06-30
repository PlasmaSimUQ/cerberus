#include "MFP_ebgeometry_stl.H"

#if AMREX_SPACEDIM > 1

    #include <AMReX.H>
    #include <AMReX_BLProfiler.H>
    #include <sstream>
    #include <string>

ReadEBGeometrySTL::ReadEBGeometrySTL() {}

ReadEBGeometrySTL::ReadEBGeometrySTL(const std::string& stl_file)
{
    BL_PROFILE("ReadEBGeometrySTL::ReadEBGeometrySTL");
    read_file(stl_file);
}

ReadEBGeometrySTL::ReadEBGeometrySTL(const std::string& stl_file, bool flip_sign) :
    m_flip_sign(flip_sign)
{
    BL_PROFILE("ReadEBGeometrySTL::ReadEBGeometrySTL");
    read_file(stl_file);
}

void ReadEBGeometrySTL::read_file(const std::string& stl_file)
{
    BL_PROFILE("ReadEBGeometrySTL::read_file");

    m_filename = stl_file;

    /*
     * This follows the AMReX_PaintEB example pattern:
     *
     *   auto mesh = EBGeometry::Parser::readIntoDCEL<T, Meta>(filename);
     *   m_sdf = std::make_shared<EBGeometry::FastCompactMeshSDF<T, Meta, BV, K>>(mesh);
     *
     * EBGeometry parses the binary/ASCII STL, builds a DCEL mesh and a compact
     * BVH-accelerated signed distance function.
     */
    auto mesh = EBGeometry::Parser::readIntoDCEL<T, Meta>(stl_file);

    m_sdf = std::make_shared<SDF>(mesh);

    if (!m_sdf) {
        amrex::Abort("ReadEBGeometrySTL failed to construct SDF from STL file '" + stl_file + "'");
    }
}

Real ReadEBGeometrySTL::query(AMREX_D_DECL(Real x, Real y, Real z)) const
{
    BL_PROFILE("ReadEBGeometrySTL::query");

    if (!m_sdf) {
        amrex::Abort("ReadEBGeometrySTL::query called before a valid STL file was loaded");
    }

    // EBGeometry is intrinsically 3D; pad the out-of-plane coordinate in 2D.
    #if AMREX_SPACEDIM == 3
    const Vec3 p(x, y, z);
    #else
    const Vec3 p(x, y, T(0));
    #endif

    const Real d = static_cast<Real>(m_sdf->value(p));

    return m_flip_sign ? -d : d;
}

const std::string ReadEBGeometrySTL::str() const
{
    BL_PROFILE("ReadEBGeometrySTL::str");

    std::stringstream ss;

    ss << "ReadEBGeometrySTL\n";
    ss << "  filename  : " << m_filename << "\n";
    ss << "  has_sdf   : " << static_cast<bool>(m_sdf) << "\n";
    ss << "  flip_sign : " << m_flip_sign << "\n";

    return ss.str();
}

void ReadEBGeometrySTL::set_flip_sign(bool flip_sign)
{
    BL_PROFILE("ReadEBGeometrySTL::set_flip_sign");
    m_flip_sign = flip_sign;
}

bool ReadEBGeometrySTL::get_flip_sign() const
{
    BL_PROFILE("ReadEBGeometrySTL::get_flip_sign");
    return m_flip_sign;
}

void ReadEBGeometrySTL::register_with_lua(sol::state& lua)
{
    BL_PROFILE("ReadEBGeometrySTL::register_with_lua");

    lua.new_usertype<ReadEBGeometrySTL>(
      "ReadEBGeometrySTL",
      sol::constructors<ReadEBGeometrySTL(const std::string&),
                        ReadEBGeometrySTL(const std::string&, bool)>(),
      "query",
      &ReadEBGeometrySTL::query,
      "set_flip_sign",
      &ReadEBGeometrySTL::set_flip_sign,
      "get_flip_sign",
      &ReadEBGeometrySTL::get_flip_sign,
      "str",
      &ReadEBGeometrySTL::str);
}

//=============================================================================
// ReadEBGeometrySTL_TriMesh (Tier 0: FastTriMeshSDF + float precision)
//=============================================================================

ReadEBGeometrySTL_TriMesh::ReadEBGeometrySTL_TriMesh() {}

ReadEBGeometrySTL_TriMesh::ReadEBGeometrySTL_TriMesh(const std::string& stl_file)
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh::ReadEBGeometrySTL_TriMesh");
    read_file(stl_file);
}

ReadEBGeometrySTL_TriMesh::ReadEBGeometrySTL_TriMesh(const std::string& stl_file, bool flip_sign) :
    m_flip_sign(flip_sign)
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh::ReadEBGeometrySTL_TriMesh");
    read_file(stl_file);
}

void ReadEBGeometrySTL_TriMesh::read_file(const std::string& stl_file)
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh::read_file");

    m_filename = stl_file;

    // Read the STL into a DCEL mesh, then build a FastTriMeshSDF. The triangle
    // SDF copies the geometry it needs (vertex positions + baked vertex/edge
    // pseudonormals) into self-contained triangles and does NOT retain the DCEL
    // connectivity: the 'mesh' local is released when it leaves scope here.
    auto mesh = EBGeometry::Parser::readIntoDCEL<T, Meta>(stl_file);

    m_sdf = std::make_shared<SDF>(mesh);

    if (!m_sdf) {
        amrex::Abort("ReadEBGeometrySTL_TriMesh failed to construct SDF from STL file '" +
                     stl_file + "'");
    }
}

Real ReadEBGeometrySTL_TriMesh::query(AMREX_D_DECL(Real x, Real y, Real z)) const
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh::query");

    if (!m_sdf) {
        amrex::Abort("ReadEBGeometrySTL_TriMesh::query called before a valid STL file was loaded");
    }

    // EBGeometry is intrinsically 3D; pad the out-of-plane coordinate in 2D.
    // Cast explicitly to the (float) geometry precision T. Build the point from
    // named locals rather than Vec3 p(T(x), T(y), T(z)): that form is parsed as
    // a function declaration (most-vexing-parse).
    const T px = static_cast<T>(x);
    const T py = static_cast<T>(y);
    #if AMREX_SPACEDIM == 3
    const T pz = static_cast<T>(z);
    #else
    const T pz = static_cast<T>(0);
    #endif

    const Vec3 p(px, py, pz);

    const Real d = static_cast<Real>(m_sdf->value(p));

    return m_flip_sign ? -d : d;
}

const std::string ReadEBGeometrySTL_TriMesh::str() const
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh::str");

    std::stringstream ss;

    ss << "ReadEBGeometrySTL_TriMesh\n";
    ss << "  filename  : " << m_filename << "\n";
    ss << "  has_sdf   : " << static_cast<bool>(m_sdf) << "\n";
    ss << "  flip_sign : " << m_flip_sign << "\n";

    return ss.str();
}

void ReadEBGeometrySTL_TriMesh::set_flip_sign(bool flip_sign)
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh::set_flip_sign");
    m_flip_sign = flip_sign;
}

bool ReadEBGeometrySTL_TriMesh::get_flip_sign() const
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh::get_flip_sign");
    return m_flip_sign;
}

void ReadEBGeometrySTL_TriMesh::register_with_lua(sol::state& lua)
{
    BL_PROFILE("ReadEBGeometrySTL_TriMesh::register_with_lua");

    lua.new_usertype<ReadEBGeometrySTL_TriMesh>(
      "ReadEBGeometrySTL_TriMesh",
      sol::constructors<ReadEBGeometrySTL_TriMesh(const std::string&),
                        ReadEBGeometrySTL_TriMesh(const std::string&, bool)>(),
      "query",
      &ReadEBGeometrySTL_TriMesh::query,
      "set_flip_sign",
      &ReadEBGeometrySTL_TriMesh::set_flip_sign,
      "get_flip_sign",
      &ReadEBGeometrySTL_TriMesh::get_flip_sign,
      "str",
      &ReadEBGeometrySTL_TriMesh::str);
}

#endif  // AMREX_SPACEDIM > 1
