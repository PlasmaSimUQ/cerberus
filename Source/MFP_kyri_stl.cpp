#include "MFP_kyri_stl.H"

#include "MFP_diagnostics.H"
#include "MFP_utility.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <math.h>
#include <regex>
#include <sstream>
#include <stdio.h>
#include <string>

ReadSTL::ReadSTL() {}

ReadSTL::ReadSTL(const std::string& stl_file)
{
    BL_PROFILE("ReadSTL::ReadSTL");

    read_file(stl_file);
}

void ReadSTL::read_file(const std::string& stl_file)
{
    BL_PROFILE("ReadSTL::read_file");

    std::filebuf fb;
    if (fb.open(stl_file, std::ios::in | std::ios::binary)) {
        std::istream is(&fb);

        // get the first line to check if we have a binary or ASCII file
        std::string line;
        std::getline(is, line);

        if (line.find("solid") == std::string::npos) {
            is.seekg(80, is.beg);

            uint32_t ntri;
            is.read(reinterpret_cast<char*>(&ntri), sizeof(uint32_t));

            while (is) {
                Array<Real, 9> facet;

                float n;
                for (int i = 0; i < 3; ++i) { is.read(reinterpret_cast<char*>(&n), sizeof(float)); }

                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        is.read(reinterpret_cast<char*>(&n), sizeof(float));
                        facet[i * 3 + j] = n;
                    }
                }

                uint16_t att;
                is.read(reinterpret_cast<char*>(&att), sizeof(uint16_t));

                facets.push_back(facet);
            }
        } else {
            while (is) {
                Array<Real, 9> facet;

                // get to the vertices
                for (int i = 0; i < 2; ++i)
                    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // grab the vertices
                for (int i = 0; i < 3; ++i) {
                    is >> line;  // grab the 'vertex' string
                    for (int j = 0; j < 3; ++j) { is >> facet[i * 3 + j]; }
                }

                // skip ending lines

                for (int i = 0; i < 3; ++i)
                    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                facets.push_back(facet);
            }
        }
        fb.close();
    } else {
        amrex::Abort("Unable to open file '" + stl_file + "' as stl");
    }

    // grow our AABB tree
    grow_tree();
}

bool ReadSTL::valid_tri(const Array<Real, 9>& facet)
{
    BL_PROFILE("ReadSTL::valid_tri");

    Real d;

    // side lengths a, b, c
    Real a = 0.0;
    Real b = 0.0;
    Real c = 0.0;
    for (int i = 0; i < 3; ++i) {
        d = facet[0 * 3 + i] - facet[1 * 3 + i];
        a += d * d;

        d = facet[1 * 3 + i] - facet[2 * 3 + i];
        b += d * d;

        d = facet[2 * 3 + i] - facet[0 * 3 + i];
        c += d * d;
    }

    a = std::sqrt(a);
    b = std::sqrt(b);
    c = std::sqrt(c);

    if ((a + b) <= c) { return false; }

    return true;
}

void ReadSTL::get_limits(const Array<Real, 9>& facet, Vector<Real>& lower, Vector<Real>& upper)
{
    BL_PROFILE("ReadSTL::get_limits");

    // get the upper and lower bounds of the facet
    lower = {facet[0], facet[1], facet[2]};
    for (int d1 = 1; d1 < 3; ++d1) {
        const Real* v = &facet[d1 * 3];
        for (int d2 = 0; d2 < 3; ++d2) {
            lower[d2] = std::min(lower[d2], v[d2]);
            upper[d2] = std::max(upper[d2], v[d2]);
        }
    }
}

void ReadSTL::grow_tree()
{
    BL_PROFILE("ReadSTL::grow_tree");

    // initialize the tree
    tree.init(3, 0.0, facets.size(), true);

    Vector<Real> lower_bound(3), upper_bound(3);

    for (size_t idx = 0; idx < facets.size(); ++idx) {
        Array<Real, 9>& f = facets[idx];

        if (!valid_tri(f)) continue;

        get_limits(f, lower_bound, upper_bound);

        // Insert the particle into the tree.
        tree.insertParticle(idx, lower_bound, upper_bound);
    }
}

Real ReadSTL::query(AMREX_D_DECL(Real x, Real y, Real z))
{
    BL_PROFILE("ReadSTL::query");

    Real d = 1.0;

    Vector<Real> xyz = {AMREX_D_DECL(x, y, z)};

    // first check if the query point is inside the root node of the tree
    if (tree.isInside(xyz)) {
        d = tree.signedDistance(xyz, facets);
        return d;
    }

    return d;
}

const std::string ReadSTL::str() const
{
    BL_PROFILE("ReadSTL::str");

    std::stringstream ss;
    for (const auto& f : facets) {
        ss << "facet : \n";
        for (int i = 0; i < 3; ++i) {
            ss << " ";
            for (int j = 0; j < 3; ++j) {
                ss << f[i * 3 + j];
                if (j < 2) {
                    ss << ", ";
                } else {
                    ss << "\n";
                }
            }
        }
    }

    return ss.str();
}

void ReadSTL::register_with_lua(sol::state& lua)
{
    BL_PROFILE("ReadSTL::register_with_lua");

    lua.new_usertype<ReadSTL>("ReadSTL",
                              sol::constructors<ReadSTL(const std::string&)>(),
                              "query",
                              &ReadSTL::query);
}

//=============================================================================
// EBGeometry-backed STL signed-distance reader
//=============================================================================

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
     * The old ReadSTL class manually parses binary/ASCII STL triangles,
     * stores facets, constructs a Cerberus AABB tree, and evaluates
     * tree.signedDistance(...). Here we delegate that work to EBGeometry.
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

    const Real d = static_cast<Real>(m_sdf->value(Vec3(AMREX_D_DECL(x, y, z))));

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
//====
