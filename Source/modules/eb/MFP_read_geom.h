#ifndef READ_GEOM_H
#define READ_GEOM_H

#ifdef AMREX_USE_EB

#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include "AABB.h"
#include "MFP_polyspline.H"

using namespace amrex;

class ReadSTL
{
public:
    ReadSTL();
    ReadSTL(const std::string &stl_file);

    void read_file(const std::string &stl_file);

    bool valid_tri(const Array<Real,9>& facet);

    void get_limits(const Array<Real,9>& facet, Vector<Real>& upper, Vector<Real>& lower);

    void grow_tree();

    Real query(AMREX_D_DECL(Real x, Real y, Real z));

    const std::string str() const;

    static void register_with_lua(sol::state& lua);

    Vector<Array<Real,9>> facets;

    aabb::Tree tree;


};


#include "rapidxml.hpp"
#include "rapidxml_print.hpp"

namespace rxml = rapidxml;

class ReadSVG
{
public:
    ReadSVG();
    ReadSVG(const std::string &path, const Real dx);

    std::map<std::string, Vector<Array<RealVect,4>>> get_elements(rxml::xml_node<> * node, const std::string& parent_name="");

    Real query(AMREX_D_DECL(Real x, Real y, Real z), std::string layer="root");

    static void register_with_lua(sol::state& lua);

private:
    std::map<std::string,PolySpline> polysplines;
};

#endif
#endif // READSTL_H
