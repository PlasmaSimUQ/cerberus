#ifdef AMREX_USE_EB
#include "MFP_read_geom.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdint>
#include <math.h>
#include <stdio.h>
#include <regex>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>

#include "MFP_utility.H"
#include "MFP_diagnostics.H"


ReadSTL::ReadSTL()
{

}

ReadSTL::ReadSTL(const std::string& stl_file)
{
    read_file(stl_file);
}

void ReadSTL::read_file(const std::string &stl_file)
{
    BL_PROFILE("ReadSTL::read_file()");

    std::filebuf fb;
    if (fb.open (stl_file, std::ios::in | std::ios::binary)) {
        std::istream is(&fb);

        // get the first line to check if we have a binary or ASCII file
        std::string line;
        std::getline(is, line);

        if (line.find("solid") == std::string::npos) {

            is.seekg(80, is.beg);

            uint32_t ntri;
            is.read(reinterpret_cast<char *>(&ntri), sizeof(uint32_t));

            while (is) {
                Array<Real,9> facet;

                float n;
                for (int i=0; i<3; ++i) {
                    is.read(reinterpret_cast<char *>(&n), sizeof(float));
                }

                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        is.read(reinterpret_cast<char *>(&n), sizeof(float));
                        facet[i*3+j] = n;
                    }
                }

                uint16_t att;
                is.read(reinterpret_cast<char *>(&att), sizeof(uint16_t));

                facets.push_back(facet);

            }
        } else {

            while (is) {
                Array<Real,9> facet;

                // get to the vertices
                for (int i=0; i<2; ++i)
                    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // grab the vertices
                for (int i=0; i<3; ++i) {
                    is >> line; // grab the 'vertex' string
                    for (int j=0; j<3; ++j) {
                        is >> facet[i*3+j];
                    }
                }

                // skip ending lines

                for (int i=0; i<3; ++i)
                    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                facets.push_back(facet);
            }
        }
        fb.close();
    } else {
        amrex::Abort("Unable to open file '"+stl_file+"' as stl");
    }

    // grow our AABB tree
    grow_tree();
}

bool ReadSTL::valid_tri(const Array<Real,9>& facet)
{

    BL_PROFILE("ReadSTL::valid_tri()");

    Real d;

    // side lengths a, b, c
    Real a = 0.0;
    Real b = 0.0;
    Real c = 0.0;
    for (int i=0; i<3; ++i) {
        d = facet[0*3 + i] - facet[1*3 + i];
        a += d*d;

        d = facet[1*3 + i] - facet[2*3 + i];
        b += d*d;

        d = facet[2*3 + i] - facet[0*3 + i];
        c += d*d;
    }

    a = std::sqrt(a);
    b = std::sqrt(b);
    c = std::sqrt(c);

    if ((a + b) <= c) {
        return false;
    }

    return true;
}

void ReadSTL::get_limits(const Array<Real,9>& facet, Vector<Real>& lower, Vector<Real>& upper)
{
    BL_PROFILE("ReadSTL::get_limits()");

    // get the upper and lower bounds of the facet
    lower = {facet[0], facet[1], facet[2]};
    for (int d1=1; d1<3; ++d1) {
        const Real* v = &facet[d1*3];
        for (int d2=0; d2<3; ++d2) {
            lower[d2] = std::min(lower[d2], v[d2]);
            upper[d2] = std::max(upper[d2], v[d2]);


        }
    }
}

void ReadSTL::grow_tree()
{
    BL_PROFILE("ReadSTL::grow_tree()");

    // initialize the tree
    tree.init(3, 0.0, facets.size(), true);

    Vector<Real> lower_bound(3), upper_bound(3);

    for (size_t idx = 0; idx < facets.size(); ++idx) {

        Array<Real,9>& f = facets[idx];

        if (!valid_tri(f)) continue;

        get_limits(f, lower_bound, upper_bound);

        // Insert the particle into the tree.
        tree.insertParticle(idx, lower_bound, upper_bound);

    }
}

Real ReadSTL::query(AMREX_D_DECL(Real x, Real y, Real z))
{
    BL_PROFILE("ReadSTL::query()");

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
    std::stringstream ss;
    for (const auto& f : facets) {

        ss << "facet : \n";
        for (int i=0; i<3; ++i) {
            ss << " ";
            for (int j=0; j<3; ++j) {
                ss << f[i*3+j];
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
    lua.new_usertype<ReadSTL>("ReadSTL", sol::constructors<ReadSTL(const std::string&)>(),
                              "query", &ReadSTL::query);
}


//=============================================================================

std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    std::stringstream ss(s+' ');
    std::string item;
    while(std::getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}

std::vector<double> get_numbers(std::string s) {
    // get all of the numbers that are present in the string
    std::smatch m;
    std::regex re("[-+]?\\d*\\.?\\d+[eE]?[-+]?\\d*");
    std::vector<double> numbers;

    while (std::regex_search (s,m,re)) {
        for (auto x:m) {
            numbers.push_back(std::stod(x));
        }
        s = m.suffix().str();
    }
    return numbers;
}

double get_number(std::string s, const double fallback=0.0) {
    // get all of the numbers that are present in the string
    std::smatch m;
    std::regex re("[-+]?\\d*\\.?\\d+[eE]?[-+]?\\d*");
    std::vector<double> numbers;

    std::regex_search (s,m,re);
    if (!m.empty()) return std::stod(m[0]);
    return fallback;
}

std::vector<double> get_n_numbers(std::vector<std::string> words, int n, size_t* i) {
    std::vector<double> numbers;

    while (numbers.size() < n) {
        *i += 1;
        std::string w = words[*i];
        std::vector<double> num = get_numbers(w);
        if (num.empty()) Abort("Not enough numbers!");
        numbers.insert(numbers.end(), num.begin(), num.end());
    }

    return numbers;
}

void insert_segment(Vector<Array<RealVect,4>>& curve, const Array<RealVect,4>& seg)
{

    if (!curve.empty()) {
        Array<RealVect,4>& prev = curve.back();

        // check that the previous segment isn't a point
        bool has_difference = false;
        for (size_t i = 0; i<3; ++i) {
            if (prev[i] != prev[i+1]) {
                has_difference = true;
            }
        }

        // if it is a point, remove it!
        if (!has_difference) {
            curve.pop_back();
        }
    }

    // if the curve is now empty, just add the segment
    if (curve.empty()) {
        curve.push_back(seg);
        return;
    }

    // otherwise, check that we aren't simply adding a duplicate
    Array<RealVect,4>& prev = curve.back();
    bool has_difference = false;
    for (size_t i = 0; i<4; ++i) {
        if (prev[i] != seg[i]) {
            has_difference = true;
            break;
        }
    }

    if (has_difference) {
        curve.push_back(seg);
        return;
    }
}

Vector<Array<RealVect,4>> get_bezier_circle(const double& cx, const double& cy, const double& r)
{
    // magic number for roundness, https://spencermortensen.com/articles/bezier-circle/
    const double c = r*0.551915024494;

    Vector<Array<RealVect,4>> circle(4);
    RealVect p0, p1, p2, p3;

    // side 1
    p0 = RealVect(AMREX_D_DECL(cx,cy+r, 0.0));
    p1 = RealVect(AMREX_D_DECL(cx+c,cy+r, 0.0));
    p2 = RealVect(AMREX_D_DECL(cx+r,cy+c, 0.0));
    p3 = RealVect(AMREX_D_DECL(cx+r,cy, 0.0));
    circle[0] = {p0, p1, p2, p3};

    // side 2
    p0 = RealVect(AMREX_D_DECL(cx+r,cy, 0.0));
    p1 = RealVect(AMREX_D_DECL(cx+r,cy-c, 0.0));
    p2 = RealVect(AMREX_D_DECL(cx+c,cy-r, 0.0));
    p3 = RealVect(AMREX_D_DECL(cx,cy-r, 0.0));
    circle[1] = {p0, p1, p2, p3};

    // side 3
    p0 = RealVect(AMREX_D_DECL(cx,cy-r, 0.0));
    p1 = RealVect(AMREX_D_DECL(cx-c,cy-r, 0.0));
    p2 = RealVect(AMREX_D_DECL(cx-r,cy-c, 0.0));
    p3 = RealVect(AMREX_D_DECL(cx-r,cy, 0.0));
    circle[2] = {p0, p1, p2, p3};

    // side 4
    p0 = RealVect(AMREX_D_DECL(cx-r,cy, 0.0));
    p1 = RealVect(AMREX_D_DECL(cx-r,cy+c, 0.0));
    p2 = RealVect(AMREX_D_DECL(cx-c,cy+r, 0.0));
    p3 = RealVect(AMREX_D_DECL(cx,cy+r, 0.0));
    circle[3] = {p0, p1, p2, p3};

    return circle;
}

float svg_angle( Real ux, Real uy, Real vx, Real vy )
{
    RealVect u(AMREX_D_DECL(ux, uy, 0.0));
    RealVect v(AMREX_D_DECL(vx, vy, 0.0));
    Real dot = u.dotProduct(v);
    Real len = u.vectorLength() * v.vectorLength();
    Real ang = std::acos( std::clamp(dot / len, -1.0, 1.0) ); //floating point precision, slightly over values appear
    if ( (u[0]*v[1] - u[1]*v[0]) < 0) ang = -ang;
    return ang;
}

RealVect elliptic_arc_point( RealVect c, RealVect r, float phi, float t )
{
     Real a = std::sin(phi);
     Real b = std::cos(phi);
    return RealVect(AMREX_D_DECL(
            c[0] + r[0] * b * std::cos(t) - r[1] * a * std::sin(t),
            c[1] + r[0] * a * std::cos(t) + r[1] * b * std::sin(t),
            0.0));
}

 RealVect elliptic_arc_derivative(RealVect r, float phi, float t )
 {
     Real a = std::sin(phi);
     Real b = std::cos(phi);
     Real st = std::sin(t);
     Real ct = std::cos(t);
     return RealVect(AMREX_D_DECL(
            -r[0] * b * st - r[1] * a * ct,
            -r[0] * a * st + r[1] * b * ct,
            0.0));
 }

Vector<Array<RealVect,4>> get_bezier_elliptical_arc(Real x1, Real y1,
                                                    Real rx, Real ry, Real phi,
                                                    Real fA, Real fS,
                                                    Real x2, Real y2)
{

    // calculate the centre of the ellipse
    // https://mortoray.com/2017/02/16/rendering-an-svg-elliptical-arc-as-bezier-curves/
    static constexpr Real pi = acos(-1);

    Real rX = std::abs(rx);
    Real rY = std::abs(ry);

    //(F.6.5.1)
    Real dx2 = (x1 - x2) / 2.0;
    Real dy2 = (y1 - y2) / 2.0;
    Real x1p = std::cos(phi)*dx2 + std::sin(phi)*dy2;
    Real y1p = -std::sin(phi)*dx2 + std::cos(phi)*dy2;

    //(F.6.5.2)
    Real rxs = rX * rX;
    Real rys = rY * rY;
    Real x1ps = x1p * x1p;
    Real y1ps = y1p * y1p;
    // check if the radius is too small `pq < 0`, when `dq > rxs * rys` (see below)
    // cr is the ratio (dq : rxs * rys)
    Real cr = x1ps/rxs + y1ps/rys;
    if (cr > 1) {
        //scale up rX,rY equally so cr == 1
        Real s = std::sqrt(cr);
        rX = s * rX;
        rY = s * rY;
        rxs = rX * rX;
        rys = rY * rY;
    }
    Real dq = (rxs * y1ps + rys * x1ps);
    Real pq = (rxs*rys - dq) / dq;
    Real q = std::sqrt( std::max(0.0,pq) ); //use Max to account for float precision
    if (fA == fS)
        q = -q;
    Real cxp = q * rX * y1p / rY;
    Real cyp = - q * rY * x1p / rX;

    //(F.6.5.3)
    Real cx = std::cos(phi)*cxp - std::sin(phi)*cyp + (x1 + x2)/2;
    Real cy = std::sin(phi)*cxp + std::cos(phi)*cyp + (y1 + y2)/2;


    //(F.6.5.5)
    Real start_angle = svg_angle( 1,0, (x1p-cxp) / rX, (y1p - cyp)/rY );
    //(F.6.5.6)
    Real delta = svg_angle((x1p - cxp)/rX, (y1p - cyp)/rY,(-x1p - cxp)/rX, (-y1p-cyp)/rY);
    delta = std::fmod(delta, pi * 2.0);
    if (!fS) delta -= 2 * pi;
//    int sign = (end_angle < start_angle) ? -1 : 1;

    // we now have the centre of the ellipse (cx, cy) as well as the angles of each of the end points
    RealVect c(AMREX_D_DECL(cx, cy, 0.0));
    RealVect r(AMREX_D_DECL(rX, rY, 0.0));


    RealVect p0, p1, p2, p3;

    p3 = RealVect(AMREX_D_DECL(x1, y1, 0.0)); // start of this curve, end of last

    const Real n_div = 4.0;

    Real angle_step = delta/n_div;
    Real alphaT = std::tan(angle_step/2.0);
    alphaT *= alphaT;

    Vector<Array<RealVect,4>> ellipse(4);
    for (size_t i=0; i<n_div; ++i) {
        Real angle0 = start_angle + i*angle_step;
        Real angle3 = start_angle + (i+1)*angle_step;

        p0 = p3; // the starting point
        p3 = elliptic_arc_point(c,r,phi,angle3); // the ending point


        Real alpha = std::sin(angle_step) * (std::sqrt(4 + 3 * alphaT)- 1) / 3.0;

        p1 = p0 + alpha*elliptic_arc_derivative(r, phi, angle0);
        p2 = p3 - alpha*elliptic_arc_derivative(r, phi, angle3);

        ellipse[i] = {p0, p1, p2, p3};

    }

    return ellipse;
}


ReadSVG::ReadSVG()
{

}

ReadSVG::ReadSVG(const std::string &path, const Real dx)
{

    // Read the file
    rxml::xml_document<> doc;
    std::ifstream theFile (path);
    std::vector<char> buffer((std::istreambuf_iterator<char>(theFile)), std::istreambuf_iterator<char>());
    buffer.push_back('\0');
    // Parse the buffer using the xml file parsing library into doc
    doc.parse<0>(&buffer[0]);

    // Find our root node
    rxml::xml_node<> * root_node;
    root_node = doc.first_node("svg");

    // get the width and height
    double width = get_number(std::string(root_node->first_attribute("width")->value()));
    double height = get_number(std::string(root_node->first_attribute("height")->value()));
    std::vector<double> viewBox = get_numbers(std::string(root_node->first_attribute("viewBox")->value()));

    RealVect doc_scale;
    doc_scale[0] = width/(viewBox[2] - viewBox[0]);
    doc_scale[1] = height/(viewBox[3] - viewBox[1]);

    std::map<std::string, Vector<Array<RealVect,4>>> paths = get_elements(root_node);

    RealVect doc_origin(0.0);

    // check for a defined origin and manipulate all of the paths to suit
    // this means that "origin" is a priveleged keyword
    for (auto& shape : paths) {
        if (shape.first.find("origin") != std::string::npos) {
            doc_origin = shape.second[0][0];
        }
    }

    // perform geometric translation to expected coordinate system
    for (auto& shape : paths) {
        if (shape.first.find("origin") != std::string::npos) continue;
        for (Array<RealVect,4>& nodes : shape.second) {
            for (RealVect& node : nodes) {
                node -= doc_origin;
                node[1] *= -1;
                node *= doc_scale;
            }
        }

        std::string path = shape.first;
        std::size_t pos = path.find_last_of("/");
        std::string layer = path.substr(0, pos);

        polysplines[layer].add_spline_element(shape.second, dx);
    }

    // do some grouping within the hierarchy

    Vector<std::string> groups;
    for (auto& poly : polysplines) {
        std::string group_name = poly.first;
        size_t n = std::count(group_name.begin(), group_name.end(), '/');
        while (n > 1) {
            n = group_name.find_last_of('/');
            group_name = group_name.substr(0,n);
            groups.push_back(group_name);
            n = std::count(group_name.begin(), group_name.end(), '/');
        }
    }

    for (const auto& group_name : groups) {
        PolySpline& poly_group = polysplines[group_name];

        for (auto& poly : polysplines) {
            if (poly.first.find(group_name) != std::string::npos) {
                poly_group.add_polyspline(poly.second);
            }
        }
    }

//    for (const auto& poly : polysplines) {
//        plot_poly_spline(poly.second, poly.first, true);
//    }

}


// useful reference for SVG parsing
// https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute

std::map<std::string, Vector<Array<RealVect,4>>> ReadSVG::get_elements(rxml::xml_node<> * node, const std::string& parent_name)
{

    std::map<std::string, Vector<Array<RealVect,4>>> curves;

    // we're probably going to use inkscape so test for inkscape labels before falling back to the generic id
    std::string id;
    rxml::xml_attribute<>* get_id = node->first_attribute("inkscape:label");
    if (get_id) {
        id = get_id->value();
    } else {
        id = node->first_attribute("id")->value();
    }

    std::string tag = node->name();
    std::string name = parent_name + "/" + id;

    if ((tag == "g") || (tag == "svg")) {
        std::string level_name;
        if (tag == "g") {
            level_name = name;
        } else {
            level_name = "";
        }
        for (rxml::xml_node<>* sub_node = node->first_node(); sub_node != NULL; sub_node = sub_node->next_sibling()) {
            std::map<std::string, Vector<Array<RealVect,4>>> sub_curves = get_elements(sub_node, level_name);
            curves.insert(sub_curves.begin(), sub_curves.end());
        }
    } else if (tag == "path") {
        std::string d = node->first_attribute("d")->value();

        std::vector<std::string> words = split(d, ' ');

        Vector<Array<RealVect,4>> curve;
        RealVect p0, p1, p2, p3;

        size_t i = 0;

        std::string last_command;
        std::vector<double> num;

        // parse the svg path commands
        while (i < words.size()) {
            std::string w = words[i];



            // check if it is a command
            if (w.size() > 1) {
                w = last_command;
                i -= 1;
            }

            if (w == "m") {
                last_command = w;
                num = get_n_numbers(words, 2, &i);
                if (curve.empty()) {
                    p0 = RealVect(AMREX_D_DECL(num[0],num[1], 0.0));
                    insert_segment(curve, {p0, p0, p0, p0});
                } else {
                    p0 = curve.back().at(3);
                    p3 = p0 +  RealVect(AMREX_D_DECL(num[0],num[1], 0.0));
                    insert_segment(curve, {p0, p0, p3, p3});
                }
            } else if (w == "M") {
                last_command = w;
                num = get_n_numbers(words, 2, &i);
                if (curve.empty()) {
                    p0 = RealVect(AMREX_D_DECL(num[0],num[1], 0.0));
                    insert_segment(curve, {p0, p0, p0, p0});
                } else {
                    p0 = curve.back().at(3);
                    p3 = RealVect(AMREX_D_DECL(num[0],num[1], 0.0));
                    insert_segment(curve, {p0, p0, p3, p3});
                }
            } else if (w == "v") {
                last_command = w;
                num = get_n_numbers(words, 1, &i);
                p0 = curve.back().at(3);
                p3 = p0 + RealVect(AMREX_D_DECL(0,num[0], 0.0));
                insert_segment(curve, {p0, p0, p3, p3});
            } else if (w == "V") {
                last_command = w;
                num = get_n_numbers(words, 1, &i);
                p3 = RealVect(AMREX_D_DECL(p0[0],num[0], 0.0));
                insert_segment(curve, {p0, p0, p3, p3});
            } else if (w == "h") {
                last_command = w;
                num = get_n_numbers(words, 1, &i);
                p0 = curve.back().at(3);
                p3 = p0 + RealVect(AMREX_D_DECL(num[0], 0.0, 0.0));
                insert_segment(curve, {p0, p0, p3, p3});
            } else if (w == "H") {
                last_command = w;
                num = get_n_numbers(words, 1, &i);
                p0 = curve.back().at(3);
                p3 = RealVect(AMREX_D_DECL(num[0], p0[1], 0));
                insert_segment(curve, {p0, p0, p3, p3});
            } else if (w == "l") {
                last_command = w;
                num = get_n_numbers(words, 2, &i);
                p0 = curve.back().at(3);
                p3 = p0 + RealVect(AMREX_D_DECL(num[0], num[1], 0.0));
                insert_segment(curve, {p0, p0, p3, p3});
            } else if (w == "L") {
                last_command = w;
                num = get_n_numbers(words, 2, &i);
                p0 = curve.back().at(3);
                p3 = RealVect(AMREX_D_DECL(num[0], num[1], 0.0));
                insert_segment(curve, {p0, p0, p3, p3});
            } else if (w == "c") {
                last_command = w;
                num = get_n_numbers(words, 6, &i);
                p0 = curve.back().at(3);
                p1 = p0 + RealVect(AMREX_D_DECL(num[0], num[1], 0.0));
                p2 = p0 + RealVect(AMREX_D_DECL(num[2], num[3], 0.0));
                p3 = p0 + RealVect(AMREX_D_DECL(num[4], num[5], 0.0));
                insert_segment(curve, {p0, p1, p2, p3});
            } else if (w == "C") {
                last_command = w;
                num = get_n_numbers(words, 6, &i);
                p0 = curve.back().at(3);
                p1 = RealVect(AMREX_D_DECL(num[0], num[1], 0.0));
                p2 = RealVect(AMREX_D_DECL(num[2], num[3], 0.0));
                p3 = RealVect(AMREX_D_DECL(num[4], num[5], 0.0));
                insert_segment(curve, {p0, p1, p2, p3});
            } else if (w == "s") {
                last_command = w;
                num = get_n_numbers(words, 4, &i);
                p0 = curve.back().at(3);
                p1 = curve.back().at(2); // previous end control point
                p1 = 2*p0 - p1; // continue smoothly on
                p2 = p0 + RealVect(AMREX_D_DECL(num[0], num[1], 0.0));
                p3 = p0 + RealVect(AMREX_D_DECL(num[2], num[3], 0.0));
                insert_segment(curve, {p0, p1, p2, p3});
            } else if (w == "S") {
                last_command = w;
                num = get_n_numbers(words, 4, &i);
                p0 = curve.back().at(3);
                p1 = curve.back().at(2); // previous end control point
                p1 = 2*p0 - p1; // continue smoothly on
                p2 = RealVect(AMREX_D_DECL(num[0], num[1], 0.0));
                p3 = RealVect(AMREX_D_DECL(num[2], num[3], 0.0));
                insert_segment(curve, {p0, p1, p2, p3});
                insert_segment(curve, {p0, p1, p2, p3});
            } else if (w == "a") {
                last_command = w;
                num = get_n_numbers(words, 7, &i);
                p0 = curve.back().at(3);
                p3 = p0 +  RealVect(AMREX_D_DECL(num[5],num[6], 0.0));
                for (const auto& seg : get_bezier_elliptical_arc(p0[0], p0[1], num[0], num[1], num[2], num[3], num[4], p3[0], p3[1])) {
                    insert_segment(curve, seg);
                }
            } else if (w == "A") {
                last_command = w;
                num = get_n_numbers(words, 7, &i);
                p0 = curve.back().at(3);
                p3 = RealVect(AMREX_D_DECL(num[5],num[6], 0.0));
                for (const auto& seg : get_bezier_elliptical_arc(p0[0], p0[1], num[0], num[1], num[2], num[3], num[4], p3[0], p3[1])) {
                    insert_segment(curve, seg);
                }
            } else if (w == "z" || w == "Z") {
                p0 = curve.back().at(3);
                p3 = curve[0].at(0);
                insert_segment(curve, {p0, p0, p3, p3});
            }
            i += 1;
            //            std::cout << curve.back()[0] << ", " << curve.back()[1] << ", " << curve.back()[2] << ", "<< curve.back()[3] << std::endl;
        }

        curves[name] = curve;
    } else if (tag == "rect") {
        double x_lo = std::stod(node->first_attribute("x")->value());
        double y_lo = std::stod(node->first_attribute("y")->value());
        double height = std::stod(node->first_attribute("height")->value());
        double width = std::stod(node->first_attribute("width")->value());

        double x_hi = x_lo + width;
        double y_hi = y_lo + height;

        Vector<Array<RealVect,4>> rect(4);
        RealVect p0, p3;

        // side 1
        p0 = RealVect(AMREX_D_DECL(x_lo,y_lo, 0.0));
        p3 = RealVect(AMREX_D_DECL(x_lo,y_hi, 0.0));
        rect[0] = {p0, p0, p3, p3};

        // side 2
        p0 = RealVect(AMREX_D_DECL(x_lo,y_hi, 0.0));
        p3 = RealVect(AMREX_D_DECL(x_hi,y_hi, 0.0));
        rect[1] = {p0, p0, p3, p3};

        // side 3
        p0 = RealVect(AMREX_D_DECL(x_hi,y_hi, 0.0));
        p3 = RealVect(AMREX_D_DECL(x_hi,y_lo, 0.0));
        rect[2] = {p0, p0, p3, p3};

        // side 4
        p0 = RealVect(AMREX_D_DECL(x_hi,y_lo, 0.0));
        p3 = RealVect(AMREX_D_DECL(x_lo,y_lo, 0.0));
        rect[3] = {p0, p0, p3, p3};

        curves[name] = rect;

    } else if (tag == "circle") {
        // approximate a circle with a collection of 4 cubic Bezier curves
        double cx = std::stod(node->first_attribute("cx")->value());
        double cy = std::stod(node->first_attribute("cy")->value());
        double r = std::stod(node->first_attribute("r")->value());

        if (id == "origin") {
            RealVect p0(AMREX_D_DECL(cx,cy, 0.0));
            curves[name].push_back({p0, p0, p0, p0});
        } else {

            curves[name] = get_bezier_circle(cx, cy, r);
        }

    }


    // check if there is a transform to apply
    rxml::xml_attribute<> *attr = node->first_attribute("transform");
    if (attr) {
        std::string transform = attr->value();

        double a=1, b=0, c=0, d=1, e=0, f=0;

        //            In the following we are applying the matrix:
        //            [a c e]
        //            [b d f]
        //            [0 0 1]
        //            which transform the points as such:
        //            newX = a * oldX + c * oldY + e
        //            newY = b * oldX + d * oldY + f

        // get all of the numbers that are present in the string

        std::vector<double> numbers = get_numbers(transform);

        // now go through all the options for transformations
        size_t found;

        // rotation
        found = transform.find("rotate");
        if (found != std::string::npos) {

            a = cos(numbers[0]);
            b = sin(numbers[0]);
            c = -b;
            d = a;

            if (numbers.size() == 3) {
                e = numbers[1] - a*numbers[1] - c*numbers[2];
                f = numbers[2] - b*numbers[1] - d*numbers[2];
            }
        }

        // translation
        found = transform.find("translate");
        if (found != std::string::npos) {
            if (numbers.size() == 2) {
                e = numbers[0];
                f = numbers[1];
            }
        }


        // scale
        found = transform.find("scale");
        if (found != std::string::npos) {
            if (numbers.size() >= 1) {
                a = numbers[0];
                d = numbers[0];
            }
            if (numbers.size() == 2) {
                d = numbers[1];
            }
        }

        // skewX
        found = transform.find("skewX");
        if (found != std::string::npos) {
            c = sin(numbers[0]);
        }

        // skewY
        found = transform.find("skewY");
        if (found != std::string::npos) {
            b = sin(numbers[0]);
        }

        // matrix
        found = transform.find("matrix");
        if (found != std::string::npos) {
            a = numbers[0];
            b = numbers[1];
            c = numbers[2];
            d = numbers[3];
            e = numbers[4];
            f = numbers[5];
        }

        // apply transformation to all shapes
        for (auto& shape : curves) {
            for (Array<RealVect,4>& nodes : shape.second) {
                for (RealVect& node : nodes) {
                    Real x = node[0];
                    Real y = node[1];
                    node[0] = a*x + c*y + e;
                    node[1] = b*x + d*y + f;
                }
            }
        }

    }

    return curves;
}

Real ReadSVG::query(AMREX_D_DECL(Real x, Real y, Real z), std::string layer)
{
    // check layer exists
    if (polysplines.find(layer) == polysplines.end()) {
        std::vector<std::string> names;
        for (auto& ps : polysplines) {
            names.push_back(ps.first);
        }
        std::stringstream msg;
        msg << "ReadSVG requested layer '" << layer << "' which does not exist. Try :" << vec2str(names);
        Abort(msg.str());
    }

    return polysplines[layer](AMREX_D_DECL(x, y, z));
}

void ReadSVG::register_with_lua(sol::state& lua)
{
    lua.new_usertype<ReadSVG>("ReadSVG", sol::constructors<ReadSVG(const std::string&, const Real)>(),
                              "query", &ReadSVG::query);
}

#endif
