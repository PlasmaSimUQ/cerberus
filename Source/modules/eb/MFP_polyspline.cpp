#ifdef AMREX_USE_EB

#include "MFP_polyspline.H"

#include <math.h>

#include "MFP_utility.H"

#include <iostream>
#include <iomanip>
#include <list>


#include "MFP_utility.H"
#include "MFP_diagnostics.H"

Bezier::Bezier(){}

Bezier::Bezier(const RealVect& pt0, const RealVect& pt1, const RealVect& pt2, const RealVect& pt3)
    : p0(pt0), p3(pt3), p1(pt1), p2(pt2)
{
    generate_coeffs();
}

void Bezier::generate_coeffs()
{
    // generate the coefficients
    coeffs[0] = -p0 + 3*p1 - 3*p2 + p3;
    coeffs[1] = 3*p0 - 6*p1 + 3*p2;
    coeffs[2] = -3*p0 + 3*p1;
    coeffs[3] = p0;
}

RealVect Bezier::at(Real t) const
{
    RealVect pt;
    pt  = coeffs[3];
    pt += coeffs[2]*t;
    pt += coeffs[1]*t*t;
    pt += coeffs[0]*t*t*t;

    return pt;
}

// sample Bezier curve with minimum arc length between sample points
Vector<RealVect> Bezier::sample(Real min_length, const int N, bool include_start) const
{
    min_length *= min_length; // deal in length squared

    std::list<std::pair<Real, RealVect>> points;

    points.push_back(std::make_pair(0.0,p0));
    points.push_back(std::make_pair(1.0,p3));

    Real l2;

    auto it = points.begin();

    while (it->first < 1.0) {
        Real ta = it->first;
        RealVect& a = it->second;
        ++it;
        Real tb = it->first;
        RealVect& b = it->second;

        RealVect v = b-a;
        l2 = v.dotProduct(v);

        if (l2 > min_length) {
            Real t = (ta + tb)/2;
            points.insert(it, std::make_pair(t, at(t)));
            --it;
            --it;
        }
    }

    Vector<RealVect> pts;

    // check if we have a zero length segment
    if (points.size() == 2 && l2 < 1e-12) {
        if (include_start) {
            pts = {p0};
        }
        return pts;
    }


    pts.reserve(points.size());
    it = points.begin();
    if (!include_start) ++it;
    while (it != points.end()) {
        pts.push_back(it->second);
        ++it;
    }

    return pts;

}


PolySpline::PolySpline()
{

}

void PolySpline::add_spline_element_lua(const sol::table& points, const Real dx)
{

    Vector<RealVect> pts;
    sol::table pt;
    int n_pts = points.size();

    if (n_pts < 4) Abort("Not enough points to define a cubic spline, need at least 4");

    Array<RealVect,4> bzp;
    size_t pt_id = 1;
    size_t offset = 0;
    while (pt_id <= n_pts) {

        for (size_t i=0; i<4; ++i) {
            pt = points[pt_id - offset];
            bzp[i] = RealVect(AMREX_D_DECL(pt[1],pt[2],pt[3]));
            ++pt_id;
        }


        Bezier b(bzp[0], bzp[1], bzp[2], bzp[3]);

        Vector<RealVect> bpts = b.sample(dx, 10, offset == 0);
        for (const auto& bpt : bpts) {
            if (pts.empty()) {
                pts.push_back(bpt);
            } else if (bpt != pts.back()) {
                pts.push_back(bpt);
            }
        }

        offset = 1;

    }

    add_shape(pts);
}

void PolySpline::add_line_element_lua(const sol::table& points)
{

    Vector<amrex::RealVect> pts;
    sol::table pt;
    int n_pts = points.size();
    int pt_id = 1;
    while (pt_id <= n_pts) {
        pt = points[pt_id]; ++pt_id;
        pts.push_back(RealVect(AMREX_D_DECL(pt[1],pt[2],pt[3])));
    }


    add_shape(pts);
}

void PolySpline::add_spline_element(const Vector<Array<RealVect,4>>& points, const Real dx)
{

    Vector<RealVect> pts;

    size_t offset = 0;
    for (const auto& bzp : points) {

        Bezier b(bzp[0], bzp[1], bzp[2], bzp[3]);

        Vector<RealVect> bpts = b.sample(dx, 20, offset == 0);
        for (const auto& bpt : bpts) {
            if (pts.empty()) {
                pts.push_back(bpt);
            } else if (bpt != pts.back()) {
                pts.push_back(bpt);
            }
        }

        offset = 1;

    }

    add_shape(pts);

}

void PolySpline::add_polyspline(const PolySpline& poly)
{
    shapes.insert(shapes.end(), poly.shapes.begin(), poly.shapes.end());
}

void PolySpline::add_shape(const Vector<RealVect> shape)
{
    size_t n_vec = shape.size() - 1;
    Vector<RealVect> vec(n_vec);
    Vector<Real> len(n_vec);

    for (size_t i=0; i<n_vec; ++i) {
        const RealVect& a = shape[i];
        const RealVect& b = shape[i+1];

        vec[i] = b-a;
        len[i] = vec[i].vectorLength();
        vec[i] /= len[i];
    }

    shapes.push_back(shape);
    vectors.push_back(vec);
    lengths.push_back(len);
}

/*
 * Calculate the signed distance function for an arbitrary number of splines in the x-y plane
 */
Real PolySpline::operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
{

//    plot_poly_spline(*this, "poly", true);


    Real d2 = std::numeric_limits<Real>::max();
    int side;

    RealVect p(AMREX_D_DECL(x, y, z));

    RealVect pa, pb;
    Real t, dn;
    int nsegs;

    for (size_t shape_idx=0; shape_idx<shapes.size(); ++shape_idx) {

        const Vector<RealVect>& shape = shapes[shape_idx];
        const Vector<RealVect>& vector = vectors[shape_idx];
        const Vector<Real>& len = lengths[shape_idx];

        Real dd = std::numeric_limits<Real>::max(); // distance squared
        int wn = 0;               // winding number (inside/outside)

        // iterate over the line segments
        nsegs = shape.size()-1;

        pb = p - shape[0];

        for (size_t i = 0; i<nsegs; ++i) {

            const RealVect& b = shape[i+1];

            // vector a -> b
            const RealVect& v = vector[i];
            const Real& l = len[i];

            pa = pb;
            pb = p - b;

            t = pa.dotProduct(v); // t-parameter of projection onto line
            dn = pa[0]*v[1] - pa[1]*v[0]; // normal distance from p to line

            // Distance to line segment
            if (t < 0) {
                dd = std::min(dd, pa.dotProduct(pa)); // distance to vertex[0] of line
            } else if (t > l) {
                dd = std::min(dd, pb.dotProduct(pb)); // distance to vertex[1] of line
            } else {
                dd = std::min(dd, dn*dn); // normal distance to line
            }

            // Is the point in the polygon?
            // See: http://geomalgorithms.com/a03-_inclusion.html
            const RealVect& a = shape[i];

            if (a[1] <= y) {
                if ((b[1] > y) and (dn < 0)) {
                    wn = wn + 1;
                }
            } else {
                if ((b[1] <= y) and (dn > 0)) {
                    wn = wn - 1;
                }
            }

        }

        if (dd < d2) {
            d2 = dd;
            side = wn;
        }

    }


    // normalise d*d to d
    d2 = std::sqrt(d2);
    if (side != 0) {
        // p is inside the polygon
        d2 *= -1;
    }

    return d2;
}

void PolySpline::register_with_lua(sol::state& lua)
{
    lua.new_usertype<PolySpline>("PolySpline", sol::constructors<PolySpline()>(),
                                 "addLine", &PolySpline::add_line_element_lua,
                                 "addSpline", &PolySpline::add_spline_element_lua,
                                 "query", &PolySpline::operator());
}

#endif
