#ifndef POLYSPLINE_H
#define POLYSPLINE_H
#ifdef AMREX_USE_EB

    #include "sol.hpp"

    #include <AMReX_REAL.H>
    #include <AMReX_RealVect.H>
    #include <list>

using namespace amrex;

class Bezier
{
  public:
    RealVect p0, p1, p2, p3;
    Bezier();
    Bezier(const RealVect& pt0, const RealVect& pt1, const RealVect& pt2, const RealVect& pt3);
    void generate_coeffs();
    Vector<RealVect> sample(Real min_length, const int N, bool include_start) const;
    RealVect at(Real t) const;

  private:
    Array<RealVect, 4> coeffs;
};

class PolySpline
{
  public:
    PolySpline();

    std::vector<RealVect> grab_points(const sol::table& points);

    void add_spline_element_lua(const sol::table& pts, const Real dx);
    void add_line_element_lua(const sol::table& pts);
    void add_spline_element(const Vector<Array<RealVect, 4>>& points, const Real dx);
    void add_polyspline(const PolySpline& poly);
    void add_shape(const Vector<RealVect> shape);

    Real operator()(AMREX_D_DECL(Real x, Real y = 0.0, Real z = 0.0)) const noexcept;

    static void register_with_lua(sol::state& lua);

    Vector<Vector<RealVect>> shapes;
    Vector<Vector<RealVect>> vectors;
    Vector<Vector<Real>> lengths;
};

#endif
#endif  // POLYSPLINE_H
