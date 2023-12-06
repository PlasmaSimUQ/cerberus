#ifndef MFP_EB_H
#define MFP_EB_H
#ifdef AMREX_USE_EB

    #include "sol.hpp"

    #include <AMReX_Array.H>
    #include <AMReX_EB2_IF_Base.H>
    #include <AMReX_EB2_IF_Spline.H>
    #include <AMReX_REAL.H>
    #include <algorithm>
    #include <limits>

using namespace amrex;

// For all implicit functions, >0: body; =0: boundary; <0: fluid

typedef std::function<Real(AMREX_D_DECL(Real, Real, Real))> lua_IF;

class FunctionIF : GPUable
{
  public:
    FunctionIF() {}

    FunctionIF(lua_IF* f, bool inside = true) { plus(f, inside); }

    FunctionIF(const FunctionIF& rhs) noexcept
    {
        m_f = rhs.m_f;
        action = rhs.action;
        invert = rhs.invert;
    }

    FunctionIF(FunctionIF&& rhs) noexcept = default;
    FunctionIF& operator=(const FunctionIF& rhs) = delete;
    FunctionIF& operator=(FunctionIF&& rhs) = delete;

    AMREX_GPU_HOST_DEVICE inline Real
      operator()(AMREX_D_DECL(Real x, Real y = 0.0, Real z = 0.0)) const noexcept
    {
        BL_PROFILE("FunctionIF::operator()");
        Real val = invert[0] * m_f[0]->operator()(AMREX_D_DECL(x, y, z));
        for (int i = 1; i < m_f.size(); ++i) {
            if (action[i] < 0) {
                val = max(val, invert[i] * m_f[i]->operator()(AMREX_D_DECL(x, y, z)));
            } else {
                val = min(val, invert[i] * m_f[i]->operator()(AMREX_D_DECL(x, y, z)));
            }
        }

        return val;
    }

    inline Real operator()(const RealArray& p) const noexcept
    {
        return this->operator()(AMREX_D_DECL(p[0], p[1], p[2]));
    }

    void plus(lua_IF* f, bool inside = true)
    {
        m_f.push_back(f);
        action.push_back(1);
        invert.push_back(inside ? 1.0 : -1.0);
    }

    void minus(lua_IF* f, bool inside = true)
    {
        m_f.push_back(f);
        action.push_back(-1);
        invert.push_back(inside ? 1.0 : -1.0);
    }

    bool empty() { return m_f.empty(); }

  protected:
    Vector<lua_IF*> m_f;
    Vector<int> action;
    Vector<int> invert;
};

// ============================================================================

class LuaSplineIF : public EB2::SplineIF
{
  public:
    LuaSplineIF() {};

    std::vector<amrex::RealVect> grab_points(const sol::table& points)
    {
        std::vector<amrex::RealVect> pts;

        int n_pts = points.size();
        pts.resize(n_pts);
        for (int pt_id = 0; pt_id < n_pts; ++pt_id) {
            sol::table pt = points[pt_id + 1];
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { pts[pt_id][d] = pt[d + 1].get<Real>(); }
        }

        return pts;
    }

    void addSpline(const sol::table& pts) { addSplineElement(grab_points(pts)); }

    void addLine(const sol::table& pts) { addLineElement(grab_points(pts)); }

    Real operator()(AMREX_D_DECL(Real x, Real y = 0.0, Real z = 0.0)) const noexcept
    {
        BL_PROFILE("LuaSplineIF::operator()");
        return EB2::SplineIF::operator()(RealArray({AMREX_D_DECL(x, y, z)}));
    }

    static void register_with_lua(sol::state& lua)
    {
        lua.new_usertype<LuaSplineIF>("SplineIF",
                                      sol::constructors<LuaSplineIF()>(),
                                      "addSpline",
                                      &LuaSplineIF::addSpline,
                                      "addLine",
                                      &LuaSplineIF::addLine,
                                      "query",
                                      &LuaSplineIF::operator());
    }
};

#endif
#endif  // MFP_EB_H
