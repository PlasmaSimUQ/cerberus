#include "MFP_optional_func.H"

#include <AMReX_BLProfiler.H>

void Optional3D1VFunction::set_value(const Real& val)
{
    BL_PROFILE("Optional3D1VFunction::set_value");

    value = val;
    valid = true;
}
void Optional3D1VFunction::set_func(const opt_func& fun)
{
    BL_PROFILE("Optional3D1VFunction::set_func");

    f = fun;
    valid = true;
}

bool Optional3D1VFunction::is_valid() const
{
    BL_PROFILE("Optional3D1VFunction::is_valid");

    return valid;
}

bool Optional3D1VFunction::has_func() const
{
    BL_PROFILE("Optional3D1VFunction::has_func");

    return f.operator bool();
}

Real Optional3D1VFunction::operator()(
  Real x, Real y, Real z, Real t, const Vector<std::string>& names, const Vector<Real>& data) const
{
    BL_PROFILE("Optional3D1VFunction::operator()");

    if (f) {
        std::map<std::string, Real> in = {
          {"x", x},
          {"y", y},
          {"z", z},
          {"t", t},
        };

        for (int i = 0; i < names.size(); ++i) { in[names[i]] = data[i]; }

        return f(in);
    } else {
        return value;
    }
}
Real Optional3D1VFunction::operator()(const std::map<std::string, Real>& data) const
{
    BL_PROFILE("Optional3D1VFunction::operator()");

    if (f) {
        return f(data);
    } else {
        return value;
    }
}

std::ostream& operator<<(std::ostream& os, const Optional3D1VFunction& f)
{
    if (f.f) {
        os << "f(" << &f.f << ")";
    } else if (f.valid) {
        os << f.value;
    } else {
        os << "*";
    }
    return os;
}

bool get_udf(const sol::object& obj, Optional3D1VFunction& v, const Real fallback)
{
    if (!obj.valid()) {
        v.set_value(fallback);
        return false;
    }

    sol::type tp = obj.get_type();

    switch (tp) {
    case sol::type::lua_nil: v.set_value(fallback); return false;
    case sol::type::number: v.set_value(obj.as<Real>()); return true;
    case sol::type::function: v.set_func(obj.as<opt_func>()); return true;
    default: return false;
    }
}

Optional3D1VFunction get_udf(const sol::object& obj)
{
    Optional3D1VFunction v;

    if (!obj.valid()) { return v; }

    sol::type tp = obj.get_type();

    switch (tp) {
    case sol::type::number: v.set_value(obj.as<Real>()); return v;
    case sol::type::function: v.set_func(obj.as<opt_func>()); return v;
    default: return v;
    }
}
