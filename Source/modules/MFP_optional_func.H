#ifndef MFP_TYPES_H
#define MFP_TYPES_H

#include "sol.hpp"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <functional>
#include <iostream>
#include <map>

using namespace amrex;

// Optional3D1VFunction provides a container for holding either a constant
// value or a function.

typedef std::function<Real(std::map<std::string, Real>)> opt_func;

class Optional3D1VFunction
{
  public:
    bool is_valid() const;
    void set_value(const Real& val);
    void set_func(const opt_func& fun);
    bool has_func() const;
    Real operator()(Real x = 0,
                    Real y = 0,
                    Real z = 0,
                    Real t = 0,
                    const Vector<std::string>& names = {},
                    const Vector<Real>& data = {}) const;
    Real operator()(const std::map<std::string, Real>& data = {}) const;
    friend std::ostream& operator<<(std::ostream& os, const Optional3D1VFunction& f);

  private:
    bool valid = false;
    Real value = 0.0;
    opt_func f;
};

std::ostream& operator<<(std::ostream& os, const Optional3D1VFunction& f);

bool get_udf(const sol::object& obj, Optional3D1VFunction& v, const Real fallback);
Optional3D1VFunction get_udf(const sol::object& obj);

#endif  // MFP_TYPES_H
