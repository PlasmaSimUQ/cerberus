#ifndef MFP_LUA_H
#define MFP_LUA_H

#include "sol.hpp"

template <typename T> void set_values(const sol::object& in, std::vector<T>& out)
{
    if (!in.valid()) return;
    if (in.get_type() == sol::type::table) {
        const sol::table m = in.as<sol::table>();
        out.resize(m.size());
        for (int i = 0; i < m.size(); ++i) { out[i] = m[i + 1]; }
    } else {
        out.push_back(in.as<T>());
    }
}

template <typename T> void get_numbered(const sol::table& in, std::vector<T>& out)
{
    if (!in.valid()) return;

    for (const auto& v : in) {
        sol::object key = v.first;
        sol::object value = v.second;

        if (key.get_type() == sol::type::number) {
            const int i = key.as<int>() - 1;  // account for 1 based lua indexing
            if (i >= out.size()) { out.resize(i + 1); }
            out[i] = value.as<T>();
        }
    }
}

#endif  // MFP_LUA_H
