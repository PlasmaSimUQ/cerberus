#ifndef MFP_LUA_H
#define MFP_LUA_H

#include "sol.hpp"

template <typename T, size_t N>
void expand(const sol::object &in, std::array<T,N> &out) {
    if (in.get_type() == sol::type::number) {
        out.fill(in.as<T>());
    } else if (in.get_type() == sol::type::table) {
        const sol::table m = in.as<sol::table>();
        for (int i=0; i<N; ++i) {
            out[i] = m[i+1];
        }
    }
}

template <typename T>
void get_numbered(const sol::table &in, std::vector<T> &out)
{
    if (!in.valid())
        return;

    for (const auto& v : in) {
        sol::object key = v.first;
        sol::object value = v.second;

        if (key.get_type() == sol::type::number) {
            const int i = key.as<int>() - 1; // account for 1 based lua indexing
            if (i >= out.size()) {
                out.resize(i+1);
            }
            out[i] = value.as<T>();
        }
    }
}

#endif // MFP_LUA_H
