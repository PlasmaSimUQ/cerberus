
#ifndef MFP_UTILITY_H
#define MFP_UTILITY_H

#include <algorithm>
#include <array>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <cassert>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>

// === PI ===
static const amrex::Real PI = 3.141592653589793;

//! Convert a number to a string
template <typename T>
std::string num2str(const T Number, const int precision=16) {
  std::ostringstream ss;
  ss.precision(precision);
  ss << Number;
  return ss.str();
}

// === MAX/MIN ===

template <typename T, size_t N>
T max_arr(const std::array<T, N> &in) {
    T v = in[0];
    for (size_t i = 0; i < in.size(); ++i) {
        v = std::max(in[i], v);
    }
    return v;
}

template <typename T>
T max_vec(const std::vector<T> &in) {
    T v = in[0];
    for (size_t i = 0; i < in.size(); ++i) {
        v = std::max(in[i], v);
    }
    return v;
}

template <typename T, size_t N>
T min_arr(const std::array<T, N> &in) {
    T v = in[0];
    for (size_t i = 0; i < in.size(); ++i) {
        v = std::min(in[i], v);
    }
    return v;
}

template <typename T>
T min_vec(const std::vector<T> &in) {
    T v = in[0];
    for (size_t i = 0; i < in.size(); ++i) {
        v = std::min(in[i], v);
    }
    return v;
}

// === VECTOR TO STRING ===

template <typename T>
std::string vec2str(const std::vector<T> &in, const std::string &left = "[",
                    const std::string &right = "]",
                    const std::string &delim = ", ") {
    std::string s = left;

    for (size_t i = 0; i < in.size(); ++i) {
        s += num2str(in[i]);
        if (i < in.size() - 1) {
            s += delim;
        }
    }
    s += right;

    return s;
}

// === ARRAY TO STRING ===

template <typename T, size_t N>
std::string arr2str(const std::array<T, N> &in, const std::string &left = "[",
                    const std::string &right = "]",
                    const std::string &delim = ", ") {
    std::string s = left;

    for (size_t i = 0; i < in.size(); ++i) {
        s += num2str(in[i]);
        if (i < in.size() - 1) {
            s += delim;
        }
    }
    s += right;

    return s;
}

// === ARRAY TO VECTOR ===
template <typename T, size_t N>
amrex::Vector<T> arr2vec(const amrex::Array<T, N> &in) {

    amrex::Vector<T> vec(N);
    for (size_t i = 0; i < N; ++i) {
        vec[i] = in[i];
    }

    return vec;
}

// === APPEND VECTOR ===

template <typename T>
void append(std::vector<T>& a, const std::vector<T>& b)
{
    a.reserve(a.size() + b.size());
    a.insert(a.end(), b.begin(), b.end());
}

// === CLAMP NUMBER WITHIN RANGE ===

template <typename T>
T clamp(const T& n, const T& lower, const T& upper) {
    return std::max(lower, std::min(n, upper));
}

template <typename T>
T circular_clamp(const T& n, const T& lower, const T& upper) {
    if (n < lower) {
        return upper;
    } else if (n > upper) {
       return lower;
    } else {
        return n;
    }
}

// === SIGN OF NUMBER ===

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T> T sign(T A, T B) {
    if (B >= T(0)) {
        return (T)std::abs(A);
    } else {
        return -(T)std::abs(A);
    }
}

// === FIND INDEX IN VECTOR ===
// https://thispointer.com/c-how-to-find-an-element-in-vector-and-get-its-index/

template < typename T>
std::pair<bool, int > findInVector(const std::vector<T>  & vecOfElements, const T  & element)
{
    std::pair<bool, int > result;

    // Find given element in vector
    auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);

    if (it != vecOfElements.end())
    {
        result.second = distance(vecOfElements.begin(), it);
        result.first = true;
    }
    else
    {
        result.first = false;
        result.second = -1;
    }

    return result;
}

template < typename T, class UnaryPredicate>
std::pair<bool, int > findInVector(const std::vector<T>  & vecOfElements, UnaryPredicate pred)
{
    std::pair<bool, int > result;

    // Find given element in vector
    auto it = std::find_if(vecOfElements.begin(), vecOfElements.end(), pred);

    if (it != vecOfElements.end())
    {
        result.second = distance(vecOfElements.begin(), it);
        result.first = true;
    }
    else
    {
        result.first = false;
        result.second = -1;
    }

    return result;
}

// === REMOVE FIRST ELEMENT OF VECTOR ===
template <typename V>
void pop_front(V & v)
{
    assert(!v.empty());
    v.erase(v.begin());
}

// === INDEXING ===

template <size_t N>
void multi_dim_index(const std::array<int,N> &lo,
                     const std::array<int,N> &hi,
                     std::vector<std::array<int,N>> &index, bool origin=true)
{
    // initialize the list
    if (index.size() == 0) {
        index.push_back(lo);
    }

    std::array<int,N> nxt = index.back();

    // exclude origin i.e. (0,0,0,...,0)
    if (!origin) {
        bool zero = true;
        for (int n=0; n<N; ++n) {
            if (nxt[n] != 0) {
                zero = false;
                break;
            }
        }
        if (zero) {
            index.pop_back();
        }
    }

    nxt[N-1] += 1;

    for (int n=N-1; n>0; --n) {

        if (nxt[n] > hi[n]) {
            nxt[n] = lo[n];
            nxt[n-1] += 1;
        } else {
            break;
        }
    }


    if (nxt[0] > hi[0]) {
        return;
    }

    index.push_back(nxt);

    multi_dim_index(lo, hi, index, origin);
}

// === 'SAFE' SQRT ===
template <typename V>
V safesqrt(V x)
{
    if (x > 0.0) {
        return std::sqrt(x);
    } else {
        return 0.0;
    }
}

// constexpr string
class str_const { // constexpr string
private:
    const char* const p_;
    const std::size_t sz_;
public:
    template<std::size_t N>
    constexpr str_const(const char(&a)[N]) : // ctor
        p_(a), sz_(N-1) {}
    constexpr char operator[](std::size_t n) { // []
        return n < sz_ ? p_[n] :
                         throw std::out_of_range("");
    }
    constexpr std::size_t size() { return sz_; } // size()
};

/*
 tanh(x)
Return the hyperbolic tangent of x.
Efficient tanh algorithm using an approximation of the Hyperbolic function:
Speed Improvement of the Back-Propagation on Current Generation
Workstations" D. Anguita, G. Parodi and R. Zunino. Proceedings of the World
Congress on Neural Networking, 1993.
*/
template <typename T>
T tanh_approx(T x)
{

    if (x > 1.92033) {
        return 0.96016;
    }

    if (x > 0) {
        return 0.96016 - 0.26037 * (x - 1.92033) * (x - 1.92033);
    }

    if (x <= -1.92033) {
        return -0.96016;
    }

    return 0.26037 * (x + 1.92033) * (x + 1.92033) - 0.96016;

}

// check if all elements in a vector are less than some value
template <typename T, typename N>
bool all_less(T first, T last, N test) {
    return std::all_of(first, last, [&test = test](N i) { return i < test;} );
}

// check if all elements in a vector are less than some value
template <typename T, typename N>
bool all_greater(T first, T last, N test) {
    return std::all_of(first, last, [&test = test](N i) { return i > test;} );
}

// check if all elements in a vector are less than some value
template <typename T, typename N>
bool all_equal(T first, T last, N test) {
    return std::all_of(first, last, [&test = test](N i) { return i == test;} );
}

// check if any element is equal to some value
template <typename T, typename N>
bool any_equal(T first, T last, N test) {
    return std::any_of(first, last, [&test = test](N i) { return i == test;} );
}

// conversion of enum class to underlying type
// usage = +Class::element
template <typename T>
constexpr auto operator+(T e) noexcept
-> std::enable_if_t<std::is_enum<T>::value, std::underlying_type_t<T>>
{
    return static_cast<std::underlying_type_t<T>>(e);
}



#endif
