#ifndef MFP_RK4_H
#define MFP_RK4_H

#include "AMReX_Array.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "MFP_utility.H"

using amrex::Array;
using amrex::Real;
using amrex::Vector;

template <size_t N, size_t M>
bool rk4_fixed(Real& t,
               Array<Real, N>& y,
               Array<Real, M>& data,
               int (*rhs_eval)(Real, Array<Real, N>, Array<Real, N>&, Array<Real, M>&),
               bool (*check)(Array<Real, N>&, Array<Real, M>&),
               Real h)
{
    // Declare local arrays
    Array<Real, N> k1, k2, k3, k4, f, dydt;

    // Zeroth intermediate step
    if ((*rhs_eval)(t, y, dydt, data)) return true;
    for (int j = 0; j < N; j++) {
        k1[j] = h * dydt[j];
        f[j] = y[j] + k1[j] / 2.;
    }

    // First intermediate step
    if ((*rhs_eval)(t + h / 2., f, dydt, data)) return true;
    for (int j = 0; j < N; j++) {
        k2[j] = h * dydt[j];
        f[j] = y[j] + k2[j] / 2.;
    }

    // Second intermediate step
    if ((*rhs_eval)(t + h / 2., f, dydt, data)) return true;
    for (int j = 0; j < N; j++) {
        k3[j] = h * dydt[j];
        f[j] = y[j] + k3[j];
    }

    // Third intermediate step
    if ((*rhs_eval)(t + h, f, dydt, data)) return true;
    for (int j = 0; j < N; j++) { k4[j] = h * dydt[j]; }

    // Actual step
    for (int j = 0; j < N; j++) { y[j] += k1[j] / 6. + k2[j] / 3. + k3[j] / 3. + k4[j] / 6.; }

    if ((*check)(y, data)) { return true; }

    t += h;

    return false;
}

template <size_t N, size_t M>
void rk4_adaptive(Real& t,
                  Array<Real, N>& y,
                  Array<Real, M>& data,
                  int (*rhs_eval)(Real, Array<Real, N>, Array<Real, N>&, Array<Real, M>&),
                  bool (*check)(Array<Real, N>&, Array<Real, M>&),
                  Real& h,
                  Real S,
                  int& depth,
                  int max_depth)
{
    // Declare local arrays
    Array<Real, N> y0;

    // Declare repetition counter
    static int count = 0;

    // Save initial data
    Real x0 = t;
    y0 = y;

    // Take full step
    if (rk4_fixed(t, y, data, rhs_eval, check, h)) {
        // step failed, reduce time step and try again

        // Print() << "rk4_fixed failed \n";
        h /= S;

        // Abort if step-length falls below h_min
        if (depth > max_depth) {
            amrex::Abort("exceeded time step refinement limitiation in rk4_adaptive()");
        }

        depth += 1;

        count++;
        t = x0;
        y = y0;

        rk4_adaptive(t, y, data, rhs_eval, check, h, S, depth, max_depth);

        depth -= 1;
    }

    return;
}

#endif  // MFP_RK4_H
