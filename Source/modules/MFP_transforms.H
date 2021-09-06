#ifndef COORDINATES_H
#define COORDINATES_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include "AMReX_BLProfiler.H"

using namespace amrex;

template<typename T1, typename T2>
void transform_global2local(T1& U, const int idir, const T2& vi)
{
    BL_PROFILE("transform_global2local");
    // iterate over the passed in index and rotate the 3-vectors in-place

    if (idir == 0) {
        return;
    }

    Real n, t1, t2;

    int in  = (idir+0)%3;
    int it1 = (idir+1)%3;
    int it2 = (idir+2)%3;

    for (const int& v : vi) {
        n  = U[v+in];
        t1 = U[v+it1];
        t2 = U[v+it2];

        U[v+0] = n;
        U[v+1] = t1;
        U[v+2] = t2;
    }
}

template<typename T1, typename T2>
void transform_local2global(T1& U, const int idir, const T2& vi)
{
    BL_PROFILE("transform_local2global");
    if (idir == 0) {
        return;
    }

    Real n, t1, t2;

    int in  = (idir+0)%3;
    int it1 = (idir+1)%3;
    int it2 = (idir+2)%3;

    for (const int& v : vi) {

        n  = U[v+0];
        t1 = U[v+1];
        t2 = U[v+2];

        U[v+in] = n;
        U[v+it1] = t1;
        U[v+it2] = t2;
    }
}

template<typename T1, typename T2>
void transform_global2local(T1 &U,
            const Array<Array<Real,3>,3> C,
            T2 &vi)
{
    BL_PROFILE("transform_global2local");
    // iterate over the passed in index and rotate the 3-vectors in-place

    Real x, y, z;

    for (const int& v : vi) {
        x = U[v+0];
        y = U[v+1];
        z = U[v+2];

        U[v+0] = x*C[0][0] + y*C[0][1] + z*C[0][2];
        U[v+1] = x*C[1][0] + y*C[1][1] + z*C[1][2];
        U[v+2] = x*C[2][0] + y*C[2][1] + z*C[2][2];
    }
}

template<typename T1, typename T2>
void transform_local2global(T1 &U,
                const Array<Array<Real,3>,3> C,
                T2 &vi)
{
    BL_PROFILE("transform_local2global");
    Real x, y, z;

    for (const int& v : vi) {
        x = U[v+0];
        y = U[v+1];
        z = U[v+2];

        U[v+0] = x*C[0][0] + y*C[1][0] + z*C[2][0];
        U[v+1] = x*C[0][1] + y*C[1][1] + z*C[2][1];
        U[v+2] = x*C[0][2] + y*C[1][2] + z*C[2][2];
    }

}


Array<Real,3> facefrac2normal(const Array<Array<Real,2>,AMREX_SPACEDIM> &alpha);
void expand_coord(Array<Array<Real,3>,3> &C);
Array<Array<Real,3>,3> normal2coord(const Array<Real,3> &normal);

#endif // COORDINATES_H
