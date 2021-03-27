#include "MFP_transforms.H"

#include <math.h>

void transform_global2local(Vector<Real> &U, const int idir, const Vector<int> &vi)
{

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

void transform_local2global(Vector<Real> &U, const int idir, const Vector<int> &vi)
{

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

void transform_global2local(Vector<Real> &U,
            const Array<Array<Real,3>,3> C,
            const Vector<int> &vi)
{

    // iterate over the passed in index and rotate the 3-vectors in-place

    Real x, y, z;

    for (const int& v : vi) {
        x = U[v+0];
        y = U[v+1];
        z = U[v+2];

        U[v+0] = x*C[0][0]  + y*C[0][1]  + z*C[0][2];
        U[v+1] = x*C[1][0] + y*C[1][1] + z*C[1][2];
        U[v+2] = x*C[2][0] + y*C[2][1] + z*C[2][2];
    }
}

void transform_local2global(Vector<Real> &U,
                const Array<Array<Real,3>,3> C,
                const Vector<int> &vi)
{

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

Array<Real,3> facefrac2normal(const Array<Array<Real,2>,AMREX_SPACEDIM> &alpha)
{
    Array<Real,3> normal = {0,0,0};

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        normal[d] = alpha[d][0] - alpha[d][1]; // pointing to the wall
    }

    Real apnorm = sqrt(AMREX_D_TERM(normal[0]*normal[0], + normal[1]*normal[1], + normal[2]*normal[2]));

    if (apnorm == 0.0) {
        Abort("Error: wall normal length = 0");
    }

    Real apnorminv = 1.0 / apnorm;
    normal[0] *= apnorminv;
    normal[1] *= apnorminv;
    normal[2] *= apnorminv;

    return normal;
}

void expand_coord(Array<Array<Real,3>,3> &coord)
{
    // assume coord[0] is the normal
    // now get an orthogonal vector
#if AMREX_SPACEDIM == 2
    coord[1] = {-coord[0][1], coord[0][0], 0.0};
#else
    const Real g = std::copysign(1.0, coord[0][2]);
    const Real h = coord[0][2] + g;
    coord[1] = {g*h - coord[0][0]*coord[0][0], -coord[0][0]*coord[0][1], -coord[0][0]*h};
#endif

    // and complete the orthogonal basis by taking the cross product
    coord[2] = {
        coord[0][1] * coord[1][2] - coord[0][2] * coord[1][1],
        -(coord[0][0] * coord[1][2] - coord[0][2] * coord[1][0]),
        coord[0][0] * coord[1][1] - coord[0][1] * coord[1][0]
    };

}

Array<Array<Real,3>,3> normal2coord(const Array<Real,3> &normal)
{
    Array<Array<Real,3>,3> coord = {{normal,{{0.0,0.0,0.0}},{{0.0,0.0,0.0}}}};

    expand_coord(coord);

    return coord;
}

