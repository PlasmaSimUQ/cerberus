#ifndef COORDINATES_H
#define COORDINATES_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

using namespace amrex;

void transform_global2local(Vector<Real>& U, const int idir, const Vector<int>& vi);
void transform_local2global(Vector<Real>& U, const int idir, const Vector<int>& vi);

void transform_global2local(Vector<Real>& U,
            const Array<Array<Real,3>,3> C,
            const Vector<int>& vi);

void transform_local2global(Vector<Real>& U,
                const Array<Array<Real,3>,3> C,
                const Vector<int>& vi);

Array<Real,3> facefrac2normal(const Array<Array<Real,2>,AMREX_SPACEDIM> &alpha);
void expand_coord(Array<Array<Real,3>,3> &C);
Array<Array<Real,3>,3> normal2coord(const Array<Real,3> &normal);

#endif // COORDINATES_H
