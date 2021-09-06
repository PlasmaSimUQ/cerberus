#include "MFP_transforms.H"

#include <math.h>


Array<Real,3> facefrac2normal(const Array<Array<Real,2>,AMREX_SPACEDIM> &alpha)
{
    BL_PROFILE("facefrac2normal");
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
    BL_PROFILE("expand_coord");
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
    BL_PROFILE("normal2coord");
    Array<Array<Real,3>,3> coord = {{normal,{{0.0,0.0,0.0}},{{0.0,0.0,0.0}}}};

    expand_coord(coord);

    return coord;
}

