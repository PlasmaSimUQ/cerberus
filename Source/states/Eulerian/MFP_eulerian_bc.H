#ifndef EULERIANEB_H
#define EULERIANEB_H

#include "AMReX_Array.H"
#include "AMReX_Array4.H"
#include "AMReX_Vector.H"

using namespace amrex;

class EulerianBoundaryEB
{
  public:
    EulerianBoundaryEB();
    ~EulerianBoundaryEB();

    enum EBType { Hydro = 0, Field, MHD };

    EBType type;

    int state_idx;

    virtual const int get_type() const = 0;

    virtual bool get_slope() const { return false; }

    virtual void update_data_idx(Vector<size_t>& glob2dat) { state_idx = glob2dat[state_idx]; }

    virtual void solve(Array<Array<Real, 3>, 3>& wall_coord,
                       Array<Real, AMREX_SPACEDIM> wall_centre,
                       const Vector<Array4<const Real>>& all_prim,
                       const int i,
                       const int j,
                       const int k,
                       const Real* dx,
                       Array<Vector<Real>, AMREX_SPACEDIM>& F)
    {
    }

    virtual Vector<Vector<Real>> get_wall_state(const Array<Real, AMREX_SPACEDIM> wall_centre,
                                                const Array<Array<Real, 3>, 3>& wall_coord,
                                                const Real t) const
    {
        return {{}};
    }
};

#endif  // EULERIANEB_H
