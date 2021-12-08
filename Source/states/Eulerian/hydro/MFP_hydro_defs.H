#ifndef MFP_HYDRO_DEFS_H
#define MFP_HYDRO_DEFS_H

#include "MFP_utility.H"

using namespace amrex;

struct HydroDef
{
    enum class ConsIdx : int {
        Density=0,
        Xmom,
        Ymom,
        Zmom,
        Eden,
        NUM,
    };

    enum class PrimIdx : int {
        Density=0,
        Xvel,
        Yvel,
        Zvel,
        Prs,
        Temp,
        Gamma,
        SpHeat,
        NUM,
    };

    enum class VectorIdx : int {
        Cons=1,
        Prim=1
    };

    enum class WallIndex : int {
        HydroSlipWall = 0,
        HydroNoSlipWall,
        HydroDefined,
        NUM
    };
};


#endif // MFP_HYDRO_DEFS_H
