#ifndef MFP_MHD_DEFS_H
#define MFP_MHD_DEFS_H

#include "MFP_utility.H"

struct MHDDef
{
    enum class ConsIdx : int {
        Density=0,
        Xmom,
        Ymom,
        Zmom,
        Eden,
        Bx,
        By,
        Bz,
        psi,
        NUM,
    };

    enum class PrimIdx : int {
        Density=0,
        Xvel,
        Yvel,
        Zvel,
        Prs,
        Bx,
        By,
        Bz,
        psi,
        NUM,
    };

    enum class VectorIdx : int {
        Cons=1,
        Prim=1
    };

    enum class WallIndex : int {
        MHDSlipWall = 0,
        MHDDefined,
        NUM
    };
};


#endif // MFP_MHD_DEFS_H
