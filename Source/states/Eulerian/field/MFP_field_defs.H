#ifndef MFP_FIELD_DEFS_H
#define MFP_FIELD_DEFS_H

#include "MFP_utility.H"

struct FieldDef
{
    enum class ConsIdx : int {
        Dx=0,
        Dy,
        Dz,
        Bx,
        By,
        Bz,
        phi,
        psi,
        mu,
        ep,
        NUM
    };

    enum class VectorIdx : int {
        Cons=2,
    };
};

#endif // MFP_FIELD_DEFS_H
