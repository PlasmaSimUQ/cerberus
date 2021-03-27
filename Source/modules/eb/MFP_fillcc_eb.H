#ifndef FAB_FILLCC_EB_H
#define FAB_FILLCC_EB_H

#include <AMReX_Box.H>
#include <AMReX_BCRec.H>

namespace amrex {

template <typename T>
void fab_filcc_eb (Box const& bx, Array4<T> const& qn, int ncomp,
                Box const& domain, Real const* /*dx*/, Real const* /*xlo*/,
                BCRec const* bcn, const T& check)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const auto domlo = amrex::lbound(domain);
    const auto domhi = amrex::ubound(domain);

    const int ilo = domlo.x;
    const int ihi = domhi.x;

#if AMREX_SPACEDIM >= 2
    const int jlo = domlo.y;
    const int jhi = domhi.y;
#endif

#if AMREX_SPACEDIM == 3
    const int klo = domlo.z;
    const int khi = domhi.z;
#endif

    for (int n = 0; n < ncomp; ++n)
    {
        Array4<T> q(qn,n);
        BCRec const& bc = bcn[n];

        if (lo.x < ilo) {
           const int imin = lo.x;
           const int imax = ilo-1;
           if (bc.lo(0) == BCType::reflect_even) {
               for (int k = lo.z; k <= hi.z; ++k) {
               for (int j = lo.y; j <= hi.y; ++j) {
               for (int i = imin; i <= imax; ++i) {
                   if (q(ilo+(ilo-i)-1,j,k) != check)
                   q(i,j,k) = q(ilo+(ilo-i)-1,j,k);
               }}}
           }
        }

        if (hi.x > ihi) {
            const int imin = ihi+1;
            const int imax = hi.x;

            if (bc.hi(0) == BCType::reflect_even) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = imin; i <= imax; ++i) {
                    if (q(ihi-(i-ihi)+1,j,k) != check)
                    q(i,j,k) = q(ihi-(i-ihi)+1,j,k);
                }}}
            }
        }

#if AMREX_SPACEDIM >= 2

        if (lo.y < jlo) {
            const int jmin = lo.y;
            const int jmax = jlo-1;
            if (bc.lo(1) == BCType::reflect_even) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (q(i,jlo+(jlo-j)-1,k) != check)
                    q(i,j,k) = q(i,jlo+(jlo-j)-1,k);
                }}}
            }
        }

        if (hi.y > jhi) {
            const int jmin = jhi+1;
            const int jmax = hi.y;
            if (bc.hi(1) == BCType::reflect_even) {
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = jmin; j <= jmax; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (q(i,jhi-(j-jhi)+1,k) != check)
                    q(i,j,k) = q(i,jhi-(j-jhi)+1,k);
                }}}
            }
        }
#endif

#if AMREX_SPACEDIM == 3

        if (lo.z < klo) {
            const int kmin = lo.z;
            const int kmax = klo-1;
            if (bc.lo(2) == BCType::reflect_even) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (q(i,j,klo+(klo-k)-1) != check)
                    q(i,j,k) = q(i,j,klo+(klo-k)-1);
                }}}
            }
        }

        if (hi.z > khi) {
            const int kmin = khi+1;
            const int kmax = hi.z;
            if (bc.hi(2) == BCType::reflect_even) {
                for (int k = kmin; k <= kmax; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (q(i,j,khi-(k-khi)+1) != check)
                    q(i,j,k) = q(i,j,khi-(k-khi)+1);
                }}}
            }
        }
#endif

    }
}

}

#endif // FAB_FILLCC_EB_H
