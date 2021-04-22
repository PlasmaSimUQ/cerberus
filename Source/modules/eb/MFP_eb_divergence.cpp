#ifdef AMREX_USE_EB
#include "MFP_eb_divergence.H"

#include <AMReX_YAFluxRegister.H>
using CellType = YAFluxRegister::CellType;

#include "MFP_state.H"
#include "MFP_global.H"
#include "sol.hpp"

#ifdef PYTHON
#include "matplotlibcpp.h"
#include "MFP_diagnostics.H"
namespace plt = matplotlibcpp;
#endif

using GD = GlobalData;

DivergenceEB::DivergenceEB()
{
}

DivergenceEB::~DivergenceEB()
{
}

void DivergenceEB::calc_eb_divergence(const Box& box,
                                      const FArrayBox &cons,
                                      Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                      Array<FArrayBox, AMREX_SPACEDIM> &wall_fluxes,
                                      FArrayBox& du,
                                      const EBCellFlagFab& flag,
                                      const FArrayBox& vfrac,
                                      const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                                      const Array<const FArrayBox *, AMREX_SPACEDIM> &fcent,
                                      int as_crse,
                                      int as_fine,
                                      const IArrayBox *rrflag_as_crse,
                                      const IArrayBox &levmsk,
                                      FArrayBox *rr_drho_crse,
                                      FArrayBox &dm_as_fine,
                                      const Real *dx,
                                      const Real dt) const
{
    // do nothing
}

void DivergenceEB::merge_cells(const Box& box,
                               FArrayBox &cons,
                               FArrayBox& du,
                               const EBCellFlagFab& flag,
                               const FArrayBox& vfrac,
                               const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                               int as_fine,
                               FArrayBox &dm_as_fine,
                               const IArrayBox& levmsk) const
{
    // do nothing
}

bool DivergenceEB::is_inside(const int i,const int j, const int k, const Dim3 &lo, const Dim3 &hi)
{
    return  i >= lo.x && i <= hi.x &&
            j >= lo.y && j <= hi.y &&
            k >= lo.z && k <= hi.z;
}

PhysicsFactory<DivergenceEB>& GetDivergenceEBBuilder()
{
    static PhysicsFactory<DivergenceEB> F;
    return F;
}

//=============================================================================

std::string RedistributeEB::tag = "redistribute";
bool RedistributeEB::registered = GetDivergenceEBBuilder().Register(RedistributeEB::tag, DivergenceEBBuilder<RedistributeEB>);

Vector<std::string> RedistributeEB::options = {"uniform", "volume", "density", "energy"};

RedistributeEB::RedistributeEB()
{
}

RedistributeEB::RedistributeEB(const sol::table &def)
{
    global_idx = def["global_idx"];

    // get the redistribution strategy

    const sol::table &div_def = def["eb_divergence"];

    std::string redist = div_def.get_or<std::string>("strategy", "volume");

    const auto found = findInVector(options,redist);
    if (found.first) {
        redistribution_strategy = (RedistributionEB)found.second;
    } else {
        Abort("Invalid redistribution option '"+redist+"', valid options are "+vec2str(options));

    }

    if (redist == "density") {
        State &istate = GD::get_state(global_idx);
        if (istate.get_cons_density_idx() < 0)
            Abort("Invalid redistribution option '"+redist+"' for state "+istate.name+" as it does not have a density");
    }

    // get the reredistribution threshold
    reredistribution_threshold = div_def.get_or("reredistribution_threshold", 1.0e-14);

}

RedistributeEB::~RedistributeEB()
{
}

void RedistributeEB::calc_eb_divergence(const Box& box,
                                        const FArrayBox &cons,
                                        Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                        Array<FArrayBox, AMREX_SPACEDIM> &wall_fluxes,
                                        FArrayBox& du,
                                        const EBCellFlagFab& flag,
                                        const FArrayBox& vfrac,
                                        const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                                        const Array<const FArrayBox *, AMREX_SPACEDIM> &fcent,
                                        int as_crse, int as_fine,
                                        const IArrayBox *rrflag_as_crse,
                                        const IArrayBox &levmsk,
                                        FArrayBox *rr_drho_crse,
                                        FArrayBox &dm_as_fine,
                                        const Real *dx,
                                        const Real dt) const
{
    BL_PROFILE("RedistributeEB::calc_eb_divergence");
    State &istate = GD::get_state(global_idx);

    // make sure arrays are empty
    du.setVal(0.0);
    dm_as_fine.setVal(0.0);

    int nc = istate.n_cons();

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<int,3> index = {0,0,0};

    Array4<const Real> const& cons4 = cons.array();
    Array4<Real> const& du4 = du.array();

    Array4<const int> const& rr_flag_crse4 = rrflag_as_crse->array();
    Array4<const int> const& levmsk4 = levmsk.array();
    Array4<Real> const& rr_drho_crse4 = rr_drho_crse->array();
    Array4<Real> const& dm_as_fine4 = dm_as_fine.array();

    Array<Array4<Real>,AMREX_SPACEDIM> flux4, wflux4;
    Array4<const EBCellFlag> const& flag4 = flag.array();
    Array4<const Real> const& vfrac4 = vfrac.array();
    Array<Array4<const Real>,AMREX_SPACEDIM> afrac4;
    Array<Array4<const Real>,AMREX_SPACEDIM> fcent4;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        flux4[d] = fluxes[d].array();
        wflux4[d] = wall_fluxes[d].array();
        afrac4[d] = afrac[d]->array();
        fcent4[d] = fcent[d]->array();
    }

    FArrayBox divc(grow(box,2),nc); // check size
    Array4<Real> const& divc4 = divc.array();

    FArrayBox rediswgt(grow(box,2)); // check size
    Array4<Real> const& rediswgt4 = rediswgt.array();

    Real dxinv[AMREX_SPACEDIM] = {AMREX_D_DECL(dt/dx[0],dt/dx[1],dt/dx[2])};

    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);


    for     (int k = lo.z-AMREX_D_PICK(0,0,2); k <= hi.z+AMREX_D_PICK(0,0,2); ++k) {
        for   (int j = lo.y-AMREX_D_PICK(0,2,2); j <= hi.y+AMREX_D_PICK(0,2,2); ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x-2; i <= hi.x+2; ++i) {

                const EBCellFlag &cflag = flag4(i,j,k);
                Real vf = vfrac4(i,j,k);

                if (vf <= 0.0) {
                    for (int n=0; n<nc; ++n) {
                        divc4(i,j,k,n) = 0.0;
                    }

                } else if (vf >= 1.0) {

                    for (int n=0; n<nc; ++n) {
                        Real div = 0.0;
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            index[d] = 1;
                            div += dxinv[d]*(flux4[d](i,j,k,n) - flux4[d](i+index[0], j+index[1], k+index[2], n));
                            index[d] = 0;
                        }
                        divc4(i,j,k,n) = div;
                    }

                } else {

                    // contribution from the standard fluxes

                    for (int n=0; n<nc; ++n) {


                        /* Use a weighted combination of fluxes from the surrounding faces
                             * to calculate the flux for this face, rather than just use the
                             * area-fraction of the flux.
                             * e.g. in 3D, for face (i,j) that has a cut in the top right we use a combination
                             * of the fluxes on faces to the left (i-1,j), below-left (i-1,j-1), and
                             * below (i,j-1)
                             *
                             *  -------- -------
                             * |        |    \**|
                             * | i-1,j  | i,j \*|
                             * |        |       |
                             *  -------- -------
                             * |        |       |
                             * |i-1,j-1 | i,j-1 |
                             * |        |       |
                             *  -------- -------
                             */

                        Array<Array<Real,2>,AMREX_SPACEDIM> flux;

                        for (int lh=0; lh<2; ++lh) {
#if AMREX_SPACEDIM == 2
                            // x
                            if (fcent4[0](i+lh,j,k,0) <= 0.0) {
                                Real fracy = -fcent4[0](i+lh,j,k,0)*cflag.isConnected(0, -1, 0);
                                flux[0][lh] = fracy*flux4[0](i+lh,j-1,k,n) + (1.0-fracy)*flux4[0](i+lh,j,k,n);
                            } else {
                                Real fracy = fcent4[0](i+lh,j,k,0)*cflag.isConnected(0, 1, 0);
                                flux[0][lh] = fracy*flux4[0](i+lh,j+1,k,n) + (1.0-fracy)*flux4[0](i+lh,j,k,n);
                            }

                            // y
                            if (fcent4[1](i,j+lh,k,0) <= 0.0) {
                                Real fracx = -fcent4[1](i,j+lh,k,0)*cflag.isConnected(-1, 0, 0);
                                flux[1][lh] = fracx*flux4[1](i-1,j+lh,k,n) + (1.0-fracx)*flux4[1](i,j+lh,k,n);

                            } else {
                                Real fracx = fcent4[1](i,j+lh,k,0)*cflag.isConnected(1, 0, 0);
                                flux[1][lh] = fracx*flux4[1](i+1,j+lh,k,n) + (1.0-fracx)*flux4[1](i,j+lh,k,n);

                            }

#elif AMREX_SPACEDIM == 3
                            // x
                            if (fcent4[0](i+lh,j,k,0) <= 0.0) {
                                Real fracy = -fcent4[0](i+lh,j,k,0)*cflag.isConnected(0, -1, 0);
                                if (fcent4[0](i+lh,j,k,1) <= 0.0) {
                                    Real fracz = -fcent4[0](i+lh,j,k,1)*cflag.isConnected(0, 0, -1);
                                    flux[0][lh] = (1.0-fracz)*(fracy*flux4[0](i+lh,j-1,k  ,n) + (1.0-fracy)*flux4[0](i+lh,j  ,k  ,n)) +
                                            fracz *(fracy*flux4[0](i+lh,j-1,k-1,n) + (1.0-fracy)*flux4[0](i+lh,j  ,k-1,n));
                                } else {
                                    Real fracz = fcent4[0](i+lh,j,k,1)*cflag.isConnected(0, 0, 1);
                                    flux[0][lh] = (1.0-fracz)*(fracy*flux4[0](i+lh,j-1,k  ,n) + (1.0-fracy)*flux4[0](i+lh,j  ,k  ,n)) +
                                            fracz *(fracy*flux4[0](i+lh,j-1,k+1,n) + (1.0-fracy)*flux4[0](i+lh,j  ,k+1,n));
                                }
                            } else {
                                Real fracy = fcent4[0](i+lh,j,k,0)*cflag.isConnected(0, 1, 0);
                                if (fcent4[0](i+lh,j,k,1) <= 0.0) {
                                    Real fracz = -fcent4[0](i+lh,j,k,1)*cflag.isConnected(0, 0, -1);
                                    flux[0][lh] = (1.0-fracz)*(fracy*flux4[0](i+lh,j+1,k  ,n) + (1.0-fracy)*flux4[0](i+lh,j  ,k  ,n)) +
                                            fracz *(fracy*flux4[0](i+lh,j+1,k-1,n) + (1.0-fracy)*flux4[0](i+lh,j  ,k-1,n));
                                } else {
                                    Real fracz = fcent4[0](i+lh,j,k,1)*cflag.isConnected(0, 0, 1);
                                    flux[0][lh] = (1.0-fracz)*(fracy*flux4[0](i+lh,j+1,k  ,n) + (1.0-fracy)*flux4[0](i+lh,j  ,k  ,n)) +
                                            fracz *(fracy*flux4[0](i+lh,j+1,k+1,n) + (1.0-fracy)*flux4[0](i+lh,j  ,k+1,n));
                                }
                            }

                            // y
                            if (fcent4[1](i,j+lh,k,0) <= 0.0) {
                                Real fracx = -fcent4[1](i,j+lh,k,0)*cflag.isConnected(-1, 0, 0);
                                if (fcent4[1](i,j+lh,k,1) <= 0.0) {
                                    Real fracz = -fcent4[1](i,j+lh,k,1)*cflag.isConnected(0, 0, -1);
                                    flux[1][lh] = (1.0-fracz)*(fracx*flux4[1](i-1,j+lh,k  ,n) + (1.0-fracx)*flux4[1](i,j+lh,k  ,n)) +
                                            fracz *(fracx*flux4[1](i-1,j+lh,k-1,n) + (1.0-fracx)*flux4[1](i,j+lh,k-1,n));
                                } else {
                                    Real fracz = fcent4[1](i,j+lh,k,1)*cflag.isConnected(0, 0, 1);
                                    flux[1][lh] = (1.0-fracz)*(fracx*flux4[1](i-1,j+lh,k  ,n) + (1.0-fracx)*flux4[1](i,j+lh,k  ,n)) +
                                            fracz *(fracx*flux4[1](i-1,j+lh,k+1,n) + (1.0-fracx)*flux4[1](i,j+lh,k+1,n));
                                }
                            } else {
                                Real fracx = fcent4[1](i,j+lh,k,0)*cflag.isConnected(1, 0, 0);
                                if (fcent4[1](i,j+lh,k,1) <= 0.0) {
                                    Real fracz = -fcent4[1](i,j+lh,k,1)*cflag.isConnected(0, 0, -1);
                                    flux[1][lh] = (1.0-fracz)*(fracx*flux4[1](i+1,j+lh,k  ,n) + (1.0-fracx)*flux4[1](i,j+lh,k  ,n)) +
                                            fracz *(fracx*flux4[1](i+1,j+lh,k-1,n) + (1.0-fracx)*flux4[1](i,j+lh,k-1,n));
                                } else {
                                    Real fracz = fcent4[1](i,j+lh,k,1)*cflag.isConnected(0, 0, 1);
                                    flux[1][lh] = (1.0-fracz)*(fracx*flux4[1](i+1,j+lh,k  ,n) + (1.0-fracx)*flux4[1](i,j+lh,k  ,n)) +
                                            fracz *(fracx*flux4[1](i+1,j+lh,k+1,n) + (1.0-fracx)*flux4[1](i,j+lh,k+1,n));
                                }
                            }

                            // z
                            if (fcent4[2](i,j,k+lh,0) <= 0.0) {
                                Real fracx = -fcent4[2](i,j,k+lh,0)*cflag.isConnected(-1, 0, 0);
                                if (fcent4[2](i,j,k+lh,1) <= 0.0) {
                                    Real fracy = -fcent4[2](i,j,k+lh,1)*cflag.isConnected(0, -1, 0);
                                    flux[2][lh] = (1.0-fracy)*(fracx*flux4[2](i-1,j  ,k+lh,n) + (1.0-fracx)*flux4[2](i,j,  k+lh,n)) +
                                            fracy *(fracx*flux4[2](i-1,j-1,k+lh,n) + (1.0-fracx)*flux4[2](i,j-1,k+lh,n));
                                } else {
                                    Real fracy = fcent4[2](i,j,k+lh,1)*cflag.isConnected(0, 1, 0);
                                    flux[2][lh] = (1.0-fracy)*(fracx*flux4[2](i-1,j  ,k+lh,n) + (1.0-fracx)*flux4[2](i,j  ,k   ,n)) +
                                            fracy *(fracx*flux4[2](i-1,j+1,k+lh,n) + (1.0-fracx)*flux4[2](i,j+1,k+lh,n));
                                }
                            } else {
                                Real fracx = fcent4[2](i,j,k+lh,0)*cflag.isConnected(1, 0, 0);
                                if (fcent4[2](i,j,k+lh,1) <= 0.0) {
                                    Real fracy = -fcent4[2](i,j,k+lh,1)*cflag.isConnected(0, -1, 0);
                                    flux[2][lh] = (1.0-fracy)*(fracx*flux4[2](i+1,j  ,k+lh,n) + (1.0-fracx)*flux4[2](i,j,  k+lh,n)) +
                                            fracy *(fracx*flux4[2](i+1,j-1,k+lh,n) + (1.0-fracx)*flux4[2](i,j-1,k+lh,n));
                                } else {
                                    Real fracy = fcent4[2](i,j,k+lh,1)*cflag.isConnected(0, 1, 0);
                                    flux[2][lh] = (1.0-fracy)*(fracx*flux4[2](i+1,j  ,k+lh,n) + (1.0-fracx)*flux4[2](i,j  ,k   ,n)) +
                                            fracy *(fracx*flux4[2](i+1,j+1,k+lh,n) + (1.0-fracx)*flux4[2](i,j+1,k+lh,n));
                                }
                            }

#endif
                        }


                        divc4(i,j,k,n) = 0.0;
                        for (int d=0; d < AMREX_SPACEDIM; ++d) {
                            index[d] = 1;


                            const Real& alpha_lo = afrac4[d](i, j, k);
                            const Real& flux_lo = flux[d][0];

                            const Real& alpha_hi = afrac4[d](i+index[0], j+index[1], k+index[2]);
                            const Real& flux_hi = flux[d][1];

                            const Real& wflux = wflux4[d](i,j,k,n);

                            divc4(i,j,k,n) -= (alpha_hi*flux_hi - alpha_lo*flux_lo)*dxinv[d]; // standard flux

                            divc4(i,j,k,n) -= (alpha_lo - alpha_hi)*wflux*dxinv[d]; // wall flux

                            index[d] = 0;
                        }

                        divc4(i,j,k,n) /= vfrac4(i,j,k);

                    }


                }


                // get the redistribution weights
                switch (redistribution_strategy) {
                case RedistributionEB::Uniform :
                    rediswgt4(i,j,k) = 1.0;
                    break;
                case RedistributionEB::VolumeFraction :
                    rediswgt4(i,j,k) = vfrac4(i,j,k);
                    break;
                case RedistributionEB::Density :
                    rediswgt4(i,j,k) = cons4(i,j,k,istate.get_cons_density_idx());
                    break;
                case RedistributionEB::Energy :
                    rediswgt4(i,j,k) = istate.get_energy_from_cons(istate.get_state_vector(cons,i,j,k));
                    break;
                }
            }
        }
    }

    FArrayBox optmp(grow(box,2));
    Array4<Real> const& optmp4 = optmp.array();

    FArrayBox delm(grow(box,1));
    Array4<Real> const& delm4 = delm.array();

    for (int n=0; n<nc; ++n) {

        //
        // compute mass loss from non-conservative flux
        //
        optmp.setVal(0.0);
        delm.setVal(0.0);
        Real vtot, divnc;
        for     (int k = lo.z-AMREX_D_PICK(0,0,1); k <= hi.z+AMREX_D_PICK(0,0,1); ++k) {
            for   (int j = lo.y-AMREX_D_PICK(0,1,1); j <= hi.y+AMREX_D_PICK(0,1,1); ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x-1; i <= hi.x+1; ++i) {

                    const EBCellFlag &cflag = flag4(i,j,k);
                    Real vf = vfrac4(i,j,k);

                    if ((vf > 0.0) && (vf < 1.0)) {
                        vtot = 0.0;
                        divnc = 0.0;

                        // compute divergence over cell as if there is no cut
                        divnc = 0.0;
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            index[d] = 1;
                            divnc += dxinv[d]*(flux4[d](i,j,k,n) - flux4[d](i+index[0], j+index[1], k+index[2], n));
                            index[d] = 0;
                        }

                        optmp4(i,j,k) = (1-vfrac4(i,j,k))*(divnc-divc4(i,j,k,n));
                        delm4(i,j,k) = -vfrac4(i,j,k)*optmp4(i,j,k);
                    }
                }
            }
        }

        //
        // redistribute
        //
        Real wtot;
        for     (int k = lo.z-AMREX_D_PICK(0,0,1); k <= hi.z+AMREX_D_PICK(0,0,1); ++k) {
            for   (int j = lo.y-AMREX_D_PICK(0,1,1); j <= hi.y+AMREX_D_PICK(0,1,1); ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x-1; i <= hi.x+1; ++i) {

                    const EBCellFlag &cflag = flag4(i,j,k);
                    Real vf = vfrac4(i,j,k);

                    if ((vf > 0.0) && (vf < 1.0)) {
                        wtot = 0.0;
                        for (const auto& g : grab) {
                            if (cflag.isConnected(g[0], g[1], g[2])) {
                                wtot += vfrac4(i+g[0],j+g[1],k+g[2])*rediswgt4(i+g[0],j+g[1],k+g[2]);
                            }
                        }

                        wtot = 1.0/(wtot + 1.e-80);

                        bool as_crse_crse_cell = false;
                        bool as_crse_covered_cell = false;
                        if (as_crse) {
                            as_crse_crse_cell = is_inside(i,j,k,lo,hi) && rr_flag_crse4(i,j,k) == CellType::crse_fine_boundary_cell;
                            as_crse_covered_cell = rr_flag_crse4(i,j,k) == CellType::fine_cell;
                        }

                        bool as_fine_valid_cell = false;  // valid cells near box boundary
                        bool as_fine_ghost_cell = false;  // ghost cells just outside valid region
                        if (as_fine) {
                            as_fine_valid_cell = is_inside(i,j,k,lo,hi);
                            as_fine_ghost_cell = levmsk4(i,j,k) == LevelMask::NotCovered; // not covered by other grids
                        }

                        for (const auto& g : grab) {
                            if (cflag.isConnected(g[0], g[1], g[2])) {

                                const int iii = i + g[0];
                                const int jjj = j + g[1];
                                const int kkk = k + g[2];

                                Real drho = delm4(i,j,k)*wtot*rediswgt4(iii,jjj,kkk);
                                optmp4(iii,jjj,kkk) += drho;

                                bool valid_dst_cell = is_inside(iii,jjj,kkk,lo,hi);

                                if (as_crse_crse_cell) {
                                    if (rr_flag_crse4(iii,jjj,kkk) == CellType::fine_cell && vfrac4(i,j,k) > reredistribution_threshold) {
                                        rr_drho_crse4(i,j,k,n) += drho*(vfrac4(iii,jjj,kkk)/vfrac4(i,j,k));
                                    }
                                }

                                if (as_crse_covered_cell) {
                                    if (valid_dst_cell) {
                                        if (rr_flag_crse4(iii,jjj,kkk) == CellType::crse_fine_boundary_cell && vfrac4(iii,jjj,kkk) > reredistribution_threshold) {
                                            // the recipient is a crse/fine boundary cell
                                            rr_drho_crse4(iii,jjj,kkk,n) -= drho;
                                        }
                                    }
                                }

                                if (as_fine_valid_cell) {
                                    if (!valid_dst_cell) {
                                        dm_as_fine4(iii,jjj,kkk,n) += drho*vfrac4(iii,jjj,kkk);
                                    }
                                }

                                if (as_fine_ghost_cell) {
                                    if (valid_dst_cell) {
                                        dm_as_fine4(i,j,k,n) -= drho*vfrac4(iii,jjj,kkk);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        //
        // push to output array
        //

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {
                    du4(i,j,k,n) = divc4(i,j,k,n) + optmp4(i,j,k);
                }
            }
        }
    }

    return;
}

std::string RedistributeEB::str() const
{
    std::stringstream msg;

    msg << get_tag() << "(";

    msg << "strategy=" << options[+redistribution_strategy];
    msg << ", reredistribution_threshold=" << reredistribution_threshold << ")";

    return msg.str();
}

//=============================================================================

std::string MergeEB::tag = "merge";
bool MergeEB::registered = GetDivergenceEBBuilder().Register(MergeEB::tag, DivergenceEBBuilder<MergeEB>);


MergeEB::MergeEB()
{
}

MergeEB::MergeEB(const sol::table &def)
{
    global_idx = def["global_idx"];

    const sol::table &div_def = def["eb_divergence"];

    // get the merge threshold

    merge_threshold = div_def.get_or("merge_threshold", 0.5);

}

MergeEB::~MergeEB()
{
}

void MergeEB::calc_eb_divergence(const Box& box,
                                 const FArrayBox &cons,
                                 Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                 Array<FArrayBox, AMREX_SPACEDIM> &wall_fluxes,
                                 FArrayBox& du,
                                 const EBCellFlagFab& flag,
                                 const FArrayBox& vfrac,
                                 const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                                 const Array<const FArrayBox *, AMREX_SPACEDIM> &fcent,
                                 int as_crse,
                                 int as_fine,
                                 const IArrayBox *rrflag_as_crse,
                                 const IArrayBox &levmsk, FArrayBox *rr_drho_crse,
                                 FArrayBox &dm_as_fine,
                                 const Real *dx,
                                 const Real dt) const
{
    BL_PROFILE("MergeEB::calc_eb_divergence");
    // make sure du is empty
    du.setVal(0.0);

    int N = du.nComp();

    const Dim3 lo = amrex::lbound(du.box());
    const Dim3 hi = amrex::ubound(du.box());

    Array<int,3> index = {0,0,0};
    Array4<Real> const& du4 = du.array();
    Real dxinv;

    Array4<const EBCellFlag> const& flag4 = flag.array();
    Array4<const Real> const& vfrac4 = vfrac.array();

    Array<Array4<const Real>,AMREX_SPACEDIM> wflux4;
    Array<Array4<const Real>,AMREX_SPACEDIM> flux4;
    Array<Array4<const Real>,AMREX_SPACEDIM> afrac4;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        flux4[d] = fluxes[d].array();
        wflux4[d] = wall_fluxes[d].array();
        afrac4[d] = afrac[d]->array();
    }

    Real wflux, vf;

    for (int n=0; n<N; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

                    const EBCellFlag &cflag = flag4(i,j,k);
                    vf = vfrac4(i,j,k);

                    if (vf <= 0.0) continue;

                    wflux = 0.0;

                    // cycle over dimensions
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        index[d] = 1;
                        dxinv = dt/dx[d];

                        if (cflag.isSingleValued()) {
                            wflux = wflux4[d](i,j,k,n);
                        }

                        Real flux_lo  = flux4[d](i,j,k,n);
                        const Real alpha_lo = afrac4[d](i,j,k);

                        Real flux_hi  = flux4[d](i+index[0], j+index[1], k+index[2], n);
                        const Real alpha_hi = afrac4[d](i+index[0], j+index[1], k+index[2]);

                        flux_lo *= alpha_lo;
                        flux_hi *= alpha_hi;
                        wflux *= (alpha_hi - alpha_lo);

                        du4(i,j,k,n) += dxinv*(flux_lo - flux_hi + wflux)/vf;
                        index[d] = 0;

                    }
                }
            }
        }
    }

    return;
}

void MergeEB::merge_cells(const Box& box,
                          FArrayBox &cons,
                          FArrayBox& du,
                          const EBCellFlagFab& flag,
                          const FArrayBox& vfrac,
                          const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                          int as_fine,
                          FArrayBox &dm_as_fine,
                          const IArrayBox& levmsk) const
{
    BL_PROFILE("MergeEB::merge_cells");
    int nc = du.nComp();

    Array<int,3> index = {0,0,0};

    Array4<Real> const& cons4 = cons.array();
    Array4<Real> const& du4 = du.array();

    Array4<const EBCellFlag> const& flag4 = flag.array();
    Array4<const Real> const& vfrac4 = vfrac.array();
    Array4<Real> const& dm_as_fine4 = dm_as_fine.array();
    Array4<const int> const& levmsk4 = levmsk.array();

    Array<Array4<const Real>,AMREX_SPACEDIM> afrac4;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        afrac4[d] = afrac[d]->array();
    }

    Real vf;

    // the sets of cells that are to be merged
    Vector<Vector<Array<int,3>>> merge;

    // a container for all of the unique cells that are part of the problem and
    // the super-cell that they are a part of
    std::map<Array<int,3>,int> cells;


    int ii, jj, kk;
    int mi, mj, mk;

    Real side_alpha, side_alpha_max;

    Box halo = grow(box,2);
    const Dim3 halo_lo = amrex::lbound(halo);
    const Dim3 halo_hi = amrex::ubound(halo);

    // expand zone over which we work so that each block sees the same modifications
    // to du without having to do inter-block communication
    Box calc = grow(box,1);
    const Dim3 lo = amrex::lbound(calc);
    const Dim3 hi = amrex::ubound(calc);


    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

                const EBCellFlag &cflag = flag4(i,j,k);

                vf = vfrac4(i,j,k);

                if (vf <= 0.0) continue;

                if (vf < merge_threshold) {

                    ii = i;
                    jj = j;
                    kk = k;

                    Array<int,3> cell_idx = {ii, jj, kk};

                    // check that this cell isn't already part of another super-cell
                    if (cells.find(cell_idx) != cells.end())
                        continue;

                    Vector<Array<int,3>> super_cell = {cell_idx};
                    bool new_super_cell = true;

                    while (vf < merge_threshold) {
                        // get the fluid fraction of each of the sides
                        // scaled by if it is a viable candidate
                        side_alpha_max = 0.0;
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            index[d] = 1;

                            // lo side

                            if (is_inside(ii-index[0], jj-index[1], kk-index[2], halo_lo, halo_hi)) {
                                side_alpha = afrac4[d](ii,jj,kk);
                                if (side_alpha > side_alpha_max) {
                                    side_alpha_max = side_alpha;
                                    mi = ii-index[0];
                                    mj = jj-index[1];
                                    mk = kk-index[2];

                                }
                            }

                            // hi side

                            if (is_inside(ii+index[0], jj+index[1], kk+index[2], halo_lo, halo_hi)) {
                                side_alpha = afrac4[d](ii+index[0], jj+index[1], kk+index[2]);
                                if (side_alpha > side_alpha_max) {
                                    side_alpha_max = side_alpha;
                                    mi = ii+index[0];
                                    mj = jj+index[1];
                                    mk = kk+index[2];

                                }
                            }

                            index[d] = 0;
                        }

                        ii = mi;
                        jj = mj;
                        kk = mk;

                        Array<int,3> cell_idx = {ii, jj, kk};

                        // check if this cell is part of another super cell
                        if (cells.find(cell_idx) == cells.end()) {
                            vf += vfrac4(ii,jj,kk);
                            super_cell.push_back(cell_idx);
                        } else {
                            int si = cells[cell_idx];

                            // register the cells with a super cell
                            for (const auto &cell : super_cell) {
                                merge[si].push_back(cell);
                                cells[cell] = si;
                            }
                            new_super_cell = false;
                            break;
                        }
                    }

                    // add the merged cell to the list
                    if (new_super_cell) {

                        for (const auto &cell : super_cell) {
                            cells[cell] = merge.size();
                        }

                        merge.push_back(super_cell);
                    }



                } else if (cflag.isSingleValued() && as_fine) {
                    for (int n=0; n<nc; ++n) {
                        dm_as_fine4(i,j,k,n) = 0.0;
                    }
                }
            }
        }
    }

    int nsuper = merge.size();

    // do update

    for (int n=0; n<nc; ++n) {

        for (int mi=0; mi<nsuper; ++mi) {
            const auto& mc = merge[mi];


            // compute the sums
            Real vf_sum   = 0.0;
            Real sum_cons = 0.0;
            Real sum_du   = 0.0;

            for (const auto& cell_idx : mc) {
                vf = vfrac4(cell_idx[0],cell_idx[1], cell_idx[2]);
                vf_sum += vf;

                sum_cons += cons4(cell_idx[0],cell_idx[1], cell_idx[2],n)*vf;
                sum_du += du4(cell_idx[0],cell_idx[1], cell_idx[2],n)*vf;
            }

            sum_cons /= vf_sum;
            sum_du /= vf_sum;

            for (const auto& cell_idx : mc) {

                const int i = cell_idx[0];
                const int j = cell_idx[1];
                const int k = cell_idx[2];

                const Real u_orig = cons4(i,j,k,n);
                const Real u_final = sum_cons + sum_du;
                const Real du_merge = u_final - u_orig;

                du4(i,j,k,n) = du_merge;

                // tell other grids how du has changed
                if (as_fine && (levmsk4(i,j,k) == LevelMask::NotCovered)) {
                    dm_as_fine4(i,j,k,n) = du_merge - du4(cell_idx[0],cell_idx[1], cell_idx[2],n);
                }
            }
        }

    }
}

std::string MergeEB::str() const
{
    std::stringstream msg;

    msg << get_tag() << "(";

    msg << "merge_threshold=" << merge_threshold << ")";

    return msg.str();
}

#endif
