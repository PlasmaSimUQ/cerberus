#include "MFP_state.H"

#include <math.h>
#include "MFP_global.H"
#include "MFP_Riemann_solvers.H"
#include "MFP_diagnostics.H"
#include "MFP_viscous.H"
#include "MFP_eb_bc.H"


using GD = GlobalData;

#ifdef AMREX_USE_EB

void State::set_eb_bc(const sol::table &bc_def)
{
    return;
}

void State::set_eb_divergence()
{

    PhysicsFactory<DivergenceEB> div_fact = GetDivergenceEBBuilder();

    sol::table state_def = GD::lua["states"][name];
    state_def["global_idx"] = global_idx;

    sol::optional<sol::table> div_def_exists = state_def["eb_divergence"];

    // set a default
    if (!div_def_exists) {
        state_def["eb_divergence"] = GD::lua["eb_divergence_default"];
    }

    sol::table div_def = state_def["eb_divergence"];

    std::string tag = div_def.get<std::string>("type");
    eb_div = div_fact.Build(tag, state_def);

    if (!eb_div)
        Abort("Invalid 'eb_divergence' type '"+tag+"'. Options are "+vec2str(div_fact.getKeys()));

    return;
}

bool State::check_covered_stencil(Array4<const EBCellFlag> const& flag, int i, int j, int k, int d, int stencil_length)
{
    BL_PROFILE("State::check_covered_stencil");
    Array<int,3> stencil_index;
    int offset = stencil_length/2;
    // cell that references a covered cell doesn't need calculating
    stencil_index.fill(0);
    for (int s=0; s<stencil_length; ++s) {
        stencil_index[d] = s - offset;
        // check if any of the stencil values are from a covered cell
        if (flag(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2]).isCovered()) {
            return true;
        }
    }
    return false;
}

Real interp2d(Real cym,
              Real cy0,
              Real cyp,
              Real czm,
              Real cz0,
              Real czp,
              Eigen::Matrix3d v)
{
    return czm*(cym*v(0,0) + cy0*v(1,0) + cyp*v(2,0))
            +     cz0*(cym*v(0,1) + cy0*v(1,1) + cyp*v(2,1))
            +     czp*(cym*v(0,2) + cy0*v(1,2) + cyp*v(2,2));
}

// The wall is assumed to be adiabatic, thus dTdn=0.
// We use no-slip boundary for velocities.
// use biquadratic interpolation for 2D and 3D (inefficient but reduces code, implement more efficient
// interpolation for 2D later).
// wall values are expected in global coordinate system
Vector<Real> State::calc_wall_normal_slopes(const Vector<int> &slope_idx,
                                            const Vector<Real> &wall_value,
                                            Array<Real,3> &wall_normal,
                                            Array<Real,AMREX_SPACEDIM> &wall_centre,
                                            Array4<const Real> const &prim4,
                                            const int i, const int j, const int k) const
{
    BL_PROFILE("State::calc_wall_normal_slopes");
    Array<int,AMREX_SPACEDIM> index;
    index.fill(0);

    int nslopes = slope_idx.size();

    Vector<Real> dat1(nslopes), dat2(nslopes);
    Eigen::Matrix3d dat;
    Real d1, d2;

    if (std::abs(wall_normal[0]) >= std::abs(wall_normal[1]) && std::abs(wall_normal[0])  >= std::abs(wall_normal[2])) {
        // y-z plane: x = const
        // the equation for the line:  x = wall_centre[0] - d*wall_normalx
        //                             y = wall_centre[1] - d*wall_normaly
        //                             z = wall_centre[2] - d*wall_normalz

        Real s = sign(1.0, -wall_normal[0]);
        int is = (int)nearbyint(s);

        //
        // the line intersects the y-z plane (x = s) at ...
        //
        d1 = (wall_centre[0] - s) * (1.00/wall_normal[0]);  // this is also the distance from wall to intersection
        Real yit = wall_centre[1] - d1*wall_normal[1];
        Real zit = wall_centre[2] - d1*wall_normal[2];
        int iyit = j + (int)nearbyint(yit);
        int izit = k + (int)nearbyint(zit);
        yit = yit - nearbyint(yit);  // shift so that the center of the nine cells are (0.,0.)
        zit = zit - nearbyint(zit);
        //
        // coefficents for quadratic interpolation
        Real cym = 0.5*yit*(yit-1.0);
        Real cy0 = 1.0-yit*yit;
        Real cyp = 0.50*yit*(yit+1.0);
        Real czm = 0.50*zit*(zit-1.0);
        Real cz0 = 1.0-zit*zit;
        Real czp = 0.50*zit*(zit+1.0);

        //
        // interploation
        for (int di = 0; di<nslopes; ++di) {
            int d = slope_idx[di];

            // grab a 3x3 stencil
            for (int ii=-1; ii<2; ++ii) {
#if AMREX_SPACEDIM == 2
                dat(ii+1,0) = prim4(i+is,j+ii,k,d);
                dat(ii+1,1) = dat(ii+1,0);
                dat(ii+1,2) = dat(ii+1,0);
#else
                for (int jj=-1; jj<2; ++jj) {
                    dat(ii+1,jj+1) = prim4(i+is,j+ii,k+jj,d);
                }
#endif
            }

            // interpolate the velocity
            dat1[di] = interp2d(cym,cy0,cyp,czm,cz0,czp,dat);
        }

        //
        // the line intersects the y-z plane (x = 2*s) at ...
        //
        d2 = (wall_centre[0] - 2.0*s) * (1.0/wall_normal[0]);
        yit = wall_centre[1] - d2*wall_normal[1];
        zit = wall_centre[2] - d2*wall_normal[2];
        iyit = j + (int)nearbyint(yit);
        izit = k + (int)nearbyint(zit);
        yit = yit - nearbyint(yit);  // shift so that the center of the nine cells are (0.,0.)
        zit = zit - nearbyint(zit);
        //
        // coefficents for quadratic interpolation
        cym = 0.5*yit*(yit-1.0);
        cy0 = 1.0-yit*yit;
        cyp = 0.5*yit*(yit+1.0);
        czm = 0.5*zit*(zit-1.0);
        cz0 = 1.0-zit*zit;
        czp = 0.5*zit*(zit+1.0);
        //
        // interploation
        for (int di = 0; di<nslopes; ++di) {
            int d = slope_idx[di];

            // grab a 3x3 stencil
            for (int ii=-1; ii<2; ++ii) {
#if AMREX_SPACEDIM == 2
                dat(ii+1,0) = prim4(i+2*is,j+ii,k,d);
                dat(ii+1,1) = dat(ii+1,0);
                dat(ii+1,2) = dat(ii+1,0);
#else
                for (int jj=-1; jj<2; ++jj) {
                    dat(ii+1,jj+1) = prim4(i+2*is,j+ii,k+jj,d);
                }
#endif
            }

            // interpolate the velocity
            dat2[di] = interp2d(cym,cy0,cyp,czm,cz0,czp,dat);
        }
    } else if (std::abs(wall_normal[1]) >= std::abs(wall_normal[0]) && std::abs(wall_normal[1]) >= std::abs(wall_normal[2])) {
        // z-x plane
        Real s = sign(1.0, -wall_normal[1]);
        int is = (int)nearbyint(s);

        d1 = (wall_centre[1] - s) * (1.0/wall_normal[1]);
        Real xit = wall_centre[0] - d1*wall_normal[0];
        Real zit = wall_centre[2] - d1*wall_normal[2];
        int ixit = i + (int)nearbyint(xit);
        int izit = k + (int)nearbyint(zit);
        xit = xit - nearbyint(xit);
        zit = zit - nearbyint(zit);

        Real cxm = 0.5*xit*(xit-1.0);
        Real cx0 = 1.0-xit*xit;
        Real cxp = 0.5*xit*(xit+1.0);
        Real czm = 0.5*zit*(zit-1.0);
        Real cz0 = 1.0-zit*zit;
        Real czp = 0.5*zit*(zit+1.0);

        for (int di = 0; di<nslopes; ++di) {
            int d = slope_idx[di];

            // grab a 3x3 stencil
            for (int ii=-1; ii<2; ++ii) {
#if AMREX_SPACEDIM == 2
                dat(ii+1,0) = prim4(i+ii,j+is,k,d);
                dat(ii+1,1) = dat(ii+1,0);
                dat(ii+1,2) = dat(ii+1,0);
#else
                for (int jj=-1; jj<2; ++jj) {
                    dat(ii+1,jj+1) = prim4(i+ii,j+is,k+jj,d);
                }
#endif
            }

            // interpolate the velocity
            dat1[di] = interp2d(cxm,cx0,cxp,czm,cz0,czp,dat);
        }

        d2 = (wall_centre[1] - 2.0*s) * (1.0/wall_normal[1]);
        xit = wall_centre[0] - d2*wall_normal[0];
        zit = wall_centre[2] - d2*wall_normal[2];
        ixit = i + (int)nearbyint(xit);
        izit = k + (int)nearbyint(zit);
        xit = xit - nearbyint(xit);
        zit = zit - nearbyint(zit);

        cxm = 0.5*xit*(xit-1.0);
        cx0 = 1.0-xit*xit;
        cxp = 0.5*xit*(xit+1.0);
        czm = 0.5*zit*(zit-1.0);
        cz0 = 1.0-zit*zit;
        czp = 0.5*zit*(zit+1.0);

        for (int di = 0; di<nslopes; ++di) {
            int d = slope_idx[di];

            // grab a 3x3 stencil
            for (int ii=-1; ii<2; ++ii) {
#if AMREX_SPACEDIM == 2
                dat(ii+1,0) = prim4(i+ii,j+2*is,k,d);
                dat(ii+1,1) = dat(ii+1,0);
                dat(ii+1,2) = dat(ii+1,0);
#else
                for (int jj=-1; jj<2; ++jj) {
                    dat(ii+1,jj+1) = prim4(i+ii,j+2*is,k+jj,d);
                }
#endif
            }

            // interpolate the velocity
            dat2[di] = interp2d(cxm,cx0,cxp,czm,cz0,czp,dat);
        }

    }

#if AMREX_SPACEDIM == 3

    else {
        // x-y plane
        Real s = sign(1.0,-wall_normal[2]);
        int is = (int)nearbyint(s);

        d1 = (wall_centre[2] - s) * (1.0/wall_normal[2]);
        Real xit = wall_centre[0] - d1*wall_normal[0];
        Real yit = wall_centre[1] - d1*wall_normal[1];
        int ixit = i + (int)nearbyint(xit);
        int iyit = j + (int)nearbyint(yit);
        xit = xit - nearbyint(xit);
        yit = yit - nearbyint(yit);

        Real cxm = 0.5*xit*(xit-1.0);
        Real cx0 = 1.0-xit*xit;
        Real cxp = 0.5*xit*(xit+1.0);
        Real cym = 0.5*yit*(yit-1.0);
        Real cy0 = 1.0-yit*yit;
        Real cyp = 0.5*yit*(yit+1.0);

        for (int di = 0; di<nslopes; ++di) {
            int d = slope_idx[di];

            // grab a 3x3 stencil
            for (int ii=-1; ii<2; ++ii) {
                for (int jj=-1; jj<2; ++jj) {
                    dat(ii+1,jj+1) = prim4(i+ii,j+jj,k+is,d);
                }
            }

            // interpolate the velocity
            dat1[di] = interp2d(cxm,cx0,cxp,cym,cy0,cyp,dat);
        }

        d2 = (wall_centre[2] - 2.0*s) * (1.0/wall_normal[2]);
        xit = wall_centre[0] - d2*wall_normal[0];
        yit = wall_centre[1] - d2*wall_normal[0];
        ixit = i + (int)nearbyint(xit);
        iyit = j + (int)nearbyint(yit);
        xit = xit - nearbyint(xit);
        yit = yit - nearbyint(yit);

        cxm = 0.5*xit*(xit-1.0);
        cx0 = 1.0-xit*xit;
        cxp = 0.5*xit*(xit+1.0);
        cym = 0.5*yit*(yit-1.0);
        cy0 = 1.0-yit*yit;
        cyp = 0.5*yit*(yit+1.0);

        for (int di = 0; di<nslopes; ++di) {
            int d = slope_idx[di];

            // grab a 3x3 stencil
            for (int ii=-1; ii<2; ++ii) {
                for (int jj=-1; jj<2; ++jj) {
                    dat(ii+1,jj+1) = prim4(i+ii,j+jj,k+2*is,d);
                }
            }

            // interpolate the velocity
            dat2[di] = interp2d(cxm,cx0,cxp,cym,cy0,cyp,dat);
        }

    }
#endif

    //
    // compute derivatives on the wall given a value on the wall.
    //
    Real ddinv = 1.0/(d1*d2*(d2-d1));
    Vector<Real> ddat(nslopes, 0.0);

    for (int d=0; d<nslopes; ++d) {
        ddat[d] = -ddinv*(d2*d2*(dat1[d] - wall_value[d])-d1*d1*(dat2[d] - wall_value[d]));
    }

    return ddat;
}

void State::calc_wall_fluxes(const Box& box,
                             const Vector<FArrayBox> &prim,
                             Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                             const EBCellFlagFab& flag,
                             const CutFab &bc_idx,
                             const FArrayBox& bcent,
                             const FArrayBox &bnorm,
                             const Array<const FArrayBox*, AMREX_SPACEDIM> &afrac,
                             const Real *dx,
                             const Real dt) const
{
    BL_PROFILE("State::calc_wall_fluxes");
    int nc = n_cons();
    int np = n_prim();

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Vector<Array4<const Real>> prim4(prim.size());
    for (size_t idx=0; idx<prim.size(); ++idx) {
        prim4[idx] = prim[idx].array();
    }

    Array4<const Real> const& p4 = prim4[global_idx];
    Array4<const Real> const& bcent4 = bcent.array();
    Array4<const Real> const& bnorm4 = bnorm.array();

    Array<Array4<Real>,AMREX_SPACEDIM> flux4;
    Array<Array4<const Real>,AMREX_SPACEDIM> afrac4;

    Array<Vector<Real>,AMREX_SPACEDIM> wall_flux;

    Array4<const EBCellFlag> const& flag4 = flag.array();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        flux4[d] = fluxes[d].array();
        afrac4[d] = afrac[d]->array();

        // zero out the flux accumulator
        fluxes[d].setVal(0.0);

        wall_flux[d].resize(nc);
    }

    const Array4<const Real>& bc_idx4 = bc_idx.array();

//    plot_FAB_2d(bc_idx, 0, "bc idx", false, true);

    Vector<Real> cell_state(np);

    Array<Array<Real,3>,3> wall_coord = {{{0,0,0},{0,0,0},{0,0,0}}};
    Array<Real,AMREX_SPACEDIM> wall_centre;

    for (int k = lo.z-AMREX_D_PICK(0,0,2); k <= hi.z+AMREX_D_PICK(0,0,2); ++k) {
            for (int j = lo.y-AMREX_D_PICK(0,2,2); j <= hi.y+AMREX_D_PICK(0,2,2); ++j) {
                AMREX_PRAGMA_SIMD
                    for (int i = lo.x-2; i <= hi.x+2; ++i) {

                const EBCellFlag &cflag = flag4(i,j,k);

                if (cflag.isSingleValued()) {


                    // grab a vector of the local state
                    cell_state = load_state_for_flux(prim4, i, j, k);

                    for (int d=0; d<AMREX_SPACEDIM; ++d) {

                        // get the wall normal
                        wall_coord[0][d] = bnorm4(i,j,k,d);

                        // get the centre of the wall
                        wall_centre[d] = bcent4(i,j,k,d);
                    }

                    // get a local coordinate system with x- aligned with the wall normal
                    expand_coord(wall_coord);

                    // the boundary condition
                    const int ebi = (int)nearbyint(bc_idx4(i,j,k));
                    const auto& bc = eb_bcs[ebi];

                    // if necessary get the normal slope
                    Vector<Real> slopes;
                    if (bc->need_slopes()) {
                        Array<Real,3> &wall_normal = wall_coord[0];
                        Vector<Real> wall_value = bc->local2global(wall_coord);
                        slopes = calc_wall_normal_slopes(bc->get_slopes(),
                                                         wall_value,
                                                         wall_normal,
                                                         wall_centre,
                                                         p4,
                                                         i, j, k);
                    }

                    // calculate the wall flux
                    bc->solve(wall_coord, cell_state, slopes, wall_flux, dx);

                    // load the flux into the fab
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        for (int n=0; n<nc; ++n) {
                            flux4[d](i,j,k,n) += wall_flux[d][n];
                        }
                    }

                }
            }
        }
    }
    return;
}


void State::get_wall_value(const Box& box,
                           std::map<BoundaryEB::EBType,FArrayBox*> bcs_data,
                           const EBCellFlagFab& flag,
                           const CutFab &bc_idx,
                           const FArrayBox& bcent,
                           const FArrayBox &bnorm,
                           const Real t,
                           const Real* dx,
                           const Real* prob_lo) const
{
    BL_PROFILE("State::get_wall_value");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& bcent4 = bcent.array();
    Array4<const Real> const& bnorm4 = bnorm.array();
    Array4<const EBCellFlag> const& flag4 = flag.array();
    const Array4<const Real>& bc_idx4 = bc_idx.array();

    std::map<BoundaryEB::EBType, Vector<Real>> wall_state;

    Array<Array<Real,3>,3> wall_coord = {{{0,0,0},{0,0,0},{0,0,0}}};
    Array<Real,AMREX_SPACEDIM> wall_centre;

    Array<Real,3> xyz;

    for     (int k = lo.z; k <= hi.z; ++k) {
        xyz[2] = prob_lo[2] + (k + 0.5)*dx[2];
        for   (int j = lo.y; j <= hi.y; ++j) {
            xyz[1] = prob_lo[1] + (j + 0.5)*dx[1];
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                xyz[0] = prob_lo[0] + (i + 0.5)*dx[0];

                const EBCellFlag &cflag = flag4(i,j,k);

                if (cflag.isSingleValued()) {

                    // the boundary condition
                    const int ebi = (int)nearbyint(bc_idx4(i,j,k));
                    const auto& bc = eb_bcs[ebi];


                    for (int d=0; d<AMREX_SPACEDIM; ++d) {

                        // get the wall normal
                        wall_coord[0][d] = bnorm4(i,j,k,d);

                        // get the centre of the wall (in cell coordinates [0,1])
                        // and convert to x,y,z coordinates
                        wall_centre[d] = xyz[d] + bcent4(i,j,k,d)*dx[d];
                    }

                    // get a local coordinate system with x- aligned with the wall normal
                    expand_coord(wall_coord);

                    // calculate the wall state
                    wall_state = bc->get_wall_state(wall_centre, wall_coord, t);

                    for (const auto &ws : wall_state) {

                        // only proceed if the boundary condition type is in the list
                        if ( bcs_data.find(ws.first) == bcs_data.end() ) {
                          continue;
                        }

                        Array4<Real> const& wall4 = bcs_data[ws.first]->array();

                        // load the wall value into the fab
                        for (size_t n = 0; n<ws.second.size(); ++n) {
                            wall4(i,j,k,n) = ws.second[n];
                        }

                    }
                }
            }
        }
    }
    return;
}

#endif
