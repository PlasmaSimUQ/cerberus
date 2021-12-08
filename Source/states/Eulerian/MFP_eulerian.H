#ifndef EULERIANSTATE_H
#define EULERIANSTATE_H

#include "MFP.H"
#include "MFP_state.H"
#include "MFP_reconstruction.H"

enum LevelMask {
    Inside = 0, // valid cells
    Covered,    // ghost cells covered by valid cells of this level
    NotCovered, // ghost cells not covered
    PhysBnd,    // outside domain
    NUM
};

class EulerianState : public State
{
public:
    EulerianState();
    ~EulerianState();

    int num_grow;
    bool reflux;
    Real effective_zero;

#ifdef AMREX_USE_EB
    Real merge_threshold;
#endif

    std::map<std::string,Optional3D1VFunction> other_functions;
    Vector<Optional3D1VFunction> functions;

    std::unique_ptr<Reconstruction> reconstructor;

    virtual StateClassification get_classification() const override {return StateClassification::Eulerian;}
    virtual int get_num_cons() const {return 0;}
    virtual const Vector<std::string>& get_cons_names() const = 0;

    virtual void get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string,FArrayBox>& out,
                                 Vector<std::string>& updated
                             #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                             #endif
                                 ) const {return;}

    virtual void set_reflux();
    virtual void set_reconstruction();
    virtual int get_num_grow() const = 0;
    virtual void set_num_grow(int n) = 0;


    virtual void init_from_lua() override;

    virtual bool is_transported() const {return true;}
    virtual bool is_viscous() const {return false;}

    virtual void calc_primitives(const Box& box,
                                 FArrayBox& cons,
                                 FArrayBox &prim,
                                 const Real* dx,
                                 const Real t,
                                 const Real* prob_lo
                             #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                             #endif
                                 ) const {}

    virtual void update_boundary_cells(const Box& box,
                                       const Geometry &geom,
                                       FArrayBox &prim,
                                   #ifdef AMREX_USE_EB
                                       const FArrayBox& vfrac,
                                   #endif
                                       const Real time) const {}

    virtual void calc_reconstruction(const Box& box,
                                     FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                 #ifdef AMREX_USE_EB
                                     ,const EBCellFlagFab &flag
                                     ,const FArrayBox &vfrac
                                 #endif
                                     ) const {}

    virtual void calc_time_averaged_faces(const Box& box,
                                          const FArrayBox &prim,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                                      #ifdef AMREX_USE_EB
                                          const EBCellFlagFab& flag,
                                      #endif
                                          const Real* dx,
                                          Real dt) const {}

    virtual void face_bc(const int dir,
                         Box const& box,
                         const FArrayBox& src,
                         FArrayBox& dest,
                         const Geometry &geom,
                     #ifdef AMREX_USE_EB
                         const EBCellFlagFab& flag,
                     #endif
                         const Real time,
                         const bool do_all) const {}

    virtual void update_face_prim(const Box& box, const Geometry& geom,
                                  Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                                  Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                              #ifdef AMREX_USE_EB
                                  const EBCellFlagFab& flag,
                              #endif
                                  const Real time, const bool do_all = false) const;

    virtual void calc_fluxes(const Box& box,
                             FArrayBox &cons,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                             Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                         #ifdef AMREX_USE_EB
                             const EBCellFlagFab& flag,
                         #endif
                             const Real *dx,
                             const Real dt) const {}

    virtual void correct_face_prim(const Box& box,
                                   Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                                   Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                                   const Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                               #ifdef AMREX_USE_EB
                                   const EBCellFlagFab &flag,
                               #endif
                                   const Real *dx,
                                   const Real dt) const {}

    virtual void calc_viscous_fluxes(const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     const FArrayBox &prim,
                                 #ifdef AMREX_USE_EB
                                     const EBCellFlagFab& flag,
                                 #endif
                                     const Real* dx) const {}

#ifdef AMREX_USE_EB
    static Real interp2d(Real cym,
                         Real cy0,
                         Real cyp,
                         Real czm,
                         Real cz0,
                         Real czp,
                         Eigen::Matrix3d v);

    // use biquadratic interpolation for 2D and 3D (inefficient but reduces code, implement more efficient
    // interpolation for 2D later).
    // wall values are expected in global coordinate system
    template<size_t N>
    static Array<Real,N> calc_wall_normal_slopes(const Array<int,N> &slope_idx,
                                                 const Array<Real,N> &wall_value,
                                                 Array<Real,3> &wall_normal,
                                                 Array<Real, AMREX_SPACEDIM> &wall_centre,
                                                 Array4<const Real> const &prim4,
                                                 const int i, const int j, const int k)
    {
        BL_PROFILE("HydroState::calc_wall_normal_slopes");
        Array<int,AMREX_SPACEDIM> index;
        index.fill(0);

        Array<Real,N> dat1, dat2;
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
            for (int di = 0; di<N; ++di) {
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
            for (int di = 0; di<N; ++di) {
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

            for (int di = 0; di<N; ++di) {
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

            for (int di = 0; di<N; ++di) {
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

            for (int di = 0; di<N; ++di) {
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

            for (int di = 0; di<N; ++di) {
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
        Array<Real,N> ddat;
        std::fill(ddat.begin(), ddat.end(), 0.0);

        for (int d=0; d<N; ++d) {
            ddat[d] = -ddinv*(d2*d2*(dat1[d] - wall_value[d])-d1*d1*(dat2[d] - wall_value[d]));
        }

        return ddat;
    }

    virtual void calc_wall_fluxes(const Box& box,
                                  const FArrayBox &prim,
                                  Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                  const EBCellFlagFab& flag,
                                  const CutFab &bc_idx,
                                  const FArrayBox& bcent,
                                  const FArrayBox &bnorm,
                                  const Array<const FArrayBox*, AMREX_SPACEDIM> &afrac,
                                  const Real *dx,
                                  const Real dt) const {}

    void calc_eb_divergence(const Box& box,
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
                            const Real dt) const;

    void merge_cells(const Box& box,
                     FArrayBox &cons,
                     FArrayBox& du,
                     const EBCellFlagFab& flag,
                     const FArrayBox& vfrac,
                     const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                     int as_fine,
                     FArrayBox &dm_as_fine,
                     const IArrayBox& levmsk) const;

    virtual void get_wall_value(const Box& box,
                                Vector<FArrayBox*> bcs_data,
                                const EBCellFlagFab& flag,
                                const CutFab &bc_idx,
                                const FArrayBox& bcent,
                                const FArrayBox &bnorm,
                                const Real t,
                                const Real* dx,
                                const Real* prob_lo) const {}

    static bool check_covered_stencil(Array4<const EBCellFlag> const& flag, int i, int j, int k, int d, int stencil_length);
#endif

    void calc_divergence(const Box& box,
                         Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                         FArrayBox& du,
                         const Real *dx,
                         const Real dt) const;

    static void calc_slope(const Box& box,
                           const FArrayBox& src,
                           FArrayBox &slope,
                       #ifdef AMREX_USE_EB
                           const EBCellFlagFab &flag,
                       #endif
                           const Real *dx,
                           const int src_idx,
                           const int slope_idx,
                           const int dim,
                           Reconstruction &reco);

    virtual void write_info(nlohmann::json& js) const override;

    static EulerianState& get_state_global(const int idx) {return static_cast<EulerianState&>(MFP::get_state(idx));}
    static EulerianState& get_state(const int idx) {return static_cast<EulerianState&>(MFP::get_state(MFP::eulerian_states[idx]));}
    static EulerianState& get_state(const std::string& name) {return static_cast<EulerianState&>(MFP::get_state(name));}

};

#endif // EULERIANSTATE_H
