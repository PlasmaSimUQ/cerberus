#include "MFP_eulerian.H"
#include "MFP.H"
#include "MFP_transforms.H"

EulerianState::EulerianState()
{
    num_grow = 0;
}

EulerianState::~EulerianState(){}

void EulerianState::set_reconstruction()
{

    ClassFactory<Reconstruction> rfact = GetReconstructionFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string rec = state_def["reconstruction"].get_or<std::string>("null");

    state_def["reconstruction"] = rec; // consistency when using default option "null"

    if (is_transported() && rec == "null")
        Abort("Reconstruction option required for state '"+name+"'. Options are "+vec2str(rfact.getKeys()));

    reconstructor = rfact.Build(rec, state_def);

    if (!reconstructor)
        Abort("Invalid reconstruction option '"+rec+"'. Options are "+vec2str(rfact.getKeys()));

    int ng = reconstructor->get_num_grow();
    set_num_grow(ng);

    return;

}

void EulerianState::set_reflux()
{
    if (AMREX_SPACEDIM < 2) {
        reflux = 0;
    } else {
        reflux = MFP::lua["states"][name]["reflux"].get_or(1);
    }
}

void EulerianState::init_from_lua()
{
    BL_PROFILE("EulerianState::init_from_lua");

    sol::state& lua = MFP::lua;
    const sol::table state_def = lua["states"][name];

    effective_zero = state_def.get_or("effective_zero", 1e-14);

    //
    // reflux
    //
    set_reflux();

    //
    // get reconstruction
    //

    set_reconstruction();

#ifdef AMREX_USE_EB
    merge_threshold = state_def.get_or("merge_threshold", 0.5);
#endif

}


void EulerianState::init_data(MFP* mfp, const Real time)
{

    const Real* dx = mfp->Geom().CellSize();
    const Real* prob_lo = mfp->Geom().ProbLo();

    MultiFab& S_new = mfp->get_data(data_idx, time);

    Vector<Real> U(n_cons());
    Vector<Real> Q(n_prim());

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        FArrayBox& S_data = S_new[mfi];

        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);
        Array4<Real> const& S4 = S_data.array();

#ifdef AMREX_USE_EB
        FabArray<EBCellFlagFab>& flags = mfp->get_eb_data(global_idx).flags;
        Array4<const EBCellFlag> const& flag4 = flags.array(mfi);
#endif

        Real x, y, z;
        for     (int k = lo.z; k <= hi.z; ++k) {
            z = prob_lo[2] + (k + 0.5)*dx[2];
            for   (int j = lo.y; j <= hi.y; ++j) {
                y = prob_lo[1] + (j + 0.5)*dx[1];
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {
                    x = prob_lo[0] + (i + 0.5)*dx[0];

#ifdef AMREX_USE_EB
                    const EBCellFlag &cflag = flag4(i,j,k);

                    if (cflag.isCovered()) {
                        for (int n=0; n<n_cons(); ++n) {
                            S4(i,j,k,n) = 0.0;
                        }
                        continue;
                    }
#endif

                    // grab the primitive variables as defined by the user functions
                    for (int icomp=0; icomp<n_prim(); ++icomp) {
                        const auto& f = functions[icomp];

                        Q[icomp] = f(x, y, z);

                    }

                    // convert primitive to conserved
                    prim2cons(Q, U);

                    // copy into array
                    for (int n=0; n<n_cons(); ++n) {
                        S4(i,j,k,n) = U[n];
                    }
                }
            }
        }
    }
}


Vector<std::string> EulerianState::get_plot_output_names() const
{
    Vector<std::string> out;
    const Vector<std::string>& cons_names = get_cons_names();
    out.insert(out.end(), cons_names.begin(), cons_names.end());

    const Vector<std::string>& prim_names = get_prim_names();
    out.insert(out.end(), prim_names.begin(), prim_names.end());
#ifdef AMREX_USE_EB
    out.push_back("vfrac");
#endif

    return out;
}

void EulerianState::get_plot_output(const Box& box,
                                 const FArrayBox& src,
                                 std::map<std::string,FArrayBox>& out,
                                 Vector<std::string>& updated
                                 #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                                 #endif
                                 ) const
{
    BL_PROFILE("MHDState::get_state_values");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
    Array4<const Real> const& vf4 = vfrac.array();
#endif

    updated.resize(0);

    const Vector<std::string>& cons_names = get_cons_names();
    const Vector<std::string>& prim_names = get_prim_names();

    // check conserved variables
    std::map<std::string,int> cons_tags;
    for (int i=0; i<n_cons(); ++i) {
        const std::string s = cons_names[i];
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        cons_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // check primitive variables
    std::map<std::string,int> prim_tags;
    for (int i=0; i<n_prim(); ++i) {
        const std::string s = prim_names[i];
        if (s == cons_names[0]) continue;
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        prim_tags[var_name] = i;
        updated.push_back(var_name);
    }

#ifdef AMREX_USE_EB
    const std::string vfrac_name = "vfrac-"+name;
    bool load_vfrac = out.find(vfrac_name) != out.end();
#endif

    std::map<std::string,Array4<Real>> out4;
    for (const std::string& s : updated) {
        out[s].resize(box, 1);
        out[s].setVal(0.0);
        out4[s] = out[s].array();
    }

    // temporary storage for retrieving the state data
    Vector<Real> S(n_cons());
    Vector<Real> Q(n_prim());

    Array4<const Real> const& src4 = src.array();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (vf4(i,j,k) == 0.0) {
                    continue;
                }
#endif

                for (int n=0; n<n_cons(); ++n) {
                    S[n] = src4(i,j,k,n);

                    if (std::isnan(S[n])) {
                        Abort();
                    }

                }


#ifdef AMREX_USE_EB
                if (load_vfrac)  out4[vfrac_name](i,j,k)  = vf4(i,j,k);
#endif

                if (!cons_tags.empty()) {
                    for (const auto& var : cons_tags) {
                        out4[var.first](i,j,k) = S[var.second];
                    }
                }

                if (!prim_tags.empty()) {
                    cons2prim(S, Q);

                    for (const auto& var : prim_tags) {
                        out4[var.first](i,j,k) = Q[var.second];
                    }
                }
            }
        }
    }


    return;
}




Real EulerianState::get_allowed_time_step(MFP* mfp) const
{
    const Real* dx = mfp->Geom().CellSize();

    MultiFab& data = mfp->get_new_data(data_idx);

    Vector<Real> U(n_cons());

    Real max_speed = std::numeric_limits<Real>::min();

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls
        const FArrayBox& vfrac = mfp->get_eb_data(global_idx).volfrac[mfi];

        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();
#endif
        Array4<const Real> const& data4 = data[mfi].array();

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) {
                        continue;
                    }
#endif

                    for (int n = 0; n<n_cons(); ++n) {
                        U[n] = data4(i,j,k,n);
                    }

                    const RealArray speed = get_speed_from_cons(U);


                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        max_speed = std::max(max_speed, speed[d]);
                    }
                }
            }
        }
    }

    Real dt = dx[0]/max_speed;
    for (int i=1; i<AMREX_SPACEDIM; ++i) {
        dt = std::min(dt, dx[i]/max_speed);
    }

    if (flux_solver->need_max_speed()) {
        ParallelDescriptor::ReduceRealMin(dt);
        flux_solver->update_max_speed(mfp->cfl*dx[0]/dt);
    }

    return dt;
}

void EulerianState::calc_primitives(const Box& box,
                                 FArrayBox& cons,
                                 FArrayBox& prim,
                                 const Real* dx,
                                 const Real t,
                                 const Real* prob_lo
                                 #ifdef AMREX_USE_EB
                                 ,const FArrayBox& vfrac
                                 #endif
                                 ) const
{
    BL_PROFILE("MHDState::calc_primitives");

    Vector<Real> U(n_cons());
    Vector<Real> Q(n_prim());

    prim.resize(box, n_prim());

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);
    Array4<Real> const& s4 = cons.array();
    Array4<Real> const& p4 = prim.array();

#ifdef AMREX_USE_EB
    Array4<const Real> const& vfrac4 = vfrac.array();

    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);
#endif

    Real x, y, z;

    for     (int k = lo.z; k <= hi.z; ++k) {
        z = prob_lo[2] + (k + 0.5)*dx[2];
        for   (int j = lo.y; j <= hi.y; ++j) {
            y = prob_lo[1] + (j + 0.5)*dx[1];
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                x = prob_lo[0] + (i + 0.5)*dx[0];

#ifdef AMREX_USE_EB
                if (vfrac4(i,j,k) == 0.0) {

                    // iterate over all neighbouring cells checking if it has valid data
                    // from these calculate the volume weighted average to populate the
                    // covered cell
                    Real vtot = 0.0;

                    std::fill(U.begin(), U.end(), 0.0);
                    for (const auto& index : grab) {

                        const int ii = i+index[0];
                        const int jj = j+index[1];
                        const int kk = k+index[2];

                        // make sure our stencil is within bounds
                        if ((lo.x > ii) || (ii > hi.x) ||
                                (lo.y > jj) || (jj > hi.y) ||
                                (lo.z > kk) || (kk > hi.z)) continue;

                        const Real vf = vfrac4(ii,jj,kk);
                        if (vf > 0.0) {
                            for (int n=0; n<n_cons(); ++n) {
                                U[n] += vf*s4(ii,jj,kk,n);
                            }

                            vtot += vf;
                        }
                    }

                    // if we were close enough to a boundary to have valid data we
                    // average out the volume fraction weighted contributions, otherwise,
                    // fill in the primitives with zeros
                    if (vtot > 0.0) {
                        for (int n=0; n<n_cons(); ++n) {
                            U[n] /= vtot;
                        }
                    } else {
                        for (int n=0; n<n_prim(); ++n) {
                            p4(i,j,k,n) = 0.0;
                        }
                        continue;
                    }
                } else {
#endif
                    // grab the conserved variables
                    //Print() << "\ncons indx 0...:\t" ; //TODO delete me 
                    for (int n=0; n<n_cons(); ++n) {
                        U[n] = s4(i,j,k,n);
                        //Print() << U[n] << " " ; //TODO delete me 
                    }
#ifdef AMREX_USE_EB
                }
#endif


                // convert to primitive
                cons2prim(U, Q);

                // copy into primitive
                for (int n=0; n<n_prim(); ++n) {
                    p4(i,j,k,n) = Q[n];
                }
            }
        }
    }

    return;
}

void EulerianState::update_boundary_cells(const Box& box,
                                       const Geometry &geom,
                                       FArrayBox &prim,
                                       #ifdef AMREX_USE_EB
                                       const FArrayBox& vfrac,
                                       #endif
                                       const Real time) const
{
    BL_PROFILE("MHDState::update_boundary_cells");

    // NOTE: This function could be skipped if a check is done at init.

    Vector<BoundaryInfo> limits = get_bc_limits(box, geom);

    if (limits.empty())
        return;

    const Box& domain = geom.Domain();
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    Real x, y, z;
    std::map<std::string, Real> Q{{"t", time}};

    Array<int,3> grab;

    Array4<Real> const& p4 = prim.array();

#ifdef AMREX_USE_EB
    Array4<const Real> const& vfrac4 = vfrac.array();
#endif

    const Vector<std::string> prim_names = get_prim_names();

    const BCRec &bc = boundary_conditions.where_is_inflow;
    for (const auto &L : limits) {
        if (((L.lo_hi == 0) && (bc.lo(L.dir) == BCType::ext_dir))||
                ((L.lo_hi == 1) && (bc.hi(L.dir) == BCType::ext_dir))) {

            if (L.lo_hi == 0) {
                grab[L.dir] = domain.smallEnd(L.dir);
            } else {
                grab[L.dir] = domain.bigEnd(L.dir);
            }

            for (int k=L.kmin; k<=L.kmax; ++k) {
                z = prob_lo[2] + (k + 0.5)*dx[2];
                Q["z"] = z;
                if (L.dir != 2) {
                    grab[2] = k;
                }
                for (int j=L.jmin; j<=L.jmax; ++j) {
                    y = prob_lo[1] + (j + 0.5)*dx[1];
                    Q["y"] = y;
                    if (L.dir != 1) {
                        grab[1] = j;
                    }
                    for (int i=L.imin; i<=L.imax; ++i) {
                        x = prob_lo[0] + (i + 0.5)*dx[0];
                        Q["x"] = x;
                        if (L.dir != 0) {
                            grab[0] = i;
                        }
#ifdef AMREX_USE_EB
                        if (vfrac4(grab[0],grab[1],grab[2]) == 0.0) {
                            continue;
                        }
#endif

                        // get data from closest internal cell
                        for (int n=0; n<prim_names.size(); ++n) {
                            Q[prim_names[n]] = p4(grab[0],grab[1],grab[2],n);
                        }

                        // update the primitives from our UDFs, but only those that are valid functions
                        for (int n=0; n<prim_names.size(); ++n) {
                            const Optional3D1VFunction &f = boundary_conditions.get(L.lo_hi, L.dir, prim_names[n]);
                            if (f.is_valid()) {
                                p4(i,j,k,n) = f(Q);
                            }
                        }
                    }
                }
            }
        }
    }
}

void EulerianState::calc_reconstruction(const Box& box,
                                     FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                     #ifdef AMREX_USE_EB
                                     ,const EBCellFlagFab &flag
                                     ,const FArrayBox &vfrac
                                     #endif
                                     ) const
{
    BL_PROFILE("MHDState::calc_reconstruction");

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& src4 = prim.array();

#ifdef AMREX_USE_EB

    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);


    Array4<const EBCellFlag> const& f4 = flag.array();
    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    Vector<Real> stencil(reconstructor->stencil_length);
    int offset = reconstructor->stencil_length/2;
    Array<int,3> stencil_index;
    Real lo_face, hi_face;


    // cycle over dimensions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        // make sure our arrays for putting lo and hi reconstructed values into
        // are the corect size
        rlo[d].resize(box, n_prim());
        rhi[d].resize(box, n_prim());

#ifdef AMREX_USE_EB
        if (check_eb) {
            rlo[d].copy(prim,box);
            rhi[d].copy(prim,box);
        }
#endif

        Array4<Real> const& lo4 = rlo[d].array();
        Array4<Real> const& hi4 = rhi[d].array();

        // cycle over all components
        for (int n = 0; n<n_cons(); ++n) {

            for     (int k = lo.z; k <= hi.z; ++k) {
                for   (int j = lo.y; j <= hi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                            for (int i = lo.x; i <= hi.x; ++i) {



#ifdef AMREX_USE_EB
                        if (check_eb) {
                            // covered cell doesn't need calculating
                            if (f4(i,j,k).isCovered() || check_covered_stencil(f4, i, j, k, d, reconstructor->stencil_length)) {
                                continue;
                            }
                        }
#endif
                        stencil_index.fill(0);
                        for (int s=0; s<reconstructor->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            stencil[s] = src4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], n);
                        }

                        // perform reconstruction
                        reconstructor->get_face_values(stencil, lo_face, hi_face);

                        lo4(i,j,k,n) = lo_face;
                        hi4(i,j,k,n) = hi_face;
                    }
                }
            }
        }
    }


    return;
}

void EulerianState::calc_time_averaged_faces(const Box& box,
                                          const FArrayBox &prim,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                          Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                                          #ifdef AMREX_USE_EB
                                          const EBCellFlagFab& flag,
                                          #endif
                                          const Real* dx,
                                          Real dt) const
{
    BL_PROFILE("MHDState::calc_time_averaged_faces");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& p4 = prim.array();

    Vector<Real> Q(n_prim());

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Real lo_face, centre, hi_face, a6;
    Real dt_2dx;

    const Real ft = 4.0/3.0;

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    continue;
                }
#endif

                for (int n=0; n<n_prim(); ++n) {
                    Q[n] = p4(i,j,k,n);
                }

                const RealArray local_c = get_speed_from_prim(Q);


                for (int d=0; d<AMREX_SPACEDIM; ++d) {

                    Array4<Real> const& lo4 = rlo[d].array();
                    Array4<Real> const& hi4 = rhi[d].array();

                    dt_2dx = 0.5*dt/dx[d];

                    // cycle over all components
                    for (int n=0; n<n_prim(); ++n) {

                        lo_face = lo4(i,j,k,n);
                        hi_face = hi4(i,j,k,n);
                        centre = p4(i,j,k,n);

                        a6 = 6.0*centre - 3*(lo_face + hi_face);

                        lo4(i,j,k,n) += local_c[d]*dt_2dx*(hi_face - lo_face + (1 - local_c[d]*ft*dt_2dx)*a6);
                        hi4(i,j,k,n) -= local_c[d]*dt_2dx*(hi_face - lo_face - (1 - local_c[d]*ft*dt_2dx)*a6);

                    }
                }
            }
        }
    }

    return;
}

/*
* the data passed to this function is indexed by cell
* but resides at the interface
*/

void EulerianState::face_bc(const int dir,
                         Box const& box,
                         const FArrayBox& src,
                         FArrayBox& dest,
                         const Geometry &geom,
                         #ifdef AMREX_USE_EB
                         const EBCellFlagFab& flag,
                         #endif
                         const Real time,
                         const bool do_all) const
{
    BL_PROFILE("MHDState::face_bc");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    const Box& domain = geom.Domain();
    const Dim3 domlo = amrex::lbound(domain);
    const Dim3 domhi = amrex::ubound(domain);

    const Real* prob_lo = geom.ProbLo();
    const Real* dx = geom.CellSize();

    Real x, y, z;
    std::map<std::string, Real> Q{{"t", time}};

    // define the limits of our operations
    Vector<BoundaryInfo> limits;

    BoundaryInfo bi;
    bi.imin = lo.x;
    bi.imax = hi.x;
    bi.jmin = lo.y;
    bi.jmax = hi.y;
    bi.kmin = lo.z;
    bi.kmax = hi.z;



    if (dir == 0) {

        int ilo = domlo.x;
        int ihi = domhi.x + 1; // account for face indexing

        if (lo.x == ilo) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.imin = lo.x;
            b.imax = lo.x;
            b.lo_hi = 0;
        }

        if (hi.x == ihi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.imin = hi.x;
            b.imax = hi.x;
            b.lo_hi = 1;
        }
    }

#if AMREX_SPACEDIM >= 2
    if (dir == 1) {
        int jlo = domlo.y;
        int jhi = domhi.y + 1; // account for face indexing

        if (lo.y == jlo) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.jmin = lo.y;
            b.jmax = lo.y;
            b.lo_hi = 0;
        }

        if (hi.y == jhi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.jmin = hi.y;
            b.jmax = hi.y;
            b.lo_hi = 1;
        }
    }
#endif

#if AMREX_SPACEDIM == 3
    if (dir == 2) {
        int klo = domlo.z;
        int khi = domhi.z + 1; // account for face indexing

        if (lo.z == klo) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.kmin = lo.z;
            b.kmax = lo.z;
            b.lo_hi = 0;
        }

        if (hi.z == khi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.kmin = hi.z;
            b.kmax = hi.z;
            b.lo_hi = 1;
        }
    }
#endif

    Array4<Real> const& d4 = dest.array();
    Array4<const Real> const& s4 = src.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    if (do_all) {
        // first do the usual boundary conditions
        for (int n=0; n<n_prim(); ++n) {
            const BCRec &bc = boundary_conditions.fill_bc[n];

            for (const auto &L : limits) {

                if (
                        ((L.lo_hi == 0) && (bc.lo(dir) == BCType::foextrap))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::foextrap))||

                        ((L.lo_hi == 0) && (bc.lo(dir) == BCType::hoextrap))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::hoextrap))||

                        ((L.lo_hi == 0) && (bc.lo(dir) == BCType::reflect_even))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::reflect_even))
                        ) {
                    for (int k=L.kmin; k<=L.kmax; ++k) {
                        for (int j=L.jmin; j<=L.jmax; ++j) {
                            for (int i=L.imin; i<=L.imax; ++i) {

#ifdef AMREX_USE_EB
                                if (f4(i,j,k).isCovered()) {
                                    continue;
                                }
#endif


                                d4(i,j,k,n) = s4(i,j,k,n);
                            }
                        }
                    }
                }
                if (((L.lo_hi == 0) && (bc.lo(dir) == BCType::reflect_odd))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::reflect_odd))) {
                    for (int k=L.kmin; k<=L.kmax; ++k) {
                        for (int j=L.jmin; j<=L.jmax; ++j) {
                            for (int i=L.imin; i<=L.imax; ++i) {

#ifdef AMREX_USE_EB
                                if (f4(i,j,k).isCovered()) {
                                    continue;
                                }
#endif


                                d4(i,j,k,n) = -s4(i,j,k,n);
                            }
                        }
                    }
                }
            }
        }
    }

    const Vector<std::string> prim_names = get_prim_names();

    // now do any dirichlet boundaries
    Array<Real,3> offset = {0.5, 0.5, 0.5};
    offset[dir] = 0.0;
    const BCRec &bc = boundary_conditions.where_is_inflow;
    for (const auto &L : limits) {
        if (((L.lo_hi == 0) && (bc.lo(dir) == BCType::ext_dir))||
                ((L.lo_hi == 1) && (bc.hi(dir) == BCType::ext_dir))) {

            for (int k=L.kmin; k<=L.kmax; ++k) {
                z = prob_lo[2] + (k+offset[2])*dx[2];
                Q["z"] = z;
                for (int j=L.jmin; j<=L.jmax; ++j) {
                    y = prob_lo[1] + (j+offset[1])*dx[1];
                    Q["y"] = y;
                    for (int i=L.imin; i<=L.imax; ++i) {
                        x = prob_lo[0] + (i+offset[0])*dx[0];
                        Q["x"] = x;

#ifdef AMREX_USE_EB
                        if (f4(i,j,k).isCovered()) {
                            continue;
                        }
#endif

                        // load data from the other side of the face
                        for (int n=0; n<prim_names.size(); ++n) {
                            Q[prim_names[n]] = s4(i,j,k,n);
                        }

                        // update the primitives from our UDFs, but only those that are valid functions
                        for (int n=0; n<prim_names.size(); ++n) {
                            const Optional3D1VFunction &f = boundary_conditions.get(L.lo_hi, dir, prim_names[n]);
                            if (f.is_valid()) {
                                d4(i,j,k,n) = f(Q);
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}


void EulerianState::update_face_prim(const Box& box, const Geometry& geom,
                                     Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                                     #ifdef AMREX_USE_EB
                                     const EBCellFlagFab& flag,
                                     #endif
                                     const Real time, const bool do_all) const
{
    BL_PROFILE("EulerianState::update_face_prim");
    const Box domain = geom.Domain();
    const int*     domainlo    = domain.loVect();
    const int*     domainhi    = domain.hiVect();



    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (geom.isPeriodic(d)) {
            continue;
        }

        if (box.smallEnd(d) <= domainlo[d]) {
            Box lo = box;
            lo.setBig(d,domainlo[d]);
            lo.setSmall(d,domainlo[d]);

            // get the left and right reconstructed primitives for all faces
            FArrayBox &L = r_hi[d];
            FArrayBox &R = r_lo[d];

            // adjust the boxes so that we can use the same index
            L.shift(d,1);

            face_bc(d,
                    lo,
                    R,
                    L,
                    geom,
        #ifdef AMREX_USE_EB
                    flag,
        #endif
                    time,
                    do_all);

            // change it back
            L.shift(d,-1);
        }

        if (box.bigEnd(d) >= domainhi[d]) {
            Box hi = box;
            hi.setBig(d,  domainhi[d]+1);
            hi.setSmall(d,domainhi[d]+1);

            // get the left and right reconstructed primitives for all faces
            FArrayBox &L = r_hi[d];
            FArrayBox &R = r_lo[d];

            // adjust the boxes so that we can use the same index
            L.shift(d,1);

            face_bc(d,
                    hi,
                    L,
                    R,
                    geom,
        #ifdef AMREX_USE_EB
                    flag,
        #endif
                    time,
                    do_all);

            // change it back
            L.shift(d,-1);
        }
    }
}


void EulerianState::calc_fluxes(const Box& box,
                             FArrayBox &cons,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                             Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                             #ifdef AMREX_USE_EB
                             const EBCellFlagFab& flag,
                             #endif
                             const Real *dx,
                             const Real dt) const
{
    BL_PROFILE("EulerianState::calc_fluxes");

    size_t ncons = n_cons();
    size_t nprim = n_prim();

    Array<int, 3> index;
    Vector<Real> L(nprim), R(nprim);
    Vector<Real> F(ncons);

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Vector<int> prim_vector_idx = get_prim_vector_idx();
    Vector<int> cons_vector_idx = get_cons_vector_idx();

    // cycle over dimensions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        FArrayBox& flux = fluxes[d];
        flux.setVal(0.0);
        Array4<Real> const& flux4 = flux.array();

        const Box &fbox = flux.box();
        const Dim3 flo = amrex::lbound(fbox);
        const Dim3 fhi = amrex::ubound(fbox);

        index.fill(0);
        index[d] = 1;

        Array4<const Real> lo4, hi4;

        lo4 = r_lo[d].array();
        hi4 = r_hi[d].array();

        for     (int k = flo.z; k <= fhi.z; ++k) {
            for   (int j = flo.y; j <= fhi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = flo.x; i <= fhi.x; ++i) {

#ifdef AMREX_USE_EB

                    // both cells have to be covered before we skip
                    if (f4(i,j,k).isCovered() && f4(i-index[0],j-index[1],k-index[2]).isCovered())
                        continue;
#endif

                    // get left and right states
                    for (size_t n=0; n<nprim; ++n) {
                        L[n] = hi4(i-index[0],j-index[1],k-index[2],n);
                        R[n] = lo4(i,j,k,n);
                    }


                    // rotate the vectors
                    transform_global2local(L, d, prim_vector_idx);
                    transform_global2local(R, d, prim_vector_idx);


                    Real shk = 0.0;

                    if (shock_detector) {
                        shk = shock_detector->solve(L, R);
                    }

                    flux_solver->solve(L, R, F, &shk);

                    // rotate the flux back to local frame
                    transform_local2global(F, d, cons_vector_idx);

                    // load the flux into the array
                    for (int n=0; n<ncons; ++n) {
                        flux4(i,j,k,n) += F[n];

                        AMREX_ASSERT(std::isfinite(flux4(i,j,k,n)));

                    }

                }
            }
        }
    }

    return;
}

void EulerianState::correct_face_prim(const Box& box,
                                   Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                                   Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                                   const Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                   #ifdef AMREX_USE_EB
                                   const EBCellFlagFab &flag,
                                   #endif
                                   const Real *dx,
                                   const Real dt) const
{
    BL_PROFILE("State::correct_face_prim");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Vector<Real> L_prim(n_prim()), R_prim(n_prim());
    Vector<Real> L_cons(n_cons()), R_cons(n_cons());

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();

    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    // the below loops selects the alternative dimensions to calculate the corrections from
    // e.g. 1D: d1 = 0, d2 = -
    // e.g. 2D: d1 = 0, d2 = 1
    //          d1 = 1, d2 = 0
    // e.g. 3D: d1 = 0, d2 = 1, 2
    //          d1 = 1, d2 = 2, 0
    //          d1 = 2, d2 = 0, 1

    Array<int, 3> idx1, idx2;

    Real hiF, loF;
    int d1, dd, d2;

    for (d1=0; d1<AMREX_SPACEDIM; ++d1) { // face direction
        idx1.fill(0);
        idx1[d1] = 1;

        Array4<Real> const &lo4 = r_lo[d1].array();
        Array4<Real> const &hi4 = r_hi[d1].array();

        for (dd=1; dd<AMREX_SPACEDIM; ++dd) {
            d2 = (d1+dd)%AMREX_SPACEDIM; // flux direction

            idx2.fill(0);
            idx2[d2] = 1;

            Array4<const Real> const &flux4 = fluxes[d2].array();

            for     (int k = lo.z; k <= hi.z + idx1[2]; ++k) {
                for   (int j = lo.y; j <= hi.y + idx1[1]; ++j) {
                    AMREX_PRAGMA_SIMD
                            for (int i = lo.x; i <= hi.x + idx1[0]; ++i) {

#ifdef AMREX_USE_EB
                        if (check_eb) {
                            // only calculate corrections for faces where the stencil of fluxes is full

                            // a) this interface has a valid flux
                            if (f4(i,j,k).isDisconnected(-idx1[0],-idx1[1],-idx1[2]))
                                continue;

                            // b) right cell has valid fluxes hi and lo
                            if (f4(i,j,k).isDisconnected( idx2[0], idx2[1], idx2[2]))
                                continue;
                            if (f4(i,j,k).isDisconnected(-idx2[0],-idx2[1],-idx2[2]))
                                continue;

                            // c) left cell has valid fluxes hi and lo
                            if (f4(i-idx1[0],j-idx1[1],k-idx1[2]).isDisconnected( idx2[0], idx2[1], idx2[2]))
                                continue;
                            if (f4(i-idx1[0],j-idx1[1],k-idx1[2]).isDisconnected(-idx2[0],-idx2[1],-idx2[2]))
                                continue;
                        }
#endif

                        // get left and right states
                        for (int n=0; n<n_prim(); ++n ) {
                            L_prim[n] = hi4(i-idx1[0],j-idx1[1],k-idx1[2],n);
                            R_prim[n] = lo4(i,j,k,n);
                        }

                        // convert to conserved
                        prim2cons(L_prim, L_cons);
                        prim2cons(R_prim, R_cons);

                        /* example for x-direction in 2D
                                   * L & R = reconstructed x- face values
                                   * hiF & loF = fluxes in y- direction on hi and lo faces
                                   *   relative to the reconstructed values
                                   * diagonal lines represent the contributing factors to L & R
                                   *
                                   * L* = L + 0.5*dt*dy*(loF_L - hiF_L)
                                   * R* = R + 0.5*dt*dy*(loF_R - hiF_R)
                                   *
                                   *
                                   * | - - hiF_L - - | - - hiF_R - - |
                                   * |        \      |     /         |
                                   * |          \    |   /           |
                                   * |             L | R             |
                                   * |          /    |   \           |
                                   * |        /      |     \         |
                                   * | - - loF_L - - | - - loF_R - - |
                                   */

                        // left side
                        for (int n=0; n<n_cons(); ++n ) {
                            loF = flux4(i-idx1[0],j-idx1[1],k-idx1[2],n);
                            hiF = flux4(i-idx1[0]+idx2[0],j-idx1[1]+idx2[1],k-idx1[2]+idx2[2],n);

                            L_cons[n] += 0.5*dt/dx[d2]*(loF - hiF);
                        }

                        // right side
                        for (int n=0; n<n_cons(); ++n ) {
                            loF = flux4(i,j,k,n);
                            hiF = flux4(i+idx2[0],j+idx2[1],k+idx2[2],n);

                            R_cons[n] += 0.5*dt/dx[d2]*(loF - hiF);
                        }

                        // convert to primitive
                        cons2prim(L_cons, L_prim);
                        cons2prim(R_cons, R_prim);

                        // update reconstruction
                        for (int n=0; n<n_prim(); ++n ) {
                            hi4(i-idx1[0],j-idx1[1],k-idx1[2],n) = L_prim[n];
                            lo4(i,j,k,n) = R_prim[n];
                        }
                    }
                }
            }
        }
    }
}



#ifdef AMREX_USE_EB
Real EulerianState::interp2d(Real cym,
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


void EulerianState::calc_wall_fluxes(const Box& box,
                                  const Vector<FArrayBox*> &all_prim,
                                  Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                  const EBCellFlagFab& flag,
                                  const CutFab &bc_idx,
                                  const FArrayBox& bcent,
                                  const FArrayBox &bnorm,
                                  const Array<const FArrayBox*, AMREX_SPACEDIM> &afrac,
                                  const Real *dx,
                                  const Real dt) const
{
    BL_PROFILE("HydroState::calc_wall_fluxes");

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Vector<Array4<const Real>> p4(all_prim.size());

    for (size_t i=0; i<all_prim.size(); ++i) {
        if (all_prim[i]) {
            p4[i] = all_prim[i]->array();
        }
    }


    Array4<const Real> const& bcent4 = bcent.array();
    Array4<const Real> const& bnorm4 = bnorm.array();

    Array<Array4<Real>,AMREX_SPACEDIM> flux4;
    Array<Array4<const Real>,AMREX_SPACEDIM> afrac4;

    Array<Vector<Real>,AMREX_SPACEDIM> wall_flux;
    for (int i=0; i<AMREX_SPACEDIM;++i) { wall_flux[i].resize(fluxes[i].nComp());}

    Array4<const EBCellFlag> const& flag4 = flag.array();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        flux4[d] = fluxes[d].array();
        afrac4[d] = afrac[d]->array();

        // zero out the flux accumulator
        fluxes[d].setVal(0.0);
    }

    const Array4<const Real>& bc_idx4 = bc_idx.array();

    Array<Array<Real,3>,3> wall_coord = {{{0,0,0},{0,0,0},{0,0,0}}};
    Array<Real,AMREX_SPACEDIM> wall_centre;

    for (int k = lo.z-AMREX_D_PICK(0,0,2); k <= hi.z+AMREX_D_PICK(0,0,2); ++k) {
        for (int j = lo.y-AMREX_D_PICK(0,2,2); j <= hi.y+AMREX_D_PICK(0,2,2); ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x-2; i <= hi.x+2; ++i) {

                const EBCellFlag &cflag = flag4(i,j,k);

                if (cflag.isSingleValued()) {

                  // the boundary condition
                  const int ebi = (int)nearbyint(bc_idx4(i,j,k));
                  EulerianBoundaryEB& bc = *eb_bcs[ebi];

                    for (int d=0; d<AMREX_SPACEDIM; ++d) {

                        // get the wall normal
                        wall_coord[0][d] = bnorm4(i,j,k,d);

                        // get the centre of the wall
                        wall_centre[d] = bcent4(i,j,k,d);
                    }

                    // get a local coordinate system with x- aligned with the wall normal
                    expand_coord(wall_coord);

                    // calculate the wall flux
                    bc.solve(wall_coord, wall_centre, p4, i, j, k, dx, wall_flux);

                    // load the flux into the fab
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        for (int n=0; n<fluxes[d].nComp(); ++n) {
                            flux4[d](i,j,k,n) += wall_flux[d][n];
                        }
                    }
                }
            }
        }
    }
    return;
}


void EulerianState::calc_eb_divergence(const Box& box,
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
    BL_PROFILE("EulerianState::calc_eb_divergence");
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

bool is_inside(const int i,const int j, const int k, const Dim3 &lo, const Dim3 &hi)
{
    return  i >= lo.x && i <= hi.x &&
            j >= lo.y && j <= hi.y &&
            k >= lo.z && k <= hi.z;
}

// https://stackoverflow.com/a/53283994
struct ArrayHasher {
    int N;
    ArrayHasher(int n){
        N = n;
    }
    int operator()(const Array<int,3> &V) const {
        return V[0] + N*(V[1] + N*V[2]);
    }
};

struct ArrayEqual {
    int operator()(const Array<int,3> &V1, const Array<int,3> &V2) const {
        return (V1[0] == V2[0]) &&  (V1[1] == V2[1]) && (V1[2] == V2[2]);
    }
};

void EulerianState::merge_cells(const Box& box,
                                FArrayBox &cons,
                                FArrayBox& du,
                                const EBCellFlagFab& flag,
                                const FArrayBox& vfrac,
                                const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                                int as_fine,
                                FArrayBox &dm_as_fine,
                                const IArrayBox& levmsk) const
{
    BL_PROFILE("EulerianState::merge_cells");
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

    // the sets of cells that are to be merged
    Vector<Vector<Array<int,3>>> merge;

    // a container for all of the unique cells that are part of the problem and
    // the super-cell that they are a part of
    int Ni = halo.size().max(); // the largest logical size of the data
    int Nc = halo.numPts()*0.1; // guess at how many cells will need to be merged
    const ArrayHasher hash(Ni);
    const ArrayEqual equal;
    std::unordered_map<Array<int,3>,int, ArrayHasher,ArrayEqual> cells(Nc, hash, equal);


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

bool EulerianState::check_covered_stencil(Array4<const EBCellFlag> const& flag, int i, int j, int k, int d, int stencil_length)
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


#endif

void EulerianState::calc_divergence(const Box& box,
                                    Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                    FArrayBox& du,
                                    const Real *dx,
                                    const Real dt) const
{
    BL_PROFILE("State::calc_divergence");
    // make sure du is empty
    du.setVal(0.0);

    int N = du.nComp();

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<int,3> index = {0,0,0};
    Array4<Real> const& du4 = du.array();
    Real dxinv;

    Array<Array4<const Real>,AMREX_SPACEDIM> flux4;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        flux4[d] = fluxes[d].array();
    }

    for (int n=0; n<N; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

                    // cycle over dimensions
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        index[d] = 1;
                        dxinv = dt/dx[d];

                        const Real& flux_lo  = flux4[d](i,j,k,n);

                        const Real& flux_hi  = flux4[d](i+index[0], j+index[1], k+index[2], n);

                        du4(i,j,k,n) += dxinv*(flux_lo - flux_hi);


                        index[d] = 0;
                    }
                }
            }
        }
    }

    return;
}

void EulerianState::calc_slope(const Box& box,
                               const FArrayBox& src,
                               FArrayBox &slope,
                               #ifdef AMREX_USE_EB
                               const EBCellFlagFab &flag,
                               #endif
                               const Real *dx,
                               const int src_idx,
                               const int slope_idx,
                               const int dim,
                               Reconstruction &reco)
{
    BL_PROFILE("EulerianState::calc_slope");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& src4 = src.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
    bool check_eb = flag.getType() != FabType::regular;

#endif

    Vector<Real> stencil(reco.stencil_length);
    int offset = reco.stencil_length/2;
    Array<int,3> stencil_index;

    Array4<Real> const& s4 = slope.array();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {


                stencil_index.fill(0);

#ifdef AMREX_USE_EB
                if (check_eb) {
                    // covered cell doesn't need calculating
                    if (f4(i,j,k).isCovered() || check_covered_stencil(f4, i, j, k, dim, reco.stencil_length)) {
                        s4(i,j,k) = 0.0;
                        continue;
                    }
                }
#endif

                for (int s=0; s<reco.stencil_length; ++s) {
                    stencil_index[dim] = s - offset;
                    stencil[s] = src4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], src_idx);
                }

                // perform reconstruction
                s4(i,j,k,slope_idx) = reco.get_slope(stencil)/dx[dim]; // account for local cell size
            }
        }
    }


    return;
}


void EulerianState::write_info(nlohmann::json &js) const
{

    State::write_info(js);

    // write out stuff that is common to all states

    js["num_grow"] = num_grow;

    js["cons_names"] = get_cons_names();

}
