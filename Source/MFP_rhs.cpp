#include "MFP.H"

#include <AMReX_Array.H>
#include <AMReX_FArrayBox.H>
#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#endif

// NOTE: face source term function expects conserved quantities but the reconstructed
// face values are primitives!!
void MFP::calc_face_source (const Box& box,
                            Vector<FArrayBox*>& src_dat,
                            Vector<Array<FArrayBox, AMREX_SPACEDIM>> &R_lo,
                            Vector<Array<FArrayBox, AMREX_SPACEDIM> > &R_hi,
                            EB_OPTIONAL(const Vector<const EBCellFlagFab*> &flag,)
                            const Real* dx,
                            Real time,
                            Real dt)
{


    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    const Real* prob_lo = geom.ProbLo();

    Real x, y, z;
    Vector<Real> y0;
    Array<Vector<Real>, AMREX_SPACEDIM> ydot_lo, ydot_hi;

#ifdef AMREX_USE_EB
    Vector<Array4<const EBCellFlag>> flag4(flag.size());
    for (int i=0; i<flag.size(); ++i) {
        flag4[i] = flag[i]->array();
    }
#endif

    int n;
    int cnt;
    int global_idx;

    // calculate the source terms
    for (const auto &ode_ptr : gd.ode_source_terms) {
        const ODESystem& ode = *ode_ptr;

        // check if we actually have anything to do here
        if (ode.has_face_src.empty()) {
            continue;
        }

        // get all our buffers to be the right size
        y0.resize(ode.n_components);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ydot_lo[d].resize(ode.n_components);
            ydot_hi[d].resize(ode.n_components);
        }

        for     (int k = lo.z; k <= hi.z; ++k) {
            z = prob_lo[2] + (k + 0.5)*dx[2];
            for   (int j = lo.y; j <= hi.y; ++j) {
                y = prob_lo[1] + (j + 0.5)*dx[1];
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {
                    x = prob_lo[0] + (i + 0.5)*dx[0];

#ifdef AMREX_USE_EB
                    // check all of the contributing states to make sure that they hold valid data
                    bool skip_cell = false;
                    for (int local_idx = 0; local_idx < ode.n_states; ++local_idx) {
                        global_idx = ode.local2global_index[local_idx]; // what the source term is operating on
                        const EBCellFlag& cflag = flag4[global_idx](i,j,k);

                        if (cflag.isCovered()) {
                            skip_cell = true;
                            break;
                        }
                    }

                    if (skip_cell)
                        continue;
#endif

                    //
                    // copy data from the source data into the aggregated y vector
                    cnt = 0;
                    for (int local_idx = 0; local_idx < ode.n_states; ++local_idx) {
                        global_idx = ode.local2global_index[local_idx]; // what the source term is operating on
                        n = src_dat[global_idx]->nComp(); // how many elements in the source term
                        Array4<Real> const& s4 = src_dat[global_idx]->array();
                        for (int h=0; h<n; ++h) {
                            y0[cnt] = s4(i,j,k,h);
                            cnt++;
                        }
                    }

                    // set all values in buffers to zero
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        std::fill(ydot_lo[d].begin(), ydot_lo[d].end(), 0);
                        std::fill(ydot_hi[d].begin(), ydot_hi[d].end(), 0);
                    }

                    // calculate all face sources
                    for (const int &src_idx : ode.has_face_src) {
                        ode.sources[src_idx]->face_src(x, y, z, time, y0, ydot_lo, ydot_hi);
                    }

                    // unpack ydot
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {

                        for (const int &src_idx : ode.has_face_src) {
                            const auto &src = *ode.sources[src_idx];

                            for (const auto &idx : src.offsets) {

                                Array4<Real> const& lo4 = R_lo[idx.global][d].array();
                                Array4<Real> const& hi4 = R_hi[idx.global][d].array();

                                n = R_lo[idx.global][d].nComp(); // how many elements in the source term

                                // half*dt as we are at a half timestep
                                for (int h=0; h<n; ++h) {
                                    lo4(i,j,k,h) += 0.5*dt*ydot_lo[d][idx.solver+h];
                                    hi4(i,j,k,h) += 0.5*dt*ydot_hi[d][idx.solver+h];
                                }

                            }

                        }
                    }
                }
            }
        }
    }
}




void MFP::calc_cell_source (const Box& box,
                            Vector<FArrayBox*>& src_dat,
                            Vector<FArrayBox*>& dst_dat,
                            EB_OPTIONAL(Vector<const EBCellFlagFab*> &flags,)
                            EB_OPTIONAL(Vector<const FArrayBox*> &volfrac,)
                            Real time, Real dt)
{


    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    Vector<Real> y0, yout;

    int n;
    int cnt;

    // have a cache of all slopes that have been calculated in order to save
    // re-calculation if it is required
    std::map<std::pair<int,int>, Array<FArrayBox, AMREX_SPACEDIM>> slopes;

    Real x, y, z;

#ifdef AMREX_USE_EB
    Vector<Array4<const EBCellFlag>> flag4(flags.size());
    for (int i=0; i<flags.size(); ++i) {
        flag4[i] = flags[i]->array();
    }
#endif

    // calculate the source terms
    for (auto &ode_ptr : gd.ode_source_terms) {
        ODESystem& ode = *ode_ptr;

        y0.resize(ode.n_components);
        yout.resize(ode.n_components);

        // calculate any slope information required by the source term
        std::map<int,Vector<FArrayBox>> slopes;

        for (int solver_idx = 0; solver_idx < ode.n_sources; ++solver_idx) {
            const auto &src = ode.sources[solver_idx];
            const int n_slopes = src->num_slopes();
            if (n_slopes <= 0) continue;

            Vector<FArrayBox>& slope = slopes[solver_idx];

            src->calc_slopes(box,
                             src_dat,
                             slope,
                             EB_OPTIONAL(flags,)
                             dx);
        }

        // link the solver to the data
        ode.solver->init_data(&y0, &yout);


        // go through cell-by-cell and calculate the contribution to the source vector
        for     (int k = lo.z; k <= hi.z; ++k) {
            z = prob_lo[2] + (k + 0.5)*dx[2];
            for   (int j = lo.y; j <= hi.y; ++j) {
                y = prob_lo[1] + (j + 0.5)*dx[1];
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {
                    x = prob_lo[0] + (i + 0.5)*dx[0];

#ifdef AMREX_USE_EB
                    // inform the source if it is holding invalid data
                    for (int local_idx = 0; local_idx < ode.n_states; ++local_idx) {
                        const int global_idx = ode.local2global_index[local_idx]; // what the source term is operating on
                        const EBCellFlag& cflag = flag4[global_idx](i,j,k);

                        ode.set_solver_state_valid(local_idx, !cflag.isCovered());
                    }
#endif

                    //
                    // copy data from the source data into the aggregated y vector
                    cnt = 0;
                    for (int local_idx = 0; local_idx < ode.n_states; ++local_idx) {
                        const int global_idx = ode.local2global_index[local_idx]; // what the source term is operating on
                        n = src_dat[global_idx]->nComp(); // how many elements in the source term
                        Array4<Real> const& s4 = src_dat[global_idx]->array();
                        for (int h=0; h<n; ++h) {
                            y0[cnt] = s4(i,j,k,h);
                            cnt++;
                        }
                    }

                    //
                    // retrieve slopes
                    for (int solver_idx = 0; solver_idx < ode.n_sources; ++solver_idx) {
                        const auto &src = ode.sources[solver_idx];
                        const int n_slopes = src->num_slopes();
                        if (n_slopes <= 0) continue;

                        Vector<FArrayBox>& slope = slopes[solver_idx];

                        // this function de-serializes the slope FABs to ensure that we can index it correctly
                        src->retrieve_slopes(slope, i, j, k);
                    }

                    //
                    // do the computation
                    // loads the final value into yout (yout = y0 + dt*src)

                    int err = ode.solver->solve(x, y, z, time, time+dt, 0);
                    if (err != 0) {
                        amrex::Abort("Failed integration: " + num2str(err));
                    }

                    // unpack y1
                    cnt = 0;
                    for (int local_idx = 0; local_idx < ode.n_states; ++local_idx) {

                        const int global_idx = ode.local2global_index[local_idx]; // what the source term is operating on
                        n = src_dat[global_idx]->nComp(); // how many elements in the source term

                        if (ode.get_solver_state_valid(local_idx)) {
                            Array4<Real> const& d4 = dst_dat[global_idx]->array();
                            Array4<Real> const& s4 = src_dat[global_idx]->array();
                            for (int h=0; h<n; ++h) {
                                d4(i,j,k,h) += yout[cnt] - s4(i,j,k,h); // dt*src = yout - y0
                                //                            d4(i,j,k,h) = yout[cnt] - s4(i,j,k,h); // output the source for visual debugging
                                cnt++;
                            }
                        } else {
                            cnt += n;
                        }
                    }

#ifdef AMREX_USE_EB
                    // reset to all valid
                    for (int solver_idx = 0; solver_idx < ode.n_states; ++solver_idx) {
                        ode.set_solver_state_valid(solver_idx, true);
                    }
#endif
                }
            }
        }

        ode.solver->clear();

    }
}

void MFP::apply_cell_sources(const Real time,
                             const Real dt) {
    BL_PROFILE("MFP::compute_source_dSdt()");

    if (gd.ode_source_terms.empty()) {
        return;
    }

    const Real* dx = geom.CellSize();

    MultiFab& cost = get_new_data(gd.Cost_Idx);

    // get the hydro and fields filled data set
    Vector<MultiFab> filled_state(gd.num_solve_state);

    MultiFab* for_iteration;
    for (int idx = 0; idx < gd.num_solve_state; ++idx) {
        State &istate = gd.get_state(idx);
        int nc = istate.n_cons();
        int ng = istate.reconstruction->num_grow;
        filled_state[idx].define(grids, dmap, nc, ng, MFInfo(), Factory());
#ifdef AMREX_USE_EB
        EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate.eb2_index));
#endif
        FillPatch(*this, filled_state[idx], ng, time, idx, 0, nc);
        for_iteration = &filled_state[idx];
    }

    Vector<FArrayBox*> filled_ptr(gd.num_solve_state,nullptr);
    Vector<FArrayBox*> updated_ptr(gd.num_solve_state,nullptr);

#ifdef AMREX_USE_EB
    Vector<const EBCellFlagFab*> fab_flags(gd.num_solve_state);
    Vector<const FArrayBox*> fab_vfrac(gd.num_solve_state);
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(*for_iteration);
             mfi.isValid(); ++mfi) {

#ifdef AMREX_USE_EB
            for (int idx=0; idx<gd.num_solve_state; ++idx) {
                const EBCellFlagFab& flag = getEBData(idx).flags[mfi];
                fab_flags[idx] = &flag;

                const FArrayBox& vfrac = getEBData(idx).volfrac[mfi];
                fab_vfrac[idx] = &vfrac;
            }
#endif

            Real wt = ParallelDescriptor::second();

            // get extents
            const Box& bx = mfi.tilebox();

            // get pointers
            for (int si = 0; si < gd.num_solve_state; ++si) {
                filled_ptr[si] = &filled_state[si][mfi];
                updated_ptr[si] = &get_new_data(si)[mfi];
            }

            // pass the data off to the function that calculate the source terms
            calc_cell_source(bx,
                             filled_ptr,
                             updated_ptr,
                             EB_OPTIONAL(fab_flags,)
                             EB_OPTIONAL(fab_vfrac,)
                             time,
                             dt);

            // update the cost function
            wt = (ParallelDescriptor::second() - wt) / bx.d_numPts();
            cost[mfi].plus(wt, bx);

        }
    }

//    for (int si = 0; si < gd.num_solve_state; ++si) {
//        State& istate = gd.get_state(si);
//        MultiFab::Subtract(filled_state[si],get_new_data(si),0,0,istate.n_cons(),istate.reconstruction->num_grow);
//        plot_FAB_2d(filled_state[si], 0, 0, istate.name, false, si == gd.num_solve_state-1);
//    }

    return;
}

