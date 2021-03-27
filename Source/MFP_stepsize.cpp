#include "MFP.H"
#include "MFP_utility.H"

#include <AMReX.H>

//#include "MFP_state.H"


using namespace amrex;


/*
 * This function calculates the estimated time step from all components of the solver
 * It calculates on a per state basis and then grabs time step info from source terms
 * that need all states.
 */
void MFP::estimate_dt (const Box& box,
                       Vector<const FArrayBox*>& src,
                       EB_OPTIONAL(Vector<const FArrayBox*>& vfrac,)
                       const Real* dx,
                       Real& dt)
{

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
    Vector<Array4<const Real>> vfrac4(src.size());
#endif


    // time step due to inviscid fluxes
    for (int idx=0; idx<src.size(); ++idx) {

        if (src[idx] == nullptr) {
            continue;
        }

#ifdef AMREX_USE_EB
        vfrac4[idx] = vfrac[idx]->array();
        Array4<const Real> const& vf4 = vfrac4[idx];
#endif

        State &istate = *gd.states[idx];

        int n_cons = istate.n_cons();

        Vector<Real> U(n_cons);

        Array4<Real const> s4 = src[idx]->array();

        Real state_dt = std::numeric_limits<Real>::max();

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {


#ifdef AMREX_USE_EB
                    const Real vf = vf4(i,j,k);
                    if (vf == 0.0) continue;
#endif

                    // collect the conserved properties
                    for (int n=0; n<n_cons; ++n) {
                        U[n] = s4(i,j,k,n);
                    }

                    // time step restriction due to inviscid wave speeds
                    if (istate.flux_solver) {
                        RealArray s = istate.get_speed_from_cons(U);
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            state_dt = std::min(state_dt, dx[d]/s[d]);
                        }
                    }
                }
            }
        }

        istate.dt = state_dt;

        dt = std::min(dt, state_dt);

    }

    // time step due to viscous fluxes
    Vector<Vector<Real>> UU;


    for (int idx=0; idx<src.size(); ++idx) {

        if (src[idx] == nullptr) {
            continue;
        }

        State &istate = *gd.states[idx];

        if (not istate.viscous)
            continue;

        Viscous &V = *istate.viscous;

        if (V.cfl <= 0)
            continue;

        // get the indexes for states required by this viscous object
        Vector<int> src_states = V.get_linked_states();

        // get the data to support the viscous objects computation
        UU.resize(src_states.size());

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

                    bool skip = false;
                    for (int isrc=0; isrc < src_states.size(); ++isrc) {
                        const int src_idx = src_states[isrc];

#ifdef AMREX_USE_EB
                        Array4<const Real> const& vf4 = vfrac4[src_idx];
                        const Real vf = vf4(i,j,k);
                        if (vf == 0.0) {
                            skip = true;
                            continue;
                        }
#endif


                        Vector<Real>& U = UU[isrc];

                        State &sstate = *gd.states[src_idx];

                        int n_cons = istate.n_cons();

                        U.resize(n_cons);

                        Array4<Real const> s4 = src[src_idx]->array();

                        // collect the conserved properties
                        for (int n=0; n<n_cons; ++n) {
                            U[n] = s4(i,j,k,n);
                        }
                    }

                    if (skip)
                        continue;

                    // now do the viscous time step calculation
                    Real v = V.get_max_speed(UU);
                    for (int d=0; d<AMREX_SPACEDIM; ++d)
                        istate.dt = std::min(istate.dt, dx[d]*dx[d]/v);



                }
            }
        }

        dt = std::min(dt, istate.dt);
    }

    // do collective states from source terms

    Vector<Real> y(gd.full_state_size); // all of the states
    bool all_there;

    for (const auto &ode_ptr : gd.ode_source_terms) {
        const ODESystem& ode = *ode_ptr;
        // check if the solver applied to this collection of ODEs has a time based stability requirement
        if (!ode.solver->has_freq())
            continue;

        for (const std::unique_ptr<SourceTerm> &source: ode.sources) {
            // check if this source term influences the time based stability
            if (!source->has_freq())
                continue;

            // THIS APPROACH POTENTIALLY DUPLICATES COPY OPERATIONS

            Vector<OffsetIndex> &index = source->offsets;

            for     (int k = lo.z; k <= hi.z; ++k) {
                for   (int j = lo.y; j <= hi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                            for (int i = lo.x; i <= hi.x; ++i) {

                        all_there = true;
                        for (const auto &idx : index) {

                            // make sure that we don't have an invalid/covered region
                            if (!src[idx.global]) {
                                all_there = false;
                                break;
                            }

#ifdef AMREX_USE_EB

                            Array4<const Real> const& vf4 = vfrac4[idx.global];
                            const Real vf = vf4(i,j,k);
                            if (vf == 0.0) {
                                all_there = false;
                                continue;
                            }
#endif

                            int n_cons = gd.states[idx.global]->n_cons();

                            Array4<Real const> s4 = src[idx.global]->array();

                            // collect the conserved properties
                            for (int n=0; n<n_cons; ++n) {
                                y[idx.solver + n] = s4(i,j,k,n);
                            }
                        }

                        // now get source term contributions to time step
                        if (all_there) {
                            Real freq = source->get_max_freq(y);
                            Real src_dt = 1/freq;
                            dt = std::min(dt, src_dt);

                            // update the source dt
                            for (const auto &idx : index) {
                                State& istate = *gd.states[idx.global];
                                istate.dt = std::min(src_dt, istate.dt);
                            }
                        }
                    }
                }
            }
        }
    }
}
