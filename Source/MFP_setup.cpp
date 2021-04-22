
#include <MFP.H>
#include <MFP_fillbc.H>
#include <MFP_utility.H>

#include <functional>


void MFP::variableSetUp() {
    BL_PROFILE("MFP::variableSetUp");
    // read in all of the user defined parameters
    read_params();

    bool state_data_extrap = false;
    bool store_in_checkpoint = true;
    Vector<std::string> comp_names;

    // select the interpolation routine
#ifdef AMREX_USE_EB
    Interpolater* interp = &eb_cell_cons_interp;
#else
    Interpolater* interp = &cell_cons_interp;
#endif

    for (int idx=0; idx<gd.num_solve_state; ++idx) {

        State &istate = *gd.states[idx];

        int nc = istate.n_cons();
        int np = istate.n_prim();
        comp_names.resize(nc);

        // get the list of bc setting functions that applies to this type of state
        const Vector<set_bc> &setbc = gd.states[idx]->get_bc_set();
        const Vector<std::string> &cons_names = gd.states[idx]->get_cons_names();

        istate.boundary_conditions.fill_bc.resize(np);

        for (int icomp=0; icomp < np; ++icomp) {
            set_bc s = setbc[icomp]; // the function that sets the bc

            // make sure our periodicity isn't being overwritten
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                if (gd.periodic[d]) {
                    istate.boundary_conditions.phys_fill_bc[icomp].setLo(d, PhysBCType::interior);
                    istate.boundary_conditions.phys_fill_bc[icomp].setHi(d, PhysBCType::interior);
                }
            }

            // grab the per component BCRec and apply the set bc function to it
            (*s)(istate.boundary_conditions.fill_bc[icomp], istate.boundary_conditions.phys_fill_bc[icomp]);
        }

        for (int icomp=0; icomp < nc; ++icomp) {
            comp_names[icomp] = cons_names[icomp] + "-" + gd.state_names[idx];
        }

        int ng = istate.num_grow;

        desc_lst.addDescriptor(idx, IndexType::TheCellType(),
                               StateDescriptor::Point, ng, nc,
                               interp, state_data_extrap,
                               store_in_checkpoint);

        desc_lst.setComponent(
                    idx, 0, comp_names, istate.boundary_conditions.fill_bc,
                    FillBC());


        if (gd.verbose >= 1) {
            Print() << istate.str();
        }

    }

    //===================
    // cost state

    gd.Cost_Idx = desc_lst.size();

    // just grab any old BCRec, its not applied anyway
    BCRec bc;
    set_scalar_bc(bc, gd.states[0]->boundary_conditions.fill_bc[0]);

    desc_lst.addDescriptor(gd.Cost_Idx, IndexType::TheCellType(),
                           StateDescriptor::Point, 0, 1, &pc_interp);
    desc_lst.setComponent(gd.Cost_Idx, 0, "Cost", bc,
                          NullFillBC());

    //===================
    // shock tracking state
    gd.Shock_Idx = -1;
    if (gd.plot_shock_detector && gd.num_shock_detector) {
        gd.Shock_Idx = desc_lst.size();

        desc_lst.addDescriptor(gd.Shock_Idx, IndexType::TheCellType(),
                               StateDescriptor::Point, 0, gd.num_shock_detector, &pc_interp);
        desc_lst.setComponent(gd.Shock_Idx, 0, "Shock", bc,
                              NullFillBC());
    }

    //===================
    gd.num_states = desc_lst.size();

    StateDescriptor::setBndryFuncThreadSafety(true);

    return;
}

void MFP::init_data(const Box& box,
                    EB_OPTIONAL(const EBCellFlagFab& flag,)
                    FArrayBox& src,
                    const Real* dx,
                    const Real* prob_lo,
                    const int idx)
{
    BL_PROFILE("MFP::init_data");

    Vector<Real> U;

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);
    Array4<Real> const& h4 = src.array();

    State &istate = *gd.states[idx];

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& flag4 = flag.array();
#endif

    Real x, y, z;
    int n_prim = istate.n_prim();
    int n_cons = istate.n_cons();
    const Vector<std::string> &prim_names = gd.states[idx]->get_prim_names();

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
                    for (int n=0; n<n_cons; ++n) {
                        h4(i,j,k,n) = 0.0;
                    }
                    continue;
                }
#endif

                U.resize(n_prim);

                // grab the primitive variables as defined by the user functions
                for (int icomp=0; icomp<n_prim; ++icomp) {
                    const std::string& prim_name = prim_names[icomp];
                    const auto& f = istate.functions[prim_name];

                    U[icomp] = f(x, y, z);

                }

                // convert primitive to conserved
                istate.prim2cons(U);

                // copy into array
                for (int n=0; n<n_cons; ++n) {
                    h4(i,j,k,n) = U[n];
                }
            }
        }
    }

    return;
}

void MFP::variableCleanUp() {
    BL_PROFILE("MFP::variableCleanUp");
#ifdef AMREX_PARTICLES
    if (gd.do_tracer_particles) {
        for (AmrTracerParticleContainer* TracerPC : particles) {
            if (TracerPC) {
                delete TracerPC;
                TracerPC = 0;
            }
        }
    }
#endif

    desc_lst.clear();
}
