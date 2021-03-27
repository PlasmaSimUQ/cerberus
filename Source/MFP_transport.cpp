
#include "MFP.H"

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

using namespace amrex;

#ifdef PYTHON
#include "matplotlibcpp.h"
#include "MFP_diagnostics.H"
namespace plt = matplotlibcpp;
#endif


void MFP::apply_cell_transport(Real time, Real dt) {


#ifdef AMREX_USE_EB
    constexpr int num_grow_eb = 2;
#else
    constexpr int num_grow_eb = 0;
#endif

    // get the maximum wave speed from the time step and cell spacing
    RealVect speed_max;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        speed_max[d] = dt/(gd.cfl*geom.CellSize(d));
    }

    // ==========================================================================
    // all of level storage

    Vector<MultiFab*> global_new(gd.num_solve_state, nullptr);
    Vector<MultiFab> local_new(gd.num_solve_state); // local 'new', gathered from whatever is most recent from the state vector
    Vector<MultiFab> update(gd.num_solve_state); // fluxes, sources
    Vector<FabType> active(gd.num_solve_state); // flag for calculation

#ifdef AMREX_USE_EB
    Vector<EBFluxRegister*> fr_as_crse(gd.num_solve_state, nullptr);
    Vector<EBFluxRegister*> fr_as_fine(gd.num_solve_state, nullptr);
#else
    Vector<YAFluxRegister*> fr_as_crse(gd.num_solve_state, nullptr);
    Vector<YAFluxRegister*> fr_as_fine(gd.num_solve_state, nullptr);
#endif

    for (int idx = 0; idx < gd.num_solve_state; ++idx) {
        State &istate = gd.get_state(idx);

        global_new[idx] = &get_new_data(idx);

        int ns = desc_lst[idx].nComp();
        int ng = istate.num_grow + num_grow_eb;

        local_new[idx].define(grids, dmap, ns, ng, MFInfo(),Factory());
        update[idx].define(grids, dmap, ns, num_grow_eb, MFInfo(),Factory());
        global_new[idx] = &get_new_data(idx);

        // get a full array of data at this level
#ifdef AMREX_USE_EB
        EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate.eb2_index));
#endif
        FillPatch(*this, local_new[idx], ng, time, idx, 0, ns);

//        plot_FAB_2d(local_new[idx], 0, 0, "cons density", false, false);
//        plot_FAB_2d(flag, "flag", true);

        if (istate.reflux && level < parent->finestLevel()) {
            MFP& fine_level = getLevel(level + 1);
            fr_as_crse[idx] = &fine_level.flux_reg[idx];
        }

        if (istate.reflux && level > 0) {
            fr_as_fine[idx] = &flux_reg[idx];
        }

        if (fr_as_crse[idx]) {
            fr_as_crse[idx]->reset();
        }
    }

    // ==========================================================================
    // per state storage

    Vector<FArrayBox*> conserved(gd.num_solve_state);
    Vector<FArrayBox> primitives(gd.num_solve_state);
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> R_lo(gd.num_solve_state), R_hi(gd.num_solve_state);
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> fluxes(gd.num_solve_state);

    Vector<FArrayBox> shock(gd.num_solve_state);
    Vector<FArrayBox*> shock_ptrs(gd.num_solve_state, nullptr);

    int nu; // per state size

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& cost = get_new_data(gd.Cost_Idx);

#ifdef AMREX_USE_EB
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> wall_fluxes(gd.num_solve_state);
    Vector<const EBCellFlagFab*> fab_flags(gd.num_solve_state);
    Vector<const FArrayBox*> fab_vfrac(gd.num_solve_state);
#endif

    MultiFab* mf_shock;
    if (gd.Shock_Idx > 0) {
        mf_shock = &get_new_data(gd.Shock_Idx);
    }


    // iterate over all of the FABs within the level performing reconstruction, flux calculation,
    // and updating the cell-centred data according to the resulting divergence

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();

        // region over which to perform reconstruction for face values
        const Box rbox = amrex::grow(box, 1 + num_grow_eb);

        // ==========================================================================
        // 1. iterate over all states to set-up the data required for flux calculation
        for (int idx=0; idx<gd.num_solve_state; ++idx) {

            State &istate = gd.get_state(idx);

            // get a pointer to the conserved quantities
            conserved[idx] = &local_new[idx][mfi];

#ifdef AMREX_USE_EB
            // get the EB data required for later calls
            const EBCellFlagFab& flag = getEBData(idx).flags[mfi]; fab_flags[idx] = &flag;
            const FArrayBox& vfrac = getEBData(idx).volfrac[mfi]; fab_vfrac[idx] = &vfrac;

//            plot_FAB_2d(flag, "flag", false);
//            plot_FAB_2d(vfrac, 0, "vfrac", false, true);

            // check if this box is active
            active[idx] = flag.getType(rbox);
            if (active[idx] == FabType::covered)
                continue;
#else
            active[idx] = FabType::regular;
#endif
            // final check to see if this state actually has any transport
            // this needs to be here (rather than higher up) as we need
            // certain variables to be defined for other uses (face sources)
            if (!istate.is_transported()) {
                active[idx] = FabType::covered;
                continue;
            }

            // region over which to get cell centered primitives for reconstruction
            const Box pbox = amrex::grow(box, istate.num_grow + num_grow_eb);

            if (is_vacuum(pbox, local_new[idx][mfi], idx)) {
                active[idx] = FabType::covered;
                continue;
            }

            int np = istate.n_prim();

            // ===============================================
            // 1.1 Calculate primitive values within each cell

            FArrayBox& cons = local_new[idx][mfi];
            FArrayBox& prim = primitives[idx];

            prim.resize(pbox, np);

            istate.calc_primitives(pbox,
                                   cons,
                                   prim,
                                   dx,
                                   time,
                                   prob_lo
                                   EB_OPTIONAL(,vfrac)
                                   );

//            plot_FAB_2d(prim, 0, "prim 0", false, true);


            // fill in any cells that need special boundary values
            istate.update_boundary_cells(pbox,
                                    geom,
                                    prim,
                                    EB_OPTIONAL(vfrac,)
                                    time);



            // =======================================
            // 1.2 Calculate reconstructed face values

            // each cell has a hi and lo side in each direction

            // calculate the reconstructed face values
            istate.calc_reconstruction(rbox,
                                       prim,
                                       R_lo[idx],
                                       R_hi[idx]
                                       EB_OPTIONAL(,flag)
                                       EB_OPTIONAL(,vfrac)
                                       );

            // ===================================================================
            // update the face values to time t+1/2 based on the local wave speeds

            // TODO: this step currently assumes a single speed for all components and
            // should be updated to calculate the correct characteristic speeds
            // for each component individually
            if (gd.do_CTU) {
                istate.calc_time_averaged_faces(rbox,
                                                prim,
                                                R_lo[idx],
                                                R_hi[idx],
                                                EB_OPTIONAL(flag,)
                                                dx,
                                                dt);
            }

            // =======================================
            // 1.3 Particle update

            // update particle locations and velocities using the reconstructed
            // face values to interpolate the local velocity for each particle
#ifdef AMREX_PARTICLES
            if (gd.do_tracer_particles && istate.particle_index > -1) {
                AmrTracerParticleContainer* pc = particles[istate.particle_index];
                if (pc) {

                    // grab the starting index for velocity
                    const int vel_idx = istate.get_prim_vector_idx()[0];

                    // grab the tile of particles
                    auto& ptile = pc->ParticlesAt(level, mfi);

                    // update the position and velocity of the tracer particles
                    // using the reconstructed primitive values on cell faces
                    push_particles(ptile,
                                   prim,
                                   R_lo[idx],
                                   R_hi[idx],
                                   vel_idx,
                                   dt
                                   EB_OPTIONAL(,flag)
                                   );

                }
            }
#endif

        }

        // ==========================================================================
        // 2. Apply face value modifications (assumes t = t + dt/2)

        if (gd.do_face_src && gd.do_CTU) {
            calc_face_source(rbox,
                             conserved,
                             R_lo,
                             R_hi,
                             EB_OPTIONAL(fab_flags,)
                             dx,
                             time+dt/2,
                             dt);
        }

        // ==========================================================================
        // 3.1 Setup for flux calculation

        // resize the flux arrays before any get used
        for (int idx=0; idx<gd.num_solve_state; ++idx) {


            if (active[idx] == FabType::covered)
                continue;


            State& istate = gd.get_state(idx);

            // do we have a shock detector?
            if (gd.Shock_Idx > 0 && istate.shock_idx > -1) {
                shock[idx].resize(rbox);
                shock[idx].setVal(0.0);
                shock_ptrs[idx] = &shock[idx];
            }

            nu = local_new[idx].nComp();
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                // get a node centered box in direction d that encloses the original box
                Box fbox = surroundingNodes(box, d);

                if (gd.do_CTU || (istate.viscous && num_grow_eb > 0)) {
                    // enlarge the fluxes for face corrections
                    // grow in all directions but that of the flux
                    IntVect expand(1);
                    expand.setVal(d, 0);
                    fbox = amrex::grow(fbox,expand);
                }

#ifdef AMREX_USE_EB
                fbox = grow(fbox, num_grow_eb);
                wall_fluxes[idx][d].resize(fbox, nu);
#endif
                fluxes[idx][d].resize(fbox, nu);
            }
        }

        // ==========================================================================
        // 3.2 Calculate fluxes

        //---
        for (int idx=0; idx<gd.num_solve_state; ++idx) {

            if (active[idx] == FabType::covered) {
                update[idx][mfi].setVal(0.0);
                continue;
            }

            State &istate = gd.get_state(idx);

            istate.update_face_prim(box,
                                    geom,
                                    R_lo[idx],
                                    R_hi[idx],
                                    EB_OPTIONAL(*fab_flags[idx],)
                                    time);
        }

        //---
        for (int idx=0; idx<gd.num_solve_state; ++idx) {


            if (active[idx] == FabType::covered) continue;


            State &istate = gd.get_state(idx);

            istate.calc_fluxes(box,
                               R_lo,
                               R_hi,
                               fluxes[idx],
                               EB_OPTIONAL(*fab_flags[idx],)
                               dx, dt,
                               shock_ptrs[idx]);

#if AMREX_SPACEDIM > 1
            if (gd.do_CTU) {
                // now that we have done the fluxes, do we need to correct them???
                // correction is done according to the following steps:
                // 1. calculate the fluxes (above)
                // 2. modify the reconstructed face states using the flux information
                // 3. recalculate the fluxes with the updated face values

                istate.correct_face_prim(grow(box,num_grow_eb),
                                         R_lo[idx],
                                         R_hi[idx],
                                         fluxes[idx],
                                         EB_OPTIONAL(*fab_flags[idx],)
                                         dx, dt);

                // following the update of the face values we need to update any boundary conditions

                istate.update_face_prim(box,
                                        geom,
                                        R_lo[idx],
                                        R_hi[idx],
                                        EB_OPTIONAL(*fab_flags[idx],)
                                        time,
                                        true);
            }
#endif
        }

        //---
#if AMREX_SPACEDIM > 1
        if (gd.do_CTU) {
            for (int idx=0; idx<gd.num_solve_state; ++idx) {


                if (active[idx] == FabType::covered) continue;


                State &istate = gd.get_state(idx);

                // recalculate the fluxes
                istate.calc_fluxes(box,
                                   R_lo,
                                   R_hi,
                                   fluxes[idx],
                                   EB_OPTIONAL(*fab_flags[idx],)
                                   dx, dt,
                                   shock_ptrs[idx]);
            }
        }
#endif

        //---
        for (int idx=0; idx<gd.num_solve_state; ++idx) {


            if (active[idx] == FabType::covered) continue;


            State &istate = gd.get_state(idx);

            // now calculate any viscous fluxes

            const Box pbox = amrex::grow(box, istate.num_grow + num_grow_eb);
            FArrayBox& prim = primitives[idx];

            istate.calc_viscous_fluxes(grow(box,num_grow_eb),
                                       fluxes[idx],
                                       pbox, primitives,
                                       EB_OPTIONAL(*fab_flags[idx],)
                                       dx);

            // shock tracking
            if (gd.Shock_Idx > 0 && istate.shock_idx > -1) {
                (*mf_shock)[mfi].copy(shock[idx], 0, istate.shock_idx);
            }

            // given all of the fluxes calculate the update to the cell centred values
            const int as_crse = (fr_as_crse[idx] != nullptr);
            const int as_fine = (fr_as_fine[idx] != nullptr);

            nu = local_new[idx].nComp();

            FArrayBox& dU = update[idx][mfi];


#ifdef AMREX_USE_EB

            Array<const FArrayBox*,AMREX_SPACEDIM> afrac, fcent;

            if (active[idx] != FabType::regular) {

                const FArrayBox &bcent = (*getEBData(idx).bndrycent)[mfi];
                const FArrayBox &bnorm = (*getEBData(idx).bndrynorm)[mfi];

                CutFab& bc_idx = getEBData(idx).bndryidx[mfi];

                afrac = {AMREX_D_DECL(&(*getEBData(idx).areafrac[0])[mfi], &(*getEBData(idx).areafrac[1])[mfi], &(*m_getEBData(idx).areafrac[2])[mfi])};


                // calculate the flux through cut cell faces

                istate.calc_wall_fluxes(box,
                                        primitives,
                                        wall_fluxes[idx],
                                        *fab_flags[idx],
                                        bc_idx,
                                        bcent,
                                        bnorm,
                                        afrac,
                                        dx,
                                        dt);

                // calculate divergence, including cut-cells

                FArrayBox dm_as_fine;
                if (as_fine) {
                    dm_as_fine.resize(amrex::grow(box,2),nu);
                    dm_as_fine.setVal(0.0);
                } else {
                    dm_as_fine.resize(Box::TheUnitBox(),nu);
                }

                FArrayBox fab_drho_as_crse(Box::TheUnitBox(),nu);
                FArrayBox* p_drho_as_crse = (fr_as_crse[idx]) ? fr_as_crse[idx]->getCrseData(mfi) : &fab_drho_as_crse;

                IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());
                const IArrayBox* p_rrflag_as_crse = (fr_as_crse[idx]) ? fr_as_crse[idx]->getCrseFlag(mfi) : &fab_rrflag_as_crse;

                fcent = {AMREX_D_DECL(&(*getEBData(idx).facecent[0])[mfi], &(*getEBData(idx).facecent[1])[mfi], &(*m_getEBData(idx).facecent[2])[mfi])};

                istate.eb_div->calc_eb_divergence(box,
                                                  *conserved[idx],
                                                  fluxes[idx],
                                                  wall_fluxes[idx],
                                                  dU,
                                                  *fab_flags[idx],
                                                  *fab_vfrac[idx],
                                                  afrac,
                                                  fcent,
                                                  as_crse,
                                                  as_fine,
                                                  p_rrflag_as_crse,
                                                  level_mask[mfi],
                                                  p_drho_as_crse,
                                                  dm_as_fine,
                                                  dx,
                                                  dt);

                istate.eb_div->merge_cells(box,
                                           *conserved[idx],
                                           dU,
                                           *fab_flags[idx],
                                           *fab_vfrac[idx],
                                           afrac,
                                           as_fine,
                                           dm_as_fine,
                                           level_mask[mfi]);


                if (as_crse) {
                    fr_as_crse[idx]->CrseAdd(mfi,
                    {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])},
                                             dx, dt,
                                             *fab_vfrac[idx],
                                             afrac, RunOn::Cpu);
                }

                if (as_fine) {
                    fr_as_fine[idx]->FineAdd(mfi,
                    {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])},
                                             dx, dt,
                                             *fab_vfrac[idx],
                                             afrac,
                                             dm_as_fine, RunOn::Cpu);
                }

            } else {

#endif
                // calculate divergence

                istate.calc_divergence(box,
                                       fluxes[idx],
                                       dU,
                                       dx,
                                       dt);

                if (fr_as_crse[idx]) {
                    fr_as_crse[idx]->CrseAdd(mfi, {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])}, dx, dt, RunOn::Cpu);
                }

                if (fr_as_fine[idx]) {
                    fr_as_fine[idx]->FineAdd(mfi, {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])}, dx, dt, RunOn::Cpu);
                }

#ifdef AMREX_USE_EB
            }
#endif

        }

        // update the 'new' state held by AmrLevel for all of the states
        // this is performed according to new = new + dt*update

        for (int idx=0; idx<gd.num_solve_state; ++idx) {

            if (active[idx] == FabType::covered) continue;

            FArrayBox& s1 = local_new[idx][mfi];
            FArrayBox& s2 = (*global_new[idx])[mfi];


            s2.linComb(s1, box, 0, update[idx][mfi], box, 0, 1.0, 1.0, box, 0, s1.nComp());

        }

        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
