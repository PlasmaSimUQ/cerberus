#include "MFP_CTU_hydro.H"
#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"
#include "Eigen"
#include "Dense"
#include "MFP_diagnostics.H"

std::string HydroCTU::tag = "CTU";
bool HydroCTU::registered = GetActionFactory().Register(HydroCTU::tag, ActionBuilder<HydroCTU>);

HydroCTU::HydroCTU(){}
HydroCTU::~HydroCTU(){}

HydroCTU::HydroCTU(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    do_CTU = def.get_or("corner_transport", true);

    const sol::table state_names = def["states"];

    for (const auto& key_value_pair : state_names) {
        std::string state_name = key_value_pair.second.as<std::string>();
        State& istate = MFP::get_state(state_name);

        switch (istate.get_classification()) {
        case State::StateClassification::Eulerian:
            states.push_back(static_cast<EulerianState*>(&istate));
            state_indexes.push_back(istate.global_idx);
            break;
        default:
            Abort("An invalid state has been defined for the CTU source "+name);
        }
    }

    return;
}

void HydroCTU::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("HydroCTU::get_data");

    Vector<Array<int,2>> options(states.size());

    for (size_t i=0; i<states.size();++i) {
        options[i] = {states[i]->global_idx, 1};
    }

    Action::get_data(mfp, options, update, time);

}

void HydroCTU::calc_spatial_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt, const Real flux_register_scale)
{
    BL_PROFILE("CTU::calc_spatial_derivative");

    const Geometry& geom = mfp->Geom();
    const int level = mfp->get_level();
    const int finest_level = mfp->get_parent()->finestLevel();


#ifdef AMREX_USE_EB
    constexpr int num_grow_eb = 2;
#else
    constexpr int num_grow_eb = 0;
#endif


    const size_t n_states = states.size();
    Vector<int> data_indexes;
    Vector<EulerianState*> data_states;

    for (const int& global_idx : state_indexes) {
        EulerianState &istate = EulerianState::get_state_global(global_idx);
        data_indexes.push_back(istate.data_idx);
        data_states.push_back(&istate);
    }


    // ==========================================================================
    // all of level storage

    Vector<MultiFab*> local_old(n_states); // local data that covers the whole level
    for (size_t i=0; i<n_states; ++i) {
        local_old[i] = &update[data_indexes[i]].U;
        update[data_indexes[i]].dU_status = UpdateData::Status::Changed;
    }
    Vector<FabType> active(n_states); // flag for calculation

#ifdef AMREX_USE_EB
    Vector<EBFluxRegister*> fr_as_crse(n_states, nullptr);
    Vector<EBFluxRegister*> fr_as_fine(n_states, nullptr);
#else
    Vector<YAFluxRegister*> fr_as_crse(n_states, nullptr);
    Vector<YAFluxRegister*> fr_as_fine(n_states, nullptr);
#endif

    for (int idx=0; idx<n_states; ++idx) {
        EulerianState &istate = *data_states[idx];
        const int data_idx = data_indexes[idx];

        if (istate.reflux && level < finest_level) {
            MFP& fine_level = mfp->getLevel(level + 1);
            fr_as_crse[idx] = &fine_level.flux_reg[data_idx];
        }

        if (istate.reflux && level > 0) {
            fr_as_fine[idx] = &mfp->flux_reg[data_idx];
        }
    }

    // ==========================================================================
    // per state storage

    Vector<FArrayBox*> conserved(n_states);
    Vector<FArrayBox> primitives(n_states);
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> R_lo(n_states), R_hi(n_states);
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> fluxes(n_states);

    int nu; // per state size

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

#ifdef AMREX_USE_EB
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> wall_fluxes(n_states);
    Vector<const EBCellFlagFab*> fab_flags(n_states);
    Vector<const FArrayBox*> fab_vfrac(n_states);
#endif

    // ==========================================================================
    // iterate over all of the FABs within the level performing reconstruction, flux calculation,
    // and updating the cell-centred data according to the resulting divergence


    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();

        // region over which to perform reconstruction for face values
        const Box rbox = amrex::grow(box, 1 + num_grow_eb);

        // ==========================================================================
        // 1. iterate over all states to set-up the data required for flux calculation
        for (int idx=0; idx<n_states; ++idx) {
            EulerianState &istate = *data_states[idx];

            // get a pointer to the conserved quantities
            conserved[idx] = &(*local_old[idx])[mfi];

#ifdef AMREX_USE_EB
            EBData& eb = mfp->get_eb_data(istate.global_idx);
            // get the EB data required for later calls
            const EBCellFlagFab& flag = eb.flags[mfi]; fab_flags[idx] = &flag;
            const FArrayBox& vfrac = eb.volfrac[mfi]; fab_vfrac[idx] = &vfrac;

            // check if this box is active
            active[idx] = flag.getType(rbox);
            if (active[idx] == FabType::covered)
                continue;
#else
            active[idx] = FabType::regular;
#endif

            // region over which to get cell centered primitives for reconstruction
            const Box pbox = amrex::grow(box, istate.num_grow + num_grow_eb);

            // ===============================================
            // 1.1 Calculate primitive values within each cell

            FArrayBox& cons = (*local_old[idx])[mfi];
            FArrayBox& prim = primitives[idx];

            istate.calc_primitives(pbox,
                                   cons,
                                   prim,
                                   dx,
                                   time,
                                   prob_lo
                       #ifdef AMREX_USE_EB
                                   ,vfrac
                       #endif
                                   );

//            plot_FAB_2d(prim, +HydroDef::PrimIdx::NUM, "prim", false, true);


            // fill in any cells that need special boundary values
            istate.update_boundary_cells(pbox,
                                         geom,
                                         prim,
                             #ifdef AMREX_USE_EB
                                         vfrac,
                             #endif
                                         time);



            // =======================================
            // 1.2 Calculate reconstructed face values

            // each cell has a hi and lo side in each direction

            // calculate the reconstructed face values
            istate.calc_reconstruction(rbox,
                                       prim,
                                       R_lo[idx],
                                       R_hi[idx]
                           #ifdef AMREX_USE_EB
                                       ,flag
                                       ,vfrac
                           #endif
                                       );

            // ===================================================================
            // update the face values to time t+1/2 based on the local wave speeds

            // TODO: this step currently assumes a single speed for all components and
            // should be updated to calculate the correct characteristic speeds
            // for each component individually
            if (do_CTU) {
                istate.calc_time_averaged_faces(rbox,
                                                prim,
                                                R_lo[idx],
                                                R_hi[idx],
                                #ifdef AMREX_USE_EB
                                                flag,
                                #endif
                                                dx,
                                                dt);
            }

        }

        // 3.1 Setup for flux calculation

        // resize the flux arrays before any get used
        for (int idx=0; idx<n_states; ++idx) {
            EulerianState &istate = *data_states[idx];

            if (active[idx] == FabType::covered)
                continue;

            nu = local_old[idx]->nComp();
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                // get a node centered box in direction d that encloses the original box
                Box fbox = surroundingNodes(box, d);

                if (do_CTU || (istate.is_viscous() && num_grow_eb > 0)) {
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
        for (int idx=0; idx<n_states; ++idx) {
            EulerianState &istate = *data_states[idx];

            if (active[idx] == FabType::covered) {
                continue;
            }

            istate.update_face_prim(box,
                                    geom,
                                    R_lo[idx],
                                    R_hi[idx],
                        #ifdef AMREX_USE_EB
                                    *fab_flags[idx],
                        #endif
                                    time);
        }

        //---
        for (int idx=0; idx<n_states; ++idx) {
            EulerianState &istate = *data_states[idx];

            if (active[idx] == FabType::covered) continue;

            istate.calc_fluxes(box,
                               *conserved[idx],
                               R_lo[idx],
                               R_hi[idx],
                               fluxes[idx],
                   #ifdef AMREX_USE_EB
                               *fab_flags[idx],
                   #endif
                               dx, dt);

#if AMREX_SPACEDIM > 1
            if (do_CTU) {
                // now that we have done the fluxes, do we need to correct them???
                // correction is done according to the following steps:
                // 1. calculate the fluxes (above)
                // 2. modify the reconstructed face states using the flux information
                // 3. recalculate the fluxes with the updated face values

                istate.correct_face_prim(grow(box,num_grow_eb),
                                         R_lo[idx],
                                         R_hi[idx],
                                         fluxes[idx],
                         #ifdef AMREX_USE_EB
                                         *fab_flags[idx],
                         #endif
                                         dx, dt);

                // following the update of the face values we need to update any boundary conditions

                istate.update_face_prim(box,
                                        geom,
                                        R_lo[idx],
                                        R_hi[idx],
                        #ifdef AMREX_USE_EB
                                        *fab_flags[idx],
                        #endif
                                        time,
                                        true);
            }
#endif
        }

        //---
#if AMREX_SPACEDIM > 1
        if (do_CTU) {
            for (int idx=0; idx<n_states; ++idx) {
                EulerianState &istate = *data_states[idx];
                if (active[idx] == FabType::covered) continue;

                // recalculate the fluxes
                istate.calc_fluxes(box,
                                   *conserved[idx],
                                   R_lo[idx],
                                   R_hi[idx],
                                   fluxes[idx],
                   #ifdef AMREX_USE_EB
                                   *fab_flags[idx],
                   #endif
                                   dx, dt);
            }
        }
#endif

        //---
        for (int idx=0; idx<n_states; ++idx) {
            EulerianState &istate = *data_states[idx];
            const int data_idx = data_indexes[idx];


            if (active[idx] == FabType::covered) continue;

            if (istate.is_viscous()) {

                // now calculate any viscous fluxes

                istate.calc_viscous_fluxes(grow(box,num_grow_eb),
                                           fluxes[idx],
                                           primitives[idx],
                               #ifdef AMREX_USE_EB
                                           *fab_flags[idx],
                               #endif
                                           dx);
            }


            // given all of the fluxes calculate the update to the cell centred values
            const int as_crse = (fr_as_crse[idx] != nullptr);
            const int as_fine = (fr_as_fine[idx] != nullptr);

            nu = local_old[idx]->nComp();

            FArrayBox& du = update[data_idx].dU[mfi];


#ifdef AMREX_USE_EB

            EBData& eb = mfp->get_eb_data(istate.global_idx);

            Array<const FArrayBox*,AMREX_SPACEDIM> afrac, fcent;

            if (active[idx] != FabType::regular) {

                const FArrayBox &bcent = (*eb.bndrycent)[mfi];
                const FArrayBox &bnorm = (*eb.bndrynorm)[mfi];

                CutFab& bc_idx = eb.bndryidx[mfi];

                afrac = {AMREX_D_DECL(&(*eb.areafrac[0])[mfi], &(*eb.areafrac[1])[mfi], &(*eb.areafrac[2])[mfi])};


                // calculate the flux through cut cell faces

                istate.calc_wall_fluxes(box,
                                        primitives[idx],
                                        wall_fluxes[idx],
                                        *fab_flags[idx],
                                        bc_idx,
                                        bcent,
                                        bnorm,
                                        afrac,
                                        dx,
                                        dt);

//                plot_FAB_2d(primitives[idx], +HydroDef::PrimIdx::Prs, "prim", false, false);

//                plot_FAB_2d(wall_fluxes[idx][0], +HydroDef::ConsIdx::Eden, "flux nrg x", false, false);
//                plot_FAB_2d(wall_fluxes[idx][1], +HydroDef::ConsIdx::Eden, "flux nrg y", false, true);

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

                fcent = {AMREX_D_DECL(&(*eb.facecent[0])[mfi], &(*eb.facecent[1])[mfi], &(*eb.facecent[2])[mfi])};

                istate.calc_eb_divergence(box,
                                          *conserved[idx],
                                          fluxes[idx],
                                          wall_fluxes[idx],
                                          du,
                                          *fab_flags[idx],
                                          *fab_vfrac[idx],
                                          afrac,
                                          fcent,
                                          as_crse,
                                          as_fine,
                                          p_rrflag_as_crse,
                                          mfp->level_mask[mfi],
                                          p_drho_as_crse,
                                          dm_as_fine,
                                          dx,
                                          dt);

                istate.merge_cells(box,
                                   *conserved[idx],
                                   du,
                                   *fab_flags[idx],
                                   *fab_vfrac[idx],
                                   afrac,
                                   as_fine,
                                   dm_as_fine,
                                   mfp->level_mask[mfi]);


                if (as_crse) {
                    fr_as_crse[idx]->CrseAdd(mfi,
                    {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])},
                                             dx, flux_register_scale,
                                             *fab_vfrac[idx],
                                             afrac, RunOn::Cpu);
                }

                if (as_fine) {
                    fr_as_fine[idx]->FineAdd(mfi,
                    {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])},
                                             dx, flux_register_scale,
                                             *fab_vfrac[idx],
                                             afrac,
                                             dm_as_fine, RunOn::Cpu);
                }

            } else {

#endif

                // calculate divergence

                istate.calc_divergence(box,
                                       fluxes[idx],
                                       du,
                                       dx,
                                       dt);

                if (fr_as_crse[idx]) {
                    fr_as_crse[idx]->CrseAdd(mfi, {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])}, dx, flux_register_scale, RunOn::Cpu);
                }

                if (fr_as_fine[idx]) {
                    fr_as_fine[idx]->FineAdd(mfi, {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])}, dx, flux_register_scale, RunOn::Cpu);
                }

#ifdef AMREX_USE_EB
            }
#endif

        }

        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }

}

