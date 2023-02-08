#include "MFP.H"
#include "MFP_state.H"
#include "MFP_bc.H"
#include <AMReX_ParmParse.H>

void MFP::variableSetUp()
{
    BL_PROFILE("MFP::variableSetUp");
    // read in all of the user defined parameters
    read_params();

    ParmParse pgm("geometry");
    Vector<int> periodic(AMREX_SPACEDIM,0);
    pgm.queryarr("is_periodic",periodic);

    for (auto& state : states) {
        state->variable_setup(periodic);
    }

    update_config_post_data_instantiation();

    //===================
    // cost state

    Cost_Idx = desc_lst.size();

    // just generate a bc
    BCRec bc;
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        bc.setLo(i,BCType::foextrap);
        bc.setHi(i,BCType::foextrap);
    }

    desc_lst.addDescriptor(Cost_Idx, IndexType::TheCellType(),
                           StateDescriptor::Point, 0, 1, &pc_interp);
    desc_lst.setComponent(Cost_Idx, 0, "Cost", bc,
                          NullFillBC());


    StateDescriptor::setBndryFuncThreadSafety(true);

    return;
}

void MFP::variableCleanUp()
{
    desc_lst.clear();

    for (auto& state : states) {
        state.reset();
    }

    states.clear();

#ifdef AMREX_USE_EB
    eb_data.clear();
#endif

}

void MFP::build_eb() {
#ifdef AMREX_USE_EB
    BL_PROFILE("MFP::build_eb()");

    // make sure dx == dy == dz
    const Real* dx = geom.CellSize();
    for (int i = 1; i < AMREX_SPACEDIM; ++i) {
        if (std::abs(dx[0] - dx[i]) > 1.e-12 * dx[0]) {
            amrex::Abort("MFP: must have dx == dy == dz\n");
        }
    }


    // get the information for each embedded boundary
    const Vector<int> ngrow = {m_eb_basic_grow_cells,m_eb_volume_grow_cells,m_eb_full_grow_cells};

    eb_data.resize(level+1);
    eb_data[level].resize(states.size());

    for (const auto& state : states) {

        EBData& eb = get_eb_data(state->global_idx);

        eb.ebfactory = makeEBFabFactory (state->eb2_index,
                                                     geom,
                                                     grids,
                                                     dmap,
                                                     ngrow,
                                                     EBSupport::full);


        const auto& flags_orig = eb.ebfactory->getMultiEBCellFlagFab();
        eb.flags.define(grids, dmap, 1, flags_orig.n_grow);
        eb.flags.ParallelCopy(flags_orig,0,0,1,flags_orig.n_grow,flags_orig.n_grow,geom.periodicity());

        const auto& vfrac_original = eb.ebfactory->getVolFrac();
        eb.volfrac.define(grids, dmap, 1, vfrac_original.n_grow);
        eb.volfrac.ParallelCopy(vfrac_original,0,0,1,vfrac_original.n_grow,vfrac_original.n_grow,geom.periodicity());

        eb.bndrycent = &(eb.ebfactory->getBndryCent());
        eb.bndrynorm = &(eb.ebfactory->getBndryNormal());
        eb.areafrac = eb.ebfactory->getAreaFrac();
        eb.facecent = eb.ebfactory->getFaceCent();

        // different boundary conditions
        eb.bndryidx.define(grids, dmap, 1, m_eb_basic_grow_cells, flags_orig);
        eb.bndryidx.setVal(-1.0);

        // update ghost cells
        auto& flags = eb.flags;
        auto& vfrac = eb.volfrac;

        for (MFIter mfi(vfrac); mfi.isValid(); ++mfi){

            //            plot_FAB_2d(flags[mfi], "flags before", false);
            //            plot_FAB_2d(vfrac[mfi], 0, "vfrac before", false,false);

            state->update_eb_vfrac(geom,vfrac[mfi]);
            state->update_eb_flags(geom,flags[mfi]);

            //            plot_FAB_2d(flags[mfi], "flags after", false);
            //            plot_FAB_2d(vfrac[mfi], 0, "vfrac after", false,true);
        }
    }



    // loop over all of the individual EB definitions and put an id into bndryidx wherever there is
    // a boundary
    for (const auto &eb : eb_def) {

        std::unique_ptr<EBFArrayBoxFactory> bc_ebfactory = makeEBFabFactory(eb.index_space,
                                                                            geom,
                                                                            grids,
                                                                            dmap,
                                                                            ngrow,
                                                                            EBSupport::full);

        const FabArray<EBCellFlagFab>& bc_flags = bc_ebfactory->getMultiEBCellFlagFab();

        for (const auto& si : eb.states) {
            EBData& ebd = get_eb_data(si.first);
            MultiCutFab& bndryidx = ebd.bndryidx;
            const FabArray<EBCellFlagFab>& state_flags = ebd.flags;

            // fill in all entries of the bndryidx with the appropriate index
            for (MFIter mfi(bc_flags); mfi.isValid(); ++mfi){
                const Box& box = mfi.growntilebox();

                const EBCellFlagFab& bc_flag = bc_flags[mfi];
                const EBCellFlagFab& state_flag = state_flags[mfi];



                // check if there is anything to do in this box
                if ((bc_flag.getType(box) != FabType::singlevalued) || (state_flag.getType(box) != FabType::singlevalued))
                    continue;

                const Dim3 lo = amrex::lbound(box);
                const Dim3 hi = amrex::ubound(box);

                CutFab& fab_idx = bndryidx[mfi];

                const Array4<const EBCellFlag>& bc_flag4 = bc_flag.array();
                const Array4<const EBCellFlag>& state_flag4 = state_flag.array();
                const Array4<Real>& fab_idx4 = fab_idx.array();

                for     (int k = lo.z; k <= hi.z; ++k) {
                    for   (int j = lo.y; j <= hi.y; ++j) {
                        AMREX_PRAGMA_SIMD
                                for (int i = lo.x; i <= hi.x; ++i) {

                            if (bc_flag4(i,j,k).isSingleValued() && state_flag4(i,j,k).isSingleValued()) {
                                fab_idx4(i,j,k) = si.second; // the index into the states list of eb boundary conditions
                            }
                        }
                    }
                }
            }
        }
    }

    level_mask.clear();
    level_mask.define(grids, dmap, 1, 2);
    level_mask.BuildMask(geom.Domain(), geom.periodicity(), level_mask_covered,
                         level_mask_notcovered, level_mask_physbnd,
                         level_mask_interior);
#endif
    return;
}

