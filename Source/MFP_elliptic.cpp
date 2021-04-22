#include "MFP.H"

//#include <AMReX_MLPoisson.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_LO_BCTYPES.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_MLEBABecLap.H>
#else
#include <AMReX_MLABecLaplacian.H>
#endif

using Location = MLLinOp::Location;

void set_field_bcs(Array<LinOpBCType, AMREX_SPACEDIM> &bc_lo,
                   Array<LinOpBCType, AMREX_SPACEDIM> &bc_hi,
                   const BCRec &bc,
                   Geometry &geom)
{
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        if (geom.isPeriodic(d)) {
            bc_lo[d] = LinOpBCType::Periodic;
            bc_hi[d] = LinOpBCType::Periodic;
        } else {
            switch (bc.lo(d)) {

            case PhysBCType::interior :
                bc_lo[d] = LinOpBCType::Periodic;
                break;

            case PhysBCType::inflow :
                bc_lo[d] = LinOpBCType::Dirichlet;
                break;

            case PhysBCType::outflow :
                bc_lo[d] = LinOpBCType::Dirichlet;
                break;

            case PhysBCType::symmetry :
                bc_lo[d] = LinOpBCType::Neumann;
                break;

            case PhysBCType::slipwall :
                bc_lo[d] = LinOpBCType::reflect_odd;
                break;

            case PhysBCType::noslipwall :
                bc_lo[d] = LinOpBCType::reflect_odd;
                break;
            }


            switch (bc.hi(d)) {

            case PhysBCType::interior :
                bc_hi[d] = LinOpBCType::Periodic;
                break;

            case PhysBCType::inflow :
                bc_hi[d] = LinOpBCType::Dirichlet;
                break;

            case PhysBCType::outflow :
                bc_hi[d] = LinOpBCType::Dirichlet;
                break;

            case PhysBCType::symmetry :
                bc_hi[d] = LinOpBCType::Neumann;
                break;

            case PhysBCType::slipwall :
                bc_hi[d] = LinOpBCType::reflect_odd;
                break;

            case PhysBCType::noslipwall :
                bc_hi[d] = LinOpBCType::reflect_odd;
                break;
            }

        }
    }
}

void MFP::solve_static_fields(const Real time)
{
    BL_PROFILE("MFP::solve_static_fields()");
#if AMREX_SPACEDIM > 1
    for (auto &istate : gd.states) {

        if (!((istate->get_type() == +StateType::isField) && (!istate->is_transported())))
            continue;

        const int state_idx = istate->global_idx;

        Real cd_scale = (gd.Larmor)/(gd.Debye*gd.Debye*gd.lightspeed);
        Real J_scale  = (gd.Larmor)/(gd.Debye*gd.Debye*gd.lightspeed*gd.lightspeed);

        // assume zero charge density
        Amr& amr = *parent;
        //            const Vector<Geometry>& amr_geom = amr.Geom();
        const Vector<BoxArray>& amr_grids = amr.boxArray(0, amr.finestLevel());
        const Vector<DistributionMapping>& amr_dmap = amr.DistributionMap(0, amr.finestLevel());
        const int nlevels = amr_grids.size();

        const Vector<Geometry>& amr_geom = amr.Geom(0, amr.finestLevel());

#ifdef AMREX_USE_EB
        Vector<EBFArrayBoxFactory const*> amr_factory(nlevels);

        for (int ilev=0; ilev<nlevels; ++ilev) {
            amr_factory[ilev] = getEBData(ilev,state_idx).ebfactory.get();
        }
#endif

        Vector<MultiFab> defined_phi(nlevels);
        Vector<MultiFab> defined_A(nlevels);
        Vector<MultiFab> defined_charge(nlevels);
        Vector<MultiFab> defined_current(nlevels);


        const int ngrow = 1;

        bool update_D = false;
        bool update_B = false;

        // go through level-by-level and calculate the charge and current density and electric potential
        for (int ilev = 0; ilev < nlevels; ++ilev) {

            defined_phi[ilev].define(amr_grids[ilev], amr_dmap[ilev], 1, ngrow, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
            defined_phi[ilev].setVal(0.0);

            defined_A[ilev].define(amr_grids[ilev], amr_dmap[ilev], 3, ngrow, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
            defined_A[ilev].setVal(0.0);

            defined_charge[ilev].define(amr_grids[ilev], amr_dmap[ilev], 1, ngrow, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
            defined_charge[ilev].setVal(0.0);

            defined_current[ilev].define(amr_grids[ilev], amr_dmap[ilev], 3, ngrow, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
            defined_current[ilev].setVal(0.0);

            AmrLevel& ilevel = amr.getLevel(ilev);
            EB_OPTIONAL(Vector<EBData>& eb_data = gd.eb_data[ilev];)

                    // get all of the source contributions
                    for (const auto& src_idx : istate->associated_sources) {
                const auto &ode = gd.ode_source_terms[src_idx.first];
                const auto &src = ode->sources[src_idx.second];

                if (!src->has_charge_density()) continue;

                // now calculate the charge and current density

                Vector<FArrayBox*> src_data(src->offsets.size());
                EB_OPTIONAL(Vector<const EBCellFlagFab*> src_flags(src->offsets.size());)
                        FArrayBox local_cd;
                FArrayBox local_J;


                // iterate over data
                for (MFIter mfi(ilevel.get_data(0, time)); mfi.isValid(); ++mfi) {
                    const Box& box = mfi.tilebox();
                    local_cd.resize(box);
                    local_J.resize(box,3);

                    // get lists of the data for each state included in the src
                    for (const auto& idx : src->offsets) {
                        src_data[idx.local] = &ilevel.get_data(idx.global, time)[mfi];
                        EB_OPTIONAL(src_flags[idx.local] = &(eb_data[idx.global].flags)[mfi];)
                    }

                    src->calc_charge_density(box,
                                             geom.ProbLo(),
                                             geom.CellSize(),
                                             parent->cumTime(),
                                             src_data,
                                             local_cd,
                                             local_J
                                             EB_OPTIONAL(,src_flags)
                                             );

                    // add the scaled charge density to the accumulator
                    defined_charge[ilev][mfi].plus(local_cd);
                    defined_current[ilev][mfi].plus(local_J);
                }
            }


            // get the embedded boundary contributions
#ifdef AMREX_USE_EB
            std::map<BoundaryEB::EBType,FArrayBox*> bcs_data;

            for (MFIter mfi(defined_phi[ilev]); mfi.isValid(); ++mfi) {

                const Box& box = mfi.validbox();
                const Box cbox = grow(box,ngrow);

                // specify here what EB bcs we are looking to get info from
                bcs_data[BoundaryEB::EBType::ScalarPotential] = &defined_phi[ilev][mfi];
                bcs_data[BoundaryEB::EBType::VectorPotential] = &defined_A[ilev][mfi];
                //                bcs_data[BoundaryEB::EBType::SurfaceCharge] = &defined_charge[ilev][mfi];
                //                bcs_data[BoundaryEB::EBType::SurfaceCurrent] = &defined_current[ilev][mfi];
                //

                const EBCellFlagFab& flag = getEBData(ilev, state_idx).flags[mfi];

                if (flag.getType(cbox) != FabType::singlevalued) continue;

                const EBData& ebd = getEBData(ilev, state_idx);

                const FArrayBox &bnorm = (*ebd.bndrynorm)[mfi];
                const FArrayBox &bcent = (*ebd.bndrycent)[mfi];

                const CutFab& bc_idx = ebd.bndryidx[mfi];

                istate->get_wall_value(cbox,
                                       bcs_data,
                                       flag,
                                       bc_idx,
                                       bcent,
                                       bnorm,
                                       time,
                                       ilevel.Geom().CellSize(),
                                       ilevel.Geom().ProbLo());
            }

#endif



            // scale the charge and current density appropriately
            defined_charge[ilev].mult(cd_scale, 0);
            defined_current[ilev].mult(J_scale, 0);

            // check if we need to update D
            if (defined_phi[ilev].norminf(0, ngrow, false, true) > 0.0) {
                update_D = true;
            }

            if (!update_D) {
                if (defined_charge[ilev].norminf(0, ngrow, false, true) > 0.0) {
                    update_D = true;
                }
            }

            // check if we need to update B
            for (int d = 0; d<3; ++d) {
                if (defined_A[ilev].norminf (d, ngrow, false, true) > 0.0) {
                    update_B = true;
                    break;
                }

                if (defined_current[ilev].norminf (d, ngrow, false, true) > 0.0) {
                    update_B = true;
                    break;
                }
            }

            //            plot_FAB_2d(defined_phi[ilev],0, "phi", false, true);

//            plot_FAB_2d(defined_A[ilev],0,0, "Ax", false, true);
//            plot_FAB_2d(defined_A[ilev],1,0, "Ay", false, true);
//            plot_FAB_2d(defined_A[ilev],2,0, "Az", false, true);

        }

        if (!update_D && !update_B) continue;

        ParallelDescriptor::Barrier();

        // now that we have the charge density over all levels do the actual solve
        const int D_idx = istate->get_cons_D_idx();
        const int B_idx = istate->get_cons_B_idx();
        const int phi_idx = istate->get_cons_phi_idx();

        //---------------------------------------------------------------------------------------



        // Linear Solver
        LPInfo lp_info;

#ifdef AMREX_USE_EB
        MLEBABecLap mlabec(amr_geom,
                           amr_grids,
                           amr_dmap,
                           lp_info,
                           amr_factory
                           );
#else
        MLABecLaplacian mlabec(amr_geom,
                               amr_grids,
                               amr_dmap,
                               lp_info);
#endif

        MLMG mlmg(mlabec);

        // relative and absolute tolerances for linear solve
        const Real tol_rel = 1.e-10;
        const Real tol_abs = 0.0;

        mlmg.setVerbose(gd.linear_solver_verbosity);

        // define array of LinOpBCType for domain boundary conditions
        Array<LinOpBCType, AMREX_SPACEDIM> bc_lo;
        Array<LinOpBCType, AMREX_SPACEDIM> bc_hi;



        Vector<const MultiFab*> rhs_ptr(nlevels, nullptr);
        Vector<Array<MultiFab,AMREX_SPACEDIM>> bcoef(nlevels);
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            // set BCoef (face centered coefficients)
            Array<MultiFab,AMREX_SPACEDIM>& bcf = bcoef[ilev];
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                bcf[idim].define(amrex::convert(amr_grids[ilev],IntVect::TheDimensionVector(idim)),
                                 amr_dmap[ilev], 1, 0, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
            }
        }



        //=====================================================================
        // phi (for D field)

        if (update_D) {

            // Boundary of the whole domain. This functions must be called,
            // and must be called before other bc functions.
            const BCRec &bc = istate->boundary_conditions.phys_fill_bc[D_idx];
            set_field_bcs(bc_lo, bc_hi, bc, geom);
            mlabec.setDomainBC(bc_lo,bc_hi);

            // scaling factors; these multiply ACoef and BCoef (see below)
            mlabec.setScalars(0.0, 1.0);

            Vector<MultiFab*> phi_ptr(nlevels, nullptr);

            for (int ilev = 0; ilev < nlevels; ++ilev) {

                MultiFab &phi_mf = defined_phi[ilev];
                phi_ptr[ilev] = &phi_mf;

                // see AMReX_MLLinOp.H for an explanation
                mlabec.setLevelBC(ilev, nullptr);

                // operator looks like (ACoef - div BCoef grad) phi = rhs
                mlabec.setACoeffs(ilev, 0.0);

                // think of this beta as the "BCoef" associated with an EB face (typically 1.0)
                EB_OPTIONAL(mlabec.setEBDirichlet(ilev, phi_mf, 1.0);)

                        // set BCoef (face centered coefficients)
                        Array<MultiFab,AMREX_SPACEDIM>& bcf = bcoef[ilev];
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    bcf[idim].setVal(1.0);
                }

                mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcf));

                // build the rhs pointers
                rhs_ptr[ilev] = &defined_charge[ilev];

            }

//            plot_FAB_2d(defined_phi[nlevels-1],0,1, "phi", false, false);
//            plot_FAB_2d(defined_charge[nlevels-1],0,1, "cd", false, true);


            // Solve linear system
            mlmg.solve(phi_ptr, rhs_ptr, tol_rel, tol_abs);

            // save phi
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                AmrLevel& ilevel = amr.getLevel(ilev);
                MultiFab& ifield = ilevel.get_data(state_idx, time);

                MultiFab::Copy(ifield, defined_phi[ilev], 0, phi_idx, 1, 0);
            }

            // Get fluxes from solver
            mlmg.getFluxes( {GetVecOfArrOfPtrs(bcoef)});  // reuse bcoef

            // solve for D field
            for (int ilev = 0; ilev < nlevels; ++ilev){
                AmrLevel& ilevel = amr.getLevel(ilev);
                MultiFab& ifield = ilevel.get_data(state_idx, time);

#ifdef AMREX_USE_EB
                EB_average_face_to_cellcenter(ifield, D_idx, GetArrOfConstPtrs(bcoef[ilev]));
#else
                average_face_to_cellcenter(ifield, D_idx, amrex::GetArrOfConstPtrs(bcoef[ilev]));
#endif
            }
        }

        //=====================================================================
        // A (for B field)
        if (update_B) {

            // Boundary of the whole domain. This functions must be called,
            // and must be called before other bc functions.
            const BCRec &bc = istate->boundary_conditions.phys_fill_bc[B_idx];
            set_field_bcs(bc_lo, bc_hi, bc, geom);
            mlabec.setDomainBC(bc_lo,bc_hi);

            // scaling factors; these multiply ACoef and BCoef (see below)
            mlabec.setScalars(0.0, 1.0);

            Vector<Array<MultiFab,3>> A(nlevels);
            Vector<MultiFab*> A_ptr(nlevels, nullptr);

            Vector<MultiFab> current_alias(nlevels);

            // zero out the B field value
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                AmrLevel& ilevel = amr.getLevel(ilev);
                MultiFab& ifield = ilevel.get_data(state_idx, time);
                ifield.setVal(0.0, B_idx, 3);

                mlabec.setLevelBC(ilev, nullptr);
                mlabec.setACoeffs(ilev, 0.0);
            }

            for (int d = 0; d<3; ++d) {
                for (int ilev = 0; ilev < nlevels; ++ilev) {

                    MultiFab &A_mf = A[ilev][d];
                    A_mf = MultiFab(defined_A[ilev], MakeType::make_alias, d, 1);
                    A_ptr[ilev] = &A_mf;

                    EB_OPTIONAL(mlabec.setEBDirichlet(ilev, A_mf, 1.0);)

                            // set BCoef (face centered coefficients)
                            Array<MultiFab,AMREX_SPACEDIM>& bcf = bcoef[ilev];
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                        bcf[idim].setVal(1.0);
                    }

                    mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcf));

                    // build the rhs pointers

                    current_alias[ilev] = MultiFab(defined_current[ilev], MakeType::make_alias, d, 1);

                    rhs_ptr[ilev] = &current_alias[ilev];

                }

//                plot_FAB_2d(*(A_ptr[nlevels-1]),0,0, "A"+num2str(d), false, false);
                //                plot_FAB_2d(*(rhs_ptr[nlevels-1]),0, "J"+num2str(d), false, true);


                // Solve linear system
                mlmg.solve(A_ptr, rhs_ptr, tol_rel, tol_abs);

//                plot_FAB_2d(*(A_ptr[nlevels-1]),0,0, "A"+num2str(d), false, true);

                // Get fluxes from solver
                mlmg.getFluxes( {GetVecOfArrOfPtrs(bcoef)});  // reuse bcoef


                // save A for debugging
                if (false) {
                    for (int ilev = 0; ilev < nlevels; ++ilev) {
                        AmrLevel& ilevel = amr.getLevel(ilev);
                        MultiFab& ifield = ilevel.get_data(state_idx, time);

                        MultiFab::Copy(ifield, A[ilev][d], 0, B_idx+d, 1, 0); // magnetic vector potential

                        //                        MultiFab::Copy(ifield, defined_current[ilev], d, B_idx+d, 1, 0); // current
                    }
                } else {
                    // solve for B field
                    for (int ilev = 0; ilev < nlevels; ++ilev) {
                        AmrLevel& ilevel = amr.getLevel(ilev);
                        MultiFab& ifield = ilevel.get_data(state_idx, time);

                        // calculate the gradients
                        MultiFab gradients(amr_grids[ilev], amr_dmap[ilev], 3, 0, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));

#ifdef AMREX_USE_EB
                        EB_average_face_to_cellcenter(gradients, 0, GetArrOfConstPtrs(bcoef[ilev]));
#else
                        average_face_to_cellcenter(gradients, 0, amrex::GetArrOfConstPtrs(bcoef[ilev]));
#endif

                        // we now have -sigma*grad(A) at the cell centre (sigma also known as B in the docs)
                        // iteratively calculate the curl and thus B, i.e. B = curl A
                        switch (d) {
                        case 0:
                            MultiFab::Add(ifield, gradients, 1, +FieldState::ConsIdx::Bz, 1, 0); //dAx/dy -> Bz
#if AMREX_SPACEDIM == 3
                            MultiFab::Subtract(ifield, gradients, 2, +FieldState::ConsIdx::By, 1, 0); //dAx/dz -> By
#endif
                            break;
                        case 1:
                            MultiFab::Subtract(ifield, gradients, 0, +FieldState::ConsIdx::Bz, 1, 0); //dAy/dx -> Bz
#if AMREX_SPACEDIM == 3
                            MultiFab::Add(ifield, gradients, 2, +FieldState::ConsIdx::By, 1, 0); //dAy/dz -> Bx
#endif
                            break;
                        case 2:
                            MultiFab::Add(ifield, gradients, 0, +FieldState::ConsIdx::By, 1, 0); //dAz/dx -> By
                            MultiFab::Subtract(ifield, gradients, 1, +FieldState::ConsIdx::Bx, 1, 0); //dAz/dy -> Bx

                            break;
                        default:
                            Abort("how did we get here?");
                        }
                    }
                }
            }

            //            plot_FAB_2d(current_alias[nlevels-1],0, "Jz", false, true);
        }
    }
#endif
}



void MFP::correct_dynamic_fields(const Real time)
{
    BL_PROFILE("MFP::correct_dynamic_fields");
    //#if AMREX_SPACEDIM > 1
    //    for (auto &istate : gd.states) {

    //        if (!((istate->has_field()) && (istate->is_transported())))
    //            continue;

    //        if (!istate->project_divergence)
    //            continue;

    //        const int state_idx = istate->global_idx;

    //        Real cd_scale = -(gd.Larmor)/(gd.Debye*gd.Debye*gd.lightspeed);

    //        Amr& amr = *parent;
    //        const Vector<BoxArray>& amr_grids = amr.boxArray(0, amr.finestLevel());
    //        const Vector<DistributionMapping>& amr_dmap = amr.DistributionMap(0, amr.finestLevel());
    //        const Vector<Geometry>& amr_geom = amr.Geom(0, amr.finestLevel());
    //        const int nlevels = amr_grids.size();

    //#ifdef AMREX_USE_EB
    //        Vector<EBFArrayBoxFactory const*> amr_factory(nlevels);

    //        for (int ilev=0; ilev<nlevels; ++ilev) {
    //            amr_factory[ilev] = getEBData(ilev,state_idx).ebfactory.get();
    //        }
    //#endif

    //        Vector<MultiFab> defined_potential(nlevels);
    //        Vector<MultiFab> defined_rhs(nlevels);

    //        //        Vector<MultiFab> charge_density(nlevels);


    //        const int ngrow = 1;


    //        // go through level-by-level and calculate the charge and current density and electric potential
    //        for (int ilev = 0; ilev < nlevels; ++ilev) {

    //            defined_potential[ilev].define(amr_grids[ilev], amr_dmap[ilev], 1, ngrow, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
    //            defined_potential[ilev].setVal(0.0);

    //            defined_rhs[ilev].define(amr_grids[ilev], amr_dmap[ilev], 1, ngrow, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
    //            defined_rhs[ilev].setVal(0.0);

    //            //            charge_density[ilev].define(amr_grids[ilev], amr_dmap[ilev], 1, ngrow, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
    //            //            charge_density[ilev].setVal(0.0);

    //            AmrLevel& ilevel = amr.getLevel(ilev);
    //            EB_OPTIONAL(Vector<EBData>& eb_data = gd.eb_data[ilev];)

    //                    // get all of the source contributions
    //                    for (const auto& src_idx : istate->associated_sources) {
    //                const auto &ode = gd.ode_source_terms[src_idx.first];
    //                const auto &src = ode->sources[src_idx.second];

    //                if (!src->has_charge_density()) continue;

    //                // now calculate the charge and current density

    //                Vector<FArrayBox*> src_data(src->offsets.size());
    //                EB_OPTIONAL(Vector<const EBCellFlagFab*> src_flags(src->offsets.size());)
    //                        FArrayBox local_cd;
    //                FArrayBox local_J; // dummy variable needed for the calc_charge_density call


    //                // iterate over data
    //                for (MFIter mfi(ilevel.get_data(0, time)); mfi.isValid(); ++mfi) {
    //                    const Box& box = mfi.tilebox();
    //                    local_cd.resize(box);
    //                    local_J.resize(box,3);

    //                    // get lists of the data for each state included in the src
    //                    for (const auto& idx : src->offsets) {
    //                        src_data[idx.local] = &ilevel.get_data(idx.global, time)[mfi];
    //                        EB_OPTIONAL(src_flags[idx.local] = &(*eb_data[idx.global].flags)[mfi];)
    //                    }

    //                    src->calc_charge_density(box,
    //                                             geom.ProbLo(),
    //                                             geom.CellSize(),
    //                                             parent->cumTime(),
    //                                             src_data,
    //                                             local_cd,
    //                                             local_J
    //                                             EB_OPTIONAL(,src_flags)
    //                                             );

    //                    // add the scaled charge density to the accumulator
    //                    defined_rhs[ilev][mfi].plus(local_cd);
    //                }
    //            }

    //            /*
    //            // get the embedded boundary contributions
    //            #ifdef AMREX_USE_EB

    //            if (false) {

    //                std::map<BoundaryEB::EBType,FArrayBox*> bcs_data;

    //                for (MFIter mfi(defined_phi[ilev]); mfi.isValid(); ++mfi) {

    //                    const Box& box = mfi.validbox();
    //                    const Box cbox = grow(box,ngrow);

    //                    // specify here what EB bcs we are looking to get info from
    //                    bcs_data[BoundaryEB::EBType::ScalarPotential] = &defined_phi[ilev][mfi];
    //                    bcs_data[BoundaryEB::EBType::VectorPotential] = &defined_A[ilev][mfi];
    //    //                bcs_data[BoundaryEB::EBType::SurfaceCharge] = &defined_charge[ilev][mfi];
    //    //                bcs_data[BoundaryEB::EBType::SurfaceCurrent] = &defined_current[ilev][mfi];
    //                    //

    //                    const EBCellFlagFab& flag = (*getEBData(ilev, state_idx).flags)[mfi];

    //                    if (flag.getType(cbox) != FabType::singlevalued) continue;

    //                    const EBData& ebd = getEBData(ilev, state_idx);

    //                    const FArrayBox &bnorm = (*ebd.bndrynorm)[mfi];

    //                    const CutFab& bc_idx = ebd.bndryidx[mfi];

    //                    istate->get_wall_value(cbox,
    //                                           bcs_data,
    //                                           flag,
    //                                           bc_idx,
    //                                           bnorm,
    //                                           time);


    //                }
    //            }
    //            #endif
    //*/


    //            // scale the charge and current density appropriately
    //            defined_rhs[ilev].mult(cd_scale, 0);


    //            //            MultiFab::Copy(charge_density[ilev], defined_rhs[ilev], 0, 0, 1, ngrow);
    //        }

    //        // now that we have the charge density over all levels do the actual solve
    //        const int D_idx = istate->get_cons_D_idx();
    //        const int B_idx = istate->get_cons_B_idx();
    //        const int phi_idx = istate->get_cons_phi_idx();
    //        const int psi_idx = istate->get_cons_psi_idx();

    //        //---------------------------------------------------------------------------------------



    //        // Linear Solver
    //        LPInfo lp_info;

    //#ifdef AMREX_USE_EB
    //        MLEBABecLap mlabec(amr_geom,
    //                           amr_grids,
    //                           amr_dmap,
    //                           lp_info
    //                           EB_OPTIONAL(,amr_factory)
    //                           );
    //#else
    //        MLABecLaplacian mlabec(amr_geom,
    //                               amr_grids,
    //                               amr_dmap,
    //                               lp_info);
    //#endif

    //        MLMG mlmg(mlabec);


    //        // relative and absolute tolerances for linear solve
    //        const Real tol_rel = 1.e-10;
    //        const Real tol_abs = 0.0;

    //        mlmg.setVerbose(gd.linear_solver_verbosity);

    //        // define array of LinOpBCType for domain boundary conditions
    //        Array<LinOpBCType, AMREX_SPACEDIM> bc_lo;
    //        Array<LinOpBCType, AMREX_SPACEDIM> bc_hi;

    //        Vector<const MultiFab*> rhs_ptr(nlevels, nullptr);
    //        Vector<Array<MultiFab,AMREX_SPACEDIM>> bcoef(nlevels);
    //        for (int ilev = 0; ilev < nlevels; ++ilev) {
    //            // set BCoef (face centered coefficients)
    //            Array<MultiFab,AMREX_SPACEDIM>& bcf = bcoef[ilev];
    //            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    //                bcf[idim].define(amrex::convert(amr_grids[ilev],IntVect::TheDimensionVector(idim)),
    //                                 amr_dmap[ilev], 1, 0, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
    //            }
    //        }


    //        //*******************************************

    //        if (0) {
    //            LPInfo lp_info2;

    //            MLNodeLaplacian matrix(amr_geom,
    //                                   amr_grids,
    //                                   amr_dmap,
    //                                   lp_info2
    //                                   EB_OPTIONAL(,amr_factory)
    //                                   );

    //            Array<LinOpBCType, AMREX_SPACEDIM> lobc;
    //            Array<LinOpBCType, AMREX_SPACEDIM> hibc;
    //            const BCRec &fbc = istate->boundary_conditions.phys_fill_bc[B_idx];
    //            set_field_bcs(lobc,hibc,fbc,geom);
    //            matrix.setDomainBC(lobc,hibc);

    //            matrix.setGaussSeidel(false);
    //            matrix.setHarmonicAverage(false);

    //            Vector<MultiFab> field(nlevels);
    //            Vector<MultiFab> rhs(nlevels);
    //            Vector<MultiFab> rhs2(nlevels);
    //            Vector<MultiFab> S_nd(nlevels);
    //            Vector<MultiFab> S_cc;
    //            Vector<MultiFab*> S_cc_ptr;
    //            Vector<MultiFab> phi(nlevels);
    //            Vector<MultiFab> phi_cc(nlevels);

    //            Vector<MultiFab*> field_ptr(nlevels, nullptr);
    //            Vector<MultiFab*> rhs_ptr1(nlevels, nullptr);
    //            Vector<const MultiFab*> rhs_ptr2(nlevels, nullptr);
    //            Vector<const MultiFab*> S_nd_ptr(nlevels, nullptr);
    //            Vector<MultiFab*> phi_ptr(nlevels, nullptr);

    //            // Cell-centered contributions to RHS
    //            // manufacture zero values if necessary

    //            S_cc.resize(nlevels);
    //            S_cc_ptr.resize(nlevels,nullptr);
    //            for (int ilev = 0; ilev < nlevels; ++ilev) {
    //                S_cc[ilev].define(amr_grids[ilev], amr_dmap[ilev], 1, 1, MFInfo());
    //                S_cc_ptr[ilev] = &S_cc[ilev];
    //                S_cc[ilev].setVal(0.0); // Set it to zero for this example
    //            }


    //            for (int ilev = 0; ilev < nlevels; ++ilev) {
    //                AmrLevel& ilevel = amr.getLevel(ilev);
    //                field[ilev].define(amr_grids[ilev], amr_dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), ilevel.Factory());
    //                field_ptr[ilev] = &field[ilev];
    //                AmrLevel::FillPatch(ilevel, field[ilev], 1, time, state_idx, B_idx, AMREX_SPACEDIM,0);
    //                const BoxArray & nd_grids = amrex::convert(amr_grids[ilev], IntVect::TheNodeVector()); // nodal grids
    //                rhs[ilev].define(nd_grids, amr_dmap[ilev], 1, 1, MFInfo());
    //                rhs2[ilev].define(grids, amr_dmap[ilev], 1, 1, MFInfo());
    //                rhs[ilev].setVal(0.0);
    //                rhs_ptr1[ilev] = &rhs[ilev];
    //                rhs_ptr2[ilev] = &rhs[ilev];
    //                S_nd[ilev].define(nd_grids, amr_dmap[ilev], 1, 1, MFInfo());
    //                S_nd_ptr[ilev] = &S_nd[ilev];
    //                S_nd[ilev].setVal(0.0); // Set it to zero for this example
    //                MultiFab mf_sigma(amr_grids[ilev], amr_dmap[ilev], 1, 1, MFInfo());
    //                mf_sigma.setVal(1.0);
    //                matrix.setSigma(ilev, mf_sigma);
    //                phi[ilev].define(nd_grids, amr_dmap[ilev], 1, 1, MFInfo());
    //                phi_ptr[ilev] = &phi[ilev];
    //                phi[ilev].setVal(0.0);
    //                phi_cc[ilev].define(grids, amr_dmap[ilev], 1, 1, MFInfo());
    //            }

    //            matrix.compRHS(rhs_ptr1, field_ptr, S_nd_ptr, S_cc_ptr);
    //            MLMG nodal_solver(matrix);
    //            nodal_solver.solve( phi_ptr, rhs_ptr2, 1.e-8, 1.e-15);


    //            //        for (int ilev = 0; ilev < nlevels; ++ilev)
    //            //        {

    //            //            average_node_to_cellcenter(rhs2[ilev], 0, *rhs_ptr2[ilev], 0, 1);
    //            //            average_node_to_cellcenter(phi_cc[ilev], 0, phi[ilev], 0, 1);

    //            //            plot_FAB_2d(rhs2[ilev],0, 1, "rhs (nodal - averaged)", false, false);
    //            //            plot_FAB_2d(phi_cc[ilev],0, 1, "phi (nodal - averaged)", false, false);

    //            //        }


    //            Vector<MultiFab> fluxes(nlevels);
    //            Vector<MultiFab*> fluxes_ptr(nlevels, nullptr);
    //            for (int ilev = 0; ilev < nlevels; ++ilev)
    //            {
    //                fluxes[ilev].define(amr_grids[ilev], amr_dmap[ilev], AMREX_SPACEDIM, 1, MFInfo());
    //                fluxes_ptr[ilev] = &fluxes[ilev];
    //                fluxes[ilev].setVal(0.0);
    //            }

    //            nodal_solver.getFluxes( fluxes_ptr );

    //            for (int ilev = 0; ilev < nlevels; ++ilev)
    //            {
    //                MultiFab::Add(field[ilev], fluxes[ilev], 0, 0, AMREX_SPACEDIM, 0);

    //                average_node_to_cellcenter(rhs2[ilev], 0, *rhs_ptr2[ilev], 0, 1);
    //                average_node_to_cellcenter(phi_cc[ilev], 0, phi[ilev], 0, 1);

    //                //            for (int d=0; d<AMREX_SPACEDIM; ++d) {
    //                //                plot_FAB_2d(fluxes[ilev], d, 1, "grad_phi (nodal) dir="+num2str(d), false, false);
    //                //            }

    //                //            plot_FAB_2d(rhs2[ilev],0, 1, "rhs (nodal - averaged)", false, false);
    //                //            plot_FAB_2d(phi_cc[ilev],0, 1, "phi (nodal - averaged)", false, false);

    //            }


    //        }
    //        //*******************************************

    //        //=====================================================================
    //        // phi (for D field)


    //        if (D_idx > -1) {

    //            // Boundary of the whole domain. This functions must be called,
    //            // and must be called before other bc functions.
    //            const BCRec &bc = istate->boundary_conditions.phys_fill_bc[D_idx];
    //            set_field_bcs(bc_lo, bc_hi, bc, geom);
    //            mlabec.setDomainBC(bc_lo,bc_hi);

    //            // scaling factors; these multiply ACoef and BCoef (see below)
    //            mlabec.setScalars(0.0, 1.0);

    //            Vector<MultiFab*> phi_ptr(nlevels, nullptr);

    //            for (int ilev = 0; ilev < nlevels; ++ilev) {

    //                MultiFab &phi_mf = defined_potential[ilev];
    //                phi_ptr[ilev] = &phi_mf;

    //                // see AMReX_MLLinOp.H for an explanation
    //                mlabec.setLevelBC(ilev, nullptr);

    //                // operator looks like (ACoef - div BCoef grad) phi = rhs
    //                mlabec.setACoeffs(ilev, 0.0);

    //                // think of this beta as the "BCoef" associated with an EB face (typically 1.0)
    //                EB_OPTIONAL(mlabec.setEBDirichlet(ilev, phi_mf, 1.0);)

    //                        // set BCoef (face centered coefficients)
    //                        Array<MultiFab,AMREX_SPACEDIM>& bcf = bcoef[ilev];
    //                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    //                    bcf[idim].setVal(1.0);
    //                }

    //                mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcf));


    //                // the rhs is given by div(D) - cd for the projection method
    //                //defined_charge[ilev] += div(D)

    //                // get filled multifab of the D field
    //                MultiFab D_field(amr_grids[ilev], amr_dmap[ilev], 3, ngrow, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));

    //#ifdef AMREX_USE_EB
    //                EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate->eb2_index));
    //#endif
    //                FillPatch(amr.getLevel(ilev), D_field, ngrow, time, state_idx, D_idx, 3, 0);

    //                const Real* dx = amr_geom[ilev].CellSize();

    //                // calculate the divergence and add it onto the charge density
    //                for (MFIter mfi(defined_rhs[ilev]); mfi.isValid(); ++mfi) {

    //                    const Box& box = mfi.validbox();
    //                    const Box cbox = grow(box,ngrow);

    //#ifdef AMREX_USE_EB
    //                    const EBCellFlagFab& flag = (*getEBData(ilev, state_idx).flags)[mfi];
    //                    if (flag.getType(cbox) == FabType::covered) continue;
    //#endif

    //                    istate->calc_divergence(box, D_field[mfi], defined_rhs[ilev][mfi], EB_OPTIONAL(flag,) dx);

    //                }



    //                // build the rhs pointers
    //                rhs_ptr[ilev] = &defined_rhs[ilev];

    //                //                plot_FAB_2d(defined_potential[ilev],0, "phi", false, false);
    //                //                plot_FAB_2d(defined_rhs[ilev], 0, "rhs", false, false);
    //                //                plot_FAB_2d(D_field, 0, "Dx", false, false);
    //                //                plot_FAB_2d(D_field, 1, "Dy", false, false);
    //                //                plot_FAB_2d(D_field, 2, "Dz", false, true);

    //                //                // set the initial guess for phi as the result of the previous solution
    //                //                AmrLevel& ilevel = amr.getLevel(ilev);
    //                //                MultiFab& ifield = ilevel.get_data(state_idx, time);
    //                //                MultiFab::Copy(phi_mf, ifield, phi_idx, 0, 1, 0);

    //            }

    //            // Solve linear system
    //            mlmg.solve(phi_ptr, rhs_ptr, tol_rel, tol_abs);

    //            // save phi
    //            for (int ilev = 0; ilev < nlevels; ++ilev) {
    //                AmrLevel& ilevel = amr.getLevel(ilev);
    //                MultiFab& ifield = ilevel.get_data(state_idx, time);

    //                MultiFab::Copy(ifield, defined_potential[ilev], 0, phi_idx, 1, 0);

    //                //                plot_FAB_2d(defined_potential[ilev],0, "phi", false, false);
    //            }

    //            mlmg.getFluxes( {GetVecOfArrOfPtrs(bcoef)});  // reuse bcoef


    //            // solve for D field
    //            for (int ilev = 0; ilev < nlevels; ++ilev){
    //                AmrLevel& ilevel = amr.getLevel(ilev);
    //                MultiFab& ifield = ilevel.get_data(state_idx, time);
    //                const Real* dx = amr_geom[ilev].CellSize();

    //                MultiFab grad_phi(amr_grids[ilev], amr_dmap[ilev], 3, 0, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
    //                grad_phi.setVal(0.0);

    //#ifdef AMREX_USE_EB
    //                EB_average_face_to_cellcenter(grad_phi, 0, GetArrOfConstPtrs(bcoef[ilev]));
    //#else
    //                average_face_to_cellcenter(grad_phi, 0, amrex::GetArrOfConstPtrs(bcoef[ilev]));
    //#endif


    //                for (int d = 0; d<AMREX_SPACEDIM; ++d) {
    //                    MultiFab::Subtract(ifield, grad_phi, d, D_idx+d,1,0);
    //                    //                    MultiFab::Saxpy(ifield, -1.0, grad_phi, d, D_idx+d, 1, 0);
    //                }
    //            }
    //        }

    //        //=====================================================================
    //        // psi (for B field)
    //        if (B_idx > -1) {

    //            // Boundary of the whole domain. This functions must be called,
    //            // and must be called before other bc functions.
    //            const BCRec &bc = istate->boundary_conditions.phys_fill_bc[B_idx];
    //            set_field_bcs(bc_lo, bc_hi, bc, geom);
    //            mlabec.setDomainBC(bc_lo,bc_hi);

    //            // scaling factors; these multiply ACoef and BCoef (see below)
    //            mlabec.setScalars(0.0, 1.0);

    //            Vector<MultiFab*> psi_ptr(nlevels, nullptr);

    //            for (int ilev = 0; ilev < nlevels; ++ilev) {

    //                MultiFab &psi_mf = defined_potential[ilev];
    //                psi_ptr[ilev] = &psi_mf;
    //                psi_mf.setVal(0.0);

    //                // see AMReX_MLLinOp.H for an explanation
    //                mlabec.setLevelBC(ilev, nullptr);

    //                // operator looks like (ACoef - div BCoef grad) phi = rhs
    //                mlabec.setACoeffs(ilev, 0.0);

    //                // think of this beta as the "BCoef" associated with an EB face (typically 1.0)
    //                EB_OPTIONAL(mlabec.setEBDirichlet(ilev, psi_mf, 1.0);)

    //                        // set BCoef (face centered coefficients)
    //                        Array<MultiFab,AMREX_SPACEDIM>& bcf = bcoef[ilev];
    //                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    //                    bcf[idim].setVal(1.0);
    //                }

    //                mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcf));

    //                // the rhs is given by div(B)  for the projection method

    //                // get filled multifab of the B field
    //                MultiFab B_field(amr_grids[ilev], amr_dmap[ilev], 3, ngrow, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));

    //#ifdef AMREX_USE_EB
    //                EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate->eb2_index));
    //#endif
    //                FillPatch(amr.getLevel(ilev), B_field, ngrow, time, state_idx, B_idx, 3, 0);

    //                const Real* dx = amr_geom[ilev].CellSize();

    //                defined_rhs[ilev].setVal(0.0);

    //                // calculate the divergence
    //                for (MFIter mfi(defined_rhs[ilev]); mfi.isValid(); ++mfi) {

    //                    const Box& box = mfi.validbox();
    //                    const Box cbox = grow(box,ngrow);

    //#ifdef AMREX_USE_EB
    //                    const EBCellFlagFab& flag = (*getEBData(ilev, state_idx).flags)[mfi];
    //                    if (flag.getType(cbox) == FabType::covered) continue;
    //#endif

    //                    istate->calc_divergence(box, B_field[mfi], defined_rhs[ilev][mfi], EB_OPTIONAL(flag,) dx);

    //                }

    //                // build the rhs pointers
    //                rhs_ptr[ilev] = &defined_rhs[ilev];

    //            }

    //            // Solve linear system
    //            mlmg.solve(psi_ptr, rhs_ptr, tol_rel, tol_abs);

    //            //            for (int ilev = 0; ilev < nlevels; ++ilev) {
    //            //                plot_FAB_2d(*rhs_ptr[ilev], 0, 1, "rhs (cell)", false, false);
    //            //                plot_FAB_2d(*psi_ptr[ilev], 0, 1, "phi (cell)", false, true);
    //            //            }

    //            // save psi
    //            for (int ilev = 0; ilev < nlevels; ++ilev) {
    //                AmrLevel& ilevel = amr.getLevel(ilev);
    //                MultiFab& ifield = ilevel.get_data(state_idx, time);

    //                MultiFab::Copy(ifield, *psi_ptr[ilev], 0, psi_idx, 1, 0);
    //            }



    //            // Get fluxes from solver
    //            mlmg.getFluxes( {GetVecOfArrOfPtrs(bcoef)});  // reuse bcoef

    //            // solve for B field
    //            for (int ilev = 0; ilev < nlevels; ++ilev){
    //                AmrLevel& ilevel = amr.getLevel(ilev);
    //                MultiFab& ifield = ilevel.get_data(state_idx, time);

    //                MultiFab grad_psi(amr_grids[ilev], amr_dmap[ilev], 3, 0, MFInfo() EB_OPTIONAL(,*amr_factory[ilev]));
    //                grad_psi.setVal(0.0);

    //#ifdef AMREX_USE_EB
    //                EB_average_face_to_cellcenter(grad_psi, 0, GetArrOfConstPtrs(bcoef[ilev]));
    //#else
    //                average_face_to_cellcenter(grad_psi, 0, GetArrOfConstPtrs(bcoef[ilev]));
    //#endif
    //                const Real* dx = amr_geom[ilev].CellSize();
    //                for (int d = 0; d<AMREX_SPACEDIM; ++d) {
    //                    //                    MultiFab::Saxpy(ifield, -1.0, grad_psi, d, B_idx+d, 1, 0);
    //                    MultiFab::Subtract(ifield, grad_psi, d, B_idx+d,1,0);

    //                    //                    plot_FAB_2d(grad_psi, d, 1, "grad_phi (cell) dir="+num2str(d), false, false);
    //                }
    //            }
    //        }
    //    }
    //#endif
}



void MFP::project_divergence(const Real time)
{
    BL_PROFILE("MFP::project_divergence");
    // MLNodeLaplacian not implemented for 1D
#if AMREX_SPACEDIM > 1
    for (auto &istate : gd.states) {

        if (!((istate->has_field()) && (istate->is_transported())))
            continue;

        if (!istate->project_divergence)
            continue;

        //----
        // B field divergence clean
        if (istate->get_type() == +StateType::isMHD || istate->get_type() == +StateType::isField) {

            const int B_idx = istate->get_cons_B_idx();
            const int psi_idx = istate->get_cons_psi_idx();

            solve_divergence(istate->global_idx, B_idx, psi_idx, istate->boundary_conditions.phys_fill_bc[B_idx]);
        }


        //----
        // D field divergence clean
        if (istate->get_type() == +StateType::isField) {
            Real cd_scale = gd.Larmor/(gd.Debye*gd.Debye*gd.lightspeed);

            // assume zero charge density
            Amr& amr = *parent;
            //            const Vector<Geometry>& amr_geom = amr.Geom();
            const Vector<BoxArray>& amr_grids = amr.boxArray(0, amr.finestLevel());
            const Vector<DistributionMapping>& amr_dmap = amr.DistributionMap(0, amr.finestLevel());
            const int nlevels = amr_grids.size();

            Vector<MultiFab> charge_density(nlevels);

            const int ngrow = 1;

            // go through level-by-level and calculate the charge and current density
            for (int ilev = 0; ilev < nlevels; ++ilev) {

#ifdef AMREX_USE_EB
                EBFArrayBoxFactory const* amr_factory = getEBData(ilev,istate->global_idx).ebfactory.get();
#endif

                charge_density[ilev].define(amr_grids[ilev], amr_dmap[ilev], 1, ngrow, MFInfo() EB_OPTIONAL(,*amr_factory));
                charge_density[ilev].setVal(0.0);

                AmrLevel& ilevel = amr.getLevel(ilev);
                EB_OPTIONAL(Vector<EBData>& eb_data = gd.eb_data[ilev];)

                        // get all of the source contributions
                for (const auto& src_idx : istate->associated_sources) {
                    const auto &ode = gd.ode_source_terms[src_idx.first];
                    const auto &src = ode->sources[src_idx.second];

                    if (!src->has_charge_density()) continue;

                    // now calculate the charge and current density

                    Vector<FArrayBox*> src_data(src->offsets.size());
                    EB_OPTIONAL(Vector<const EBCellFlagFab*> src_flags(src->offsets.size());)
                    FArrayBox local_cd;
                    FArrayBox local_J; // dummy variable needed for the calc_charge_density call


                    // iterate over data
                    for (MFIter mfi(ilevel.get_data(0, time)); mfi.isValid(); ++mfi) {
                        const Box& box = mfi.tilebox();
                        local_cd.resize(box);
                        local_J.resize(box,3);

                        // get lists of the data for each state included in the src
                        for (const auto& idx : src->offsets) {
                            src_data[idx.local] = &ilevel.get_data(idx.global, time)[mfi];
                            EB_OPTIONAL(src_flags[idx.local] = &(eb_data[idx.global].flags)[mfi];)
                        }

                        src->calc_charge_density(box,
                                                 geom.ProbLo(),
                                                 geom.CellSize(),
                                                 parent->cumTime(),
                                                 src_data,
                                                 local_cd,
                                                 local_J
                                                 EB_OPTIONAL(,src_flags)
                                                 );

                        // add the charge density to the accumulator
                        charge_density[ilev][mfi].plus(local_cd);
                    }
                }

                // scale the charge and current density appropriately
                charge_density[ilev].mult(-cd_scale, 0);
            }

            //                plot_FABs_2d(charge_density[0],0, "charge density", false, true);

            ParallelDescriptor::Barrier();

            // now that we have the charge density over all levels do the actual solve
            const int D_idx = istate->get_cons_D_idx();
            const int phi_idx = istate->get_cons_phi_idx();
            solve_divergence(istate->global_idx, D_idx, phi_idx, istate->boundary_conditions.phys_fill_bc[D_idx], GetVecOfPtrs(charge_density));
        }
    }
#endif
}

void MFP::solve_divergence(const int state_idx,
                           const int vector_idx,
                           const int phi_idx,
                           const BCRec &bc,
                           Vector<MultiFab*> S_cc_ptr)
{
    BL_PROFILE("MFP::solve_divergence");
    Amr& amr = *parent;

    //
    // Given a cell-centered velocity (vel) field, a cell-centered
    // scalar field (sigma) field, and a source term S (either node-
    // or cell-centered) solve:
    //
    //   div( sigma * grad(phi) ) = div(vel) - S
    //
    // and then perform the projection:
    //
    //     vel = vel - sigma * grad(phi)
    //
    // NOTE: This is only approximate on the order of cell spacing i.e. O(dx^2)

    const Vector<Geometry>& amr_geom = amr.Geom(0, amr.finestLevel());
    const Vector<BoxArray>& amr_grids = amr.boxArray(0, amr.finestLevel());
    const Vector<DistributionMapping>& amr_dmap = amr.DistributionMap(0, amr.finestLevel());

    const int nlevels = amr_geom.size();

#ifdef AMREX_USE_EB
    Vector<EBFArrayBoxFactory const*> amr_factory(nlevels);

    for (int ilev=0; ilev<nlevels; ++ilev) {
        amr_factory[ilev] = getEBData(ilev,state_idx).ebfactory.get();
    }
#endif


    // Linear Solver
    LPInfo lp_info;

    MLNodeLaplacian matrix(amr_geom,
                           amr_grids,
                           amr_dmap,
                           lp_info
                           EB_OPTIONAL(,amr_factory)
                           );

    // Set boundary conditions.
    // Note that Dirichlet boundary conditions are assumed to be homogeneous (i.e. phi = 0)

    // we assume that the BCRec given for the first component of the field vector
    // is representative and thus use it for determining the overall BC for phi

    Array<LinOpBCType, AMREX_SPACEDIM> lobc;
    Array<LinOpBCType, AMREX_SPACEDIM> hibc;

    set_field_bcs(lobc,hibc,bc,geom);

    matrix.setDomainBC(lobc,hibc);

    // Set matrix attributes to be used by MLMG solver
    matrix.setGaussSeidel(false);
    matrix.setHarmonicAverage(false);

    Vector<MultiFab> field(nlevels);
    Vector<MultiFab> rhs(nlevels);
    Vector<MultiFab> S_nd(nlevels);
    Vector<MultiFab> S_cc;
    Vector<MultiFab> phi(nlevels);

    Vector<MultiFab*> field_ptr(nlevels, nullptr);
    Vector<MultiFab*> rhs_ptr(nlevels, nullptr);
    Vector<const MultiFab*> rhs_ptr2(nlevels, nullptr);
    Vector<const MultiFab*> S_nd_ptr(nlevels, nullptr);
    Vector<MultiFab*> phi_ptr(nlevels, nullptr);


    // Cell-centered contributions to RHS
    // manufacture zero values if necessary
    if (S_cc_ptr.empty()) {
        S_cc.resize(nlevels);
        S_cc_ptr.resize(nlevels,nullptr);
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            S_cc[ilev].define(amr_grids[ilev], amr_dmap[ilev], 1, 1, MFInfo());
            S_cc_ptr[ilev] = &S_cc[ilev];
            S_cc[ilev].setVal(0.0); // Set it to zero for this example
        }
    }

    const Real time = amr.cumTime();


    for (int ilev = 0; ilev < nlevels; ++ilev) {

        AmrLevel& ilevel = amr.getLevel(ilev);

        //
        //  Create the cell-centered field we want to project
        //
        field[ilev].define(amr_grids[ilev], amr_dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), ilevel.Factory());
        field_ptr[ilev] = &field[ilev];

        // Set field including ghost cells

#ifdef AMREX_USE_EB
        State& istate = GlobalData::get_state(state_idx);
        EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate.eb2_index));
#endif

        AmrLevel::FillPatch(ilevel, field[ilev], 1, time, state_idx, vector_idx, AMREX_SPACEDIM,0);

        // RHS is nodal
        const BoxArray & nd_grids = amrex::convert(amr_grids[ilev], IntVect::TheNodeVector()); // nodal grids

        // Multifab to host RHS
        rhs[ilev].define(nd_grids, amr_dmap[ilev], 1, 1, MFInfo());
        rhs[ilev].setVal(0.0);
        rhs_ptr[ilev] = &rhs[ilev];
        rhs_ptr2[ilev] = &rhs[ilev];

        // Node-centered contributions to RHS
        S_nd[ilev].define(nd_grids, amr_dmap[ilev], 1, 1, MFInfo());
        S_nd_ptr[ilev] = &S_nd[ilev];
        S_nd[ilev].setVal(0.0); // Set it to zero for this example

        //
        // Create the cell-centered sigma field and set it to 1.0
        //
        MultiFab mf_sigma(amr_grids[ilev], amr_dmap[ilev], 1, 1, MFInfo());
        mf_sigma.setVal(1.0);

        // Set sigma
        matrix.setSigma(ilev, mf_sigma);

        //
        // Create node-centered phi
        //
        phi[ilev].define(nd_grids, amr_dmap[ilev], 1, 1, MFInfo());
        phi_ptr[ilev] = &phi[ilev];
        phi[ilev].setVal(0.0);

        //
        // Create the cell-centered boundary condition data
        //
        //        MultiFab levelbcdata(amr_grids[ilev], amr_dmap[ilev], 1, 1, MFInfo());
        //        levelbcdata.setVal(0.0);

        //        matrix.setLevelBC(ilev, &levelbcdata);
    }

    // Compute RHS -- field must be cell-centered
    matrix.compRHS(rhs_ptr, field_ptr, S_nd_ptr, S_cc_ptr);


    //
    // Setup MLMG solver
    //
    MLMG nodal_solver(matrix);

    // We can specify the maximum number of iterations
    // nodal_solver.setMaxIter(nodal_mg_maxiter);
    // nodal_solver.setCGMaxIter(nodal_mg_cg_maxiter);

    int v = parent->Verbose();
    nodal_solver.setVerbose(v);

    // Define the relative tolerance
    Real reltol = 1.e-8;

    // Define the absolute tolerance; note that this argument is optional
    Real abstol = 1.e-15;

    //
    // Solve div( sigma * grad(phi) ) = RHS
    //
    nodal_solver.solve( phi_ptr, rhs_ptr2, reltol, abstol);


    //
    // Create cell-centered multifab to hold value of -sigma*grad(phi) at cell-centers
    //
    Vector<MultiFab> fluxes(nlevels);
    Vector<MultiFab*> fluxes_ptr(nlevels, nullptr);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        fluxes[ilev].define(amr_grids[ilev], amr_dmap[ilev], AMREX_SPACEDIM, 1, MFInfo());
        fluxes_ptr[ilev] = &fluxes[ilev];
        fluxes[ilev].setVal(0.0);
    }

    // Get fluxes from solver
    nodal_solver.getFluxes( fluxes_ptr );

    //
    // Apply projection explicitly --  vel = vel - sigma * grad(phi)
    //
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        AmrLevel& ilevel = amr.getLevel(ilev);
        MultiFab& ifield = ilevel.get_data(state_idx, time);
        MultiFab::Add(ifield, fluxes[ilev], 0, vector_idx, AMREX_SPACEDIM, 0);

        average_node_to_cellcenter(ifield, phi_idx, phi[ilev], 0, 1);

    }

    return;
}
