#include <AMReX_Array.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_Vector.H>

#ifdef AMREX_USE_EB
    #include <AMReX_EB2.H>
    #include <AMReX_MLEBABecLap.H>
#else
    #include <AMReX_MLABecLaplacian.H>
#endif

using Location = amrex::MLLinOp::Location;

#include "MFP.H"
#include "MFP_diagnostics.H"
#include "MFP_elliptic.H"
#include "MFP_lorentz.H"
#include "MFP_plasma5.H"
#include "MFP_state.H"
#include "sol.hpp"

std::string Elliptic::tag = "elliptic";
bool Elliptic::registered = GetActionFactory().Register(Elliptic::tag, ActionBuilder<Elliptic>);

Elliptic::Elliptic() {}
Elliptic::~Elliptic() {}

Elliptic::Elliptic(const int idx, const sol::table& def)
{
    BL_PROFILE("Elliptic::Elliptic");

    action_idx = idx;
    name = def["name"];

    projection = def.get_or("projection", 0);

    // Define the relative tolerance
    reltol = def.get_or("relative_tolerance", 1.e-8);

    // Define the absolute tolerance; note that this argument is optional
    abstol = def.get_or("absolute_tolerance", 1.e-15);

    std::string field_state_name = def.get_or<std::string>("state", "null");

    if (field_state_name == "null")
        Abort("Action '" + name + "' requires option 'state' to be defined.");

    State& istate = MFP::get_state(field_state_name);

    select = istate.get_type();

    switch (select) {
    case State::StateType::Field:
        field = static_cast<FieldState*>(&istate);
        state_indexes.push_back(istate.global_idx);
        break;
    case State::StateType::MHD:
        mhd = static_cast<MHDState*>(&istate);
        state_indexes.push_back(istate.global_idx);
        break;
    default: Abort("An invalid state has been defined for the '" + name + "' action");
    }

    return;
}

void Elliptic::apply_change(MFP* mfp, const Real time, const Real dt)
{
    BL_PROFILE("Elliptic::apply_change");

    if (mfp->get_level() > 0) return;  // only do this once per coarse timestep

    if (projection) {
        project_divergence(mfp, time);
    } else {
        solve_static_fields(mfp, time);
    }
}

void set_field_bcs(Array<LinOpBCType, AMREX_SPACEDIM>& bc_lo,
                   Array<LinOpBCType, AMREX_SPACEDIM>& bc_hi,
                   const BCRec& bc,
                   const Geometry& geom)
{
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (geom.isPeriodic(d)) {
            bc_lo[d] = LinOpBCType::Periodic;
            bc_hi[d] = LinOpBCType::Periodic;
        } else {
            switch (bc.lo(d)) {
            case PhysBCType::interior: bc_lo[d] = LinOpBCType::Periodic; break;

            case PhysBCType::inflow: bc_lo[d] = LinOpBCType::Dirichlet; break;

            case PhysBCType::outflow: bc_lo[d] = LinOpBCType::Neumann; break;

            case PhysBCType::symmetry: bc_lo[d] = LinOpBCType::Neumann; break;

            case PhysBCType::slipwall: bc_lo[d] = LinOpBCType::reflect_odd; break;

            case PhysBCType::noslipwall: bc_lo[d] = LinOpBCType::reflect_odd; break;
            }

            switch (bc.hi(d)) {
            case PhysBCType::interior: bc_hi[d] = LinOpBCType::Periodic; break;

            case PhysBCType::inflow: bc_hi[d] = LinOpBCType::Dirichlet; break;

            case PhysBCType::outflow: bc_hi[d] = LinOpBCType::Neumann; break;

            case PhysBCType::symmetry: bc_hi[d] = LinOpBCType::Neumann; break;

            case PhysBCType::slipwall: bc_hi[d] = LinOpBCType::reflect_odd; break;

            case PhysBCType::noslipwall: bc_hi[d] = LinOpBCType::reflect_odd; break;
            }
        }
    }
}

void Elliptic::solve_static_fields(MFP* mfp, const Real time)
{
    BL_PROFILE("Elliptic::solve_static_fields");

#if AMREX_SPACEDIM > 1

    const int state_idx = field->data_idx;

    Real cd_scale = (MFP::Larmor) / (MFP::Debye * MFP::Debye * MFP::lightspeed);
    Real J_scale = (MFP::Larmor) / (MFP::Debye * MFP::Debye * MFP::lightspeed * MFP::lightspeed);

    // assume zero charge density
    Amr& amr = *mfp->get_parent();

    const Vector<BoxArray>& amr_grids = amr.boxArray(0, amr.finestLevel());
    const Vector<DistributionMapping>& amr_dmap = amr.DistributionMap(0, amr.finestLevel());
    const int nlevels = amr_grids.size();

    const Vector<Geometry>& amr_geom = amr.Geom(0, amr.finestLevel());
    const Geometry geom = mfp->Geom();

    #ifdef AMREX_USE_EB
    Vector<EBFArrayBoxFactory const*> amr_factory(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev) {
        amr_factory[ilev] = mfp->get_eb_data(ilev, state_idx).ebfactory.get();
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
        defined_phi[ilev].define(amr_grids[ilev],
                                 amr_dmap[ilev],
                                 1,
                                 ngrow,
                                 MFInfo()
    #ifdef AMREX_USE_EB
                                   ,
                                 *amr_factory[ilev]
    #endif
        );
        defined_phi[ilev].setVal(0.0);

        defined_A[ilev].define(amr_grids[ilev],
                               amr_dmap[ilev],
                               3,
                               ngrow,
                               MFInfo()
    #ifdef AMREX_USE_EB
                                 ,
                               *amr_factory[ilev]
    #endif
        );
        defined_A[ilev].setVal(0.0);

        defined_charge[ilev].define(amr_grids[ilev],
                                    amr_dmap[ilev],
                                    1,
                                    ngrow,
                                    MFInfo()
    #ifdef AMREX_USE_EB
                                      ,
                                    *amr_factory[ilev]
    #endif
        );
        defined_charge[ilev].setVal(0.0);

        defined_current[ilev].define(amr_grids[ilev],
                                     amr_dmap[ilev],
                                     3,
                                     ngrow,
                                     MFInfo()
    #ifdef AMREX_USE_EB
                                       ,
                                     *amr_factory[ilev]
    #endif
        );
        defined_current[ilev].setVal(0.0);

        AmrLevel& ilevel = amr.getLevel(ilev);

    #ifdef AMREX_USE_EB
        Vector<EBData>& eb_data = mfp->eb_data[ilev];
    #endif

        // get all of the source contributions
        for (const auto& src_idx : field->associated_actions) {
            Vector<HydroState*> hydro_states;

            switch (mfp->actions[src_idx]->get_type()) {
            case ActionType::Plasma5: {
                Plasma5& plasma = static_cast<Plasma5&>(*(mfp->actions[src_idx]));
                hydro_states = plasma.species;
                break;
            }
            case ActionType::Lorentz: {
                Lorentz& plasma = static_cast<Lorentz&>(*(mfp->actions[src_idx]));
                hydro_states = plasma.species;
                break;
            }
            default: continue;
            }

            // now calculate the charge and current density
            const size_t num_src = hydro_states.size();

            // iterate over data
            for (MFIter mfi(ilevel.get_new_data(MFP::Cost_Idx)); mfi.isValid(); ++mfi) {
                const Box& box = mfi.tilebox();

                FArrayBox& local_cd = defined_charge[ilev][mfi];
                local_cd.setVal(0.0);
                FArrayBox& local_J = defined_current[ilev][mfi];
                local_J.setVal(0.0);

                // get lists of the data for each state included in the src
                for (size_t src_idx = 0; src_idx < num_src; ++src_idx) {
                    HydroState& hydro = *hydro_states[src_idx];
                    const FArrayBox& cons = ilevel.get_new_data(hydro.data_idx)[mfi];
    #ifdef AMREX_USE_EB
                    const FArrayBox& vfrac = eb_data[hydro.global_idx].volfrac[mfi];
    #endif

                    hydro.calc_current_and_charge(box,
                                                  cons,
                                                  &local_cd,
                                                  &local_J

    #ifdef AMREX_USE_EB
                                                  ,
                                                  vfrac
    #endif
                    );
                }

                // scale by the relative permittivity and permeability

                FArrayBox& field_data = ilevel.get_new_data(field->data_idx)[mfi];

                local_cd.divide(field_data, box, +FieldDef::ConsIdx::ep, 0, 1);
                local_J.divide(field_data, box, +FieldDef::ConsIdx::mu, 0, 1);
                local_J.divide(field_data, box, +FieldDef::ConsIdx::mu, 1, 1);
                local_J.divide(field_data, box, +FieldDef::ConsIdx::mu, 2, 1);
            }
        }

        // get the embedded boundary contributions
    #ifdef AMREX_USE_EB
        Vector<FArrayBox*> bcs_data(+FieldBoundaryEB::EBType::NUM);

        for (MFIter mfi(defined_phi[ilev]); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();
            const Box cbox = grow(box, ngrow);

            // specify here what EB bcs we are looking to get info from
            bcs_data[+FieldBoundaryEB::EBType::ScalarPotential] = &defined_phi[ilev][mfi];
            bcs_data[+FieldBoundaryEB::EBType::VectorPotential] = &defined_A[ilev][mfi];

            const EBData& ebd = mfp->get_eb_data(ilev, state_idx);
            const EBCellFlagFab& flag = ebd.flags[mfi];

            if (flag.getType(cbox) != FabType::singlevalued) continue;

            const FArrayBox& bnorm = (*ebd.bndrynorm)[mfi];
            const FArrayBox& bcent = (*ebd.bndrycent)[mfi];

            const CutFab& bc_idx = ebd.bndryidx[mfi];

            //            plot_FAB_2d(bc_idx,0,"bc idx", false, true);

            field->get_wall_value(cbox,
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
        if (defined_phi[ilev].norminf(0, ngrow, false, true) > 0.0) { update_D = true; }

        if (!update_D) {
            if (defined_charge[ilev].norminf(0, ngrow, false, true) > 0.0) { update_D = true; }
        }

        // check if we need to update B
        for (int d = 0; d < 3; ++d) {
            if (defined_A[ilev].norminf(d, ngrow, false, true) > 0.0) {
                update_B = true;
                break;
            }

            if (defined_current[ilev].norminf(d, ngrow, false, true) > 0.0) {
                update_B = true;
                break;
            }
        }

        //        plot_FAB_2d(defined_phi[ilev],0,0, "phi", false, true);
        //        plot_FAB_2d(defined_A[ilev],0,0, "Ax", false, false);
        //        plot_FAB_2d(defined_A[ilev],1,0, "Ay", false, false);
        //        plot_FAB_2d(defined_A[ilev],2,0, "Az", false, true);

        //        plot_FAB_2d(defined_charge[ilev],0,0,"cd", false, false);
        //        plot_FAB_2d(defined_current[ilev],0,0,"Jx", false, false);
        //        plot_FAB_2d(defined_current[ilev],1,0,"Jy", false, false);
        //        plot_FAB_2d(defined_current[ilev],2,0,"Jz", false, true);
    }

    if (!update_D && !update_B) return;

    ParallelDescriptor::Barrier();

    // now that we have the charge density over all levels do the actual solve

    //---------------------------------------------------------------------------------------

    // Linear Solver
    LPInfo lp_info;

    #ifdef AMREX_USE_EB
    MLEBABecLap mlabec(amr_geom, amr_grids, amr_dmap, lp_info, amr_factory);
    #else
    MLABecLaplacian mlabec(amr_geom, amr_grids, amr_dmap, lp_info);
    #endif

    MLMG mlmg(mlabec);

    mlmg.setVerbose(MFP::linear_solver_verbosity);

    // define array of LinOpBCType for domain boundary conditions
    Array<LinOpBCType, AMREX_SPACEDIM> bc_lo;
    Array<LinOpBCType, AMREX_SPACEDIM> bc_hi;

    Vector<const MultiFab*> rhs_ptr(nlevels, nullptr);
    Vector<Array<MultiFab, AMREX_SPACEDIM>> bcoef(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        // set BCoef (face centered coefficients)
        Array<MultiFab, AMREX_SPACEDIM>& bcf = bcoef[ilev];
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcf[idim].define(amrex::convert(amr_grids[ilev], IntVect::TheDimensionVector(idim)),
                             amr_dmap[ilev],
                             1,
                             0,
                             MFInfo()
    #ifdef AMREX_USE_EB
                               ,
                             *amr_factory[ilev]
    #endif
            );
        }
    }

    //=====================================================================
    // phi (for D field)

    if (update_D) {
        // Boundary of the whole domain. This functions must be called,
        // and must be called before other bc functions.
        const BCRec& bc = field->boundary_conditions.phys_fill_bc[+FieldDef::ConsIdx::Dx];
        set_field_bcs(bc_lo, bc_hi, bc, geom);
        mlabec.setDomainBC(bc_lo, bc_hi);

        // scaling factors; these multiply ACoef and BCoef (see below)
        mlabec.setScalars(0.0, 1.0);

        Vector<MultiFab*> phi_ptr(nlevels, nullptr);

        for (int ilev = 0; ilev < nlevels; ++ilev) {
            MultiFab& phi_mf = defined_phi[ilev];
            phi_ptr[ilev] = &phi_mf;

            // see AMReX_MLLinOp.H for an explanation
            mlabec.setLevelBC(ilev, nullptr);

            // operator looks like (ACoef - div BCoef grad) phi = rhs
            mlabec.setACoeffs(ilev, 0.0);

            // think of this beta as the "BCoef" associated with an EB face (typically 1.0)
    #ifdef AMREX_USE_EB
            mlabec.setEBDirichlet(ilev, phi_mf, 1.0);
    #endif

            // set BCoef (face centered coefficients)
            Array<MultiFab, AMREX_SPACEDIM>& bcf = bcoef[ilev];
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) { bcf[idim].setVal(1.0); }

            mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcf));

            // build the rhs pointers
            rhs_ptr[ilev] = &defined_charge[ilev];
        }

        //        plot_FAB_2d(defined_phi[nlevels-1],0,1, "phi", false, false);
        //        plot_FAB_2d(defined_charge[nlevels-1],0,1, "cd", false, true);

        // Solve linear system
        mlmg.solve(phi_ptr, rhs_ptr, reltol, abstol);

        // save phi
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            AmrLevel& ilevel = amr.getLevel(ilev);
            MultiFab& ifield = ilevel.get_new_data(state_idx);

            MultiFab::Copy(ifield, defined_phi[ilev], 0, +FieldDef::ConsIdx::phi, 1, 0);

            //            plot_FAB_2d(ifield,+FieldDef::ConsIdx::phi,ifield.nGrow(), "phi", false,
            //            true);
        }

        // Get fluxes from solver
        mlmg.getFluxes({GetVecOfArrOfPtrs(bcoef)});  // reuse bcoef

        // solve for D field
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            AmrLevel& ilevel = amr.getLevel(ilev);
            MultiFab& ifield = ilevel.get_new_data(state_idx);

    #ifdef AMREX_USE_EB
            EB_average_face_to_cellcenter(ifield,
                                          +FieldDef::ConsIdx::Dx,
                                          GetArrOfConstPtrs(bcoef[ilev]));
    #else
            average_face_to_cellcenter(ifield,
                                       +FieldDef::ConsIdx::Dx,
                                       amrex::GetArrOfConstPtrs(bcoef[ilev]));
    #endif
        }
    }

    //=====================================================================
    // A (for B field)
    if (update_B) {
        // Boundary of the whole domain. This functions must be called,
        // and must be called before other bc functions.
        const BCRec& bc = field->boundary_conditions.phys_fill_bc[+FieldDef::ConsIdx::Bx];
        set_field_bcs(bc_lo, bc_hi, bc, geom);
        mlabec.setDomainBC(bc_lo, bc_hi);

        // scaling factors; these multiply ACoef and BCoef (see below)
        mlabec.setScalars(0.0, 1.0);

        Vector<Array<MultiFab, 3>> A(nlevels);
        Vector<MultiFab*> A_ptr(nlevels, nullptr);

        Vector<MultiFab> current_alias(nlevels);

        // zero out the B field value
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            AmrLevel& ilevel = amr.getLevel(ilev);
            MultiFab& ifield = ilevel.get_new_data(state_idx);
            ifield.setVal(0.0, +FieldDef::ConsIdx::Bx, 3);

            mlabec.setLevelBC(ilev, nullptr);
            mlabec.setACoeffs(ilev, 0.0);
        }

        for (int d = 0; d < 3; ++d) {
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                MultiFab& A_mf = A[ilev][d];
                A_mf = MultiFab(defined_A[ilev], MakeType::make_alias, d, 1);
                A_ptr[ilev] = &A_mf;

    #ifdef AMREX_USE_EB
                mlabec.setEBDirichlet(ilev, A_mf, 1.0);
    #endif

                // set BCoef (face centered coefficients)
                Array<MultiFab, AMREX_SPACEDIM>& bcf = bcoef[ilev];
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) { bcf[idim].setVal(1.0); }

                mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcf));

                // build the rhs pointers

                current_alias[ilev] = MultiFab(defined_current[ilev], MakeType::make_alias, d, 1);

                rhs_ptr[ilev] = &current_alias[ilev];
            }

            //            plot_FAB_2d(*(A_ptr[nlevels-1]),0,0, "A"+num2str(d), false, false);
            //            plot_FAB_2d(*(rhs_ptr[nlevels-1]),0,0, "J"+num2str(d), false, true);

            // Solve linear system
            mlmg.solve(A_ptr, rhs_ptr, reltol, abstol);

            //                plot_FAB_2d(*(A_ptr[nlevels-1]),0,0, "A"+num2str(d), false, true);

            // Get fluxes from solver
            mlmg.getFluxes({GetVecOfArrOfPtrs(bcoef)});  // reuse bcoef

            // save A for debugging
            if (false) {
                for (int ilev = 0; ilev < nlevels; ++ilev) {
                    AmrLevel& ilevel = amr.getLevel(ilev);
                    MultiFab& ifield = ilevel.get_new_data(state_idx);

                    MultiFab::Copy(ifield,
                                   A[ilev][d],
                                   0,
                                   +FieldDef::ConsIdx::Bx + d,
                                   1,
                                   0);  // magnetic vector potential

                    //                        MultiFab::Copy(ifield, defined_current[ilev], d,
                    //                        B_idx+d, 1, 0); // current
                }
            } else {
                // solve for B field
                for (int ilev = 0; ilev < nlevels; ++ilev) {
                    AmrLevel& ilevel = amr.getLevel(ilev);
                    MultiFab& ifield = ilevel.get_new_data(state_idx);

                    // calculate the gradients
                    MultiFab gradients(amr_grids[ilev],
                                       amr_dmap[ilev],
                                       3,
                                       0,
                                       MFInfo()
    #ifdef AMREX_USE_EB
                                         ,
                                       *amr_factory[ilev]
    #endif
                    );

    #ifdef AMREX_USE_EB
                    EB_average_face_to_cellcenter(gradients, 0, GetArrOfConstPtrs(bcoef[ilev]));
    #else
                    average_face_to_cellcenter(gradients, 0, amrex::GetArrOfConstPtrs(bcoef[ilev]));
    #endif

                    // we now have -sigma*grad(A) at the cell centre (sigma also known as B in the
                    // docs) iteratively calculate the curl and thus B, i.e. B = curl A
                    switch (d) {
                    case 0:
                        MultiFab::Add(ifield,
                                      gradients,
                                      1,
                                      +FieldDef::ConsIdx::Bz,
                                      1,
                                      0);  // dAx/dy -> Bz
    #if AMREX_SPACEDIM == 3
                        MultiFab::Subtract(ifield,
                                           gradients,
                                           2,
                                           +FieldDef::ConsIdx::By,
                                           1,
                                           0);  // dAx/dz -> By
    #endif
                        break;
                    case 1:
                        MultiFab::Subtract(ifield,
                                           gradients,
                                           0,
                                           +FieldDef::ConsIdx::Bz,
                                           1,
                                           0);  // dAy/dx -> Bz
    #if AMREX_SPACEDIM == 3
                        MultiFab::Add(ifield,
                                      gradients,
                                      2,
                                      +FieldDef::ConsIdx::By,
                                      1,
                                      0);  // dAy/dz -> Bx
    #endif
                        break;
                    case 2:
                        MultiFab::Add(ifield,
                                      gradients,
                                      0,
                                      +FieldDef::ConsIdx::By,
                                      1,
                                      0);  // dAz/dx -> By
                        MultiFab::Subtract(ifield,
                                           gradients,
                                           1,
                                           +FieldDef::ConsIdx::Bx,
                                           1,
                                           0);  // dAz/dy -> Bx

                        break;
                    default: Abort("how did we get here?");
                    }
                }
            }
        }

        //            plot_FAB_2d(current_alias[nlevels-1],0, "Jz", false, true);
    }

#endif
}

void Elliptic::project_divergence(MFP* mfp, const Real time)
{
    BL_PROFILE("Elliptic::project_divergence");

    // MLNodeLaplacian not implemented for 1D
#if AMREX_SPACEDIM > 1
    //----
    // B field divergence clean

    switch (select) {
    case State::StateType::Field: project_divergence_field(mfp, time); break;
    case State::StateType::MHD: project_divergence_mhd(mfp, time); break;
    default: Abort("An invalid state has been defined for the '" + name + "' action");
    }
#endif
}

void Elliptic::project_divergence_mhd(MFP* mfp, const Real time)
{
    BL_PROFILE("Elliptic::project_divergence_mhd");

    // MLNodeLaplacian not implemented for 1D
#if AMREX_SPACEDIM > 1
    solve_divergence(mfp,
                     mhd->data_idx,
                     +MHDDef::ConsIdx::Bx,
                     +MHDDef::ConsIdx::psi,
                     mhd->boundary_conditions.phys_fill_bc[+MHDDef::ConsIdx::Bx]);
#endif
}

void Elliptic::project_divergence_field(MFP* mfp, const Real time)
{
    BL_PROFILE("Elliptic::project_divergence_field");

    // MLNodeLaplacian not implemented for 1D
#if AMREX_SPACEDIM > 1

    //----
    // B field divergence clean

    solve_divergence(mfp,
                     field->data_idx,
                     +FieldDef::ConsIdx::Bx,
                     +FieldDef::ConsIdx::psi,
                     field->boundary_conditions.phys_fill_bc[+FieldDef::ConsIdx::Bx]);

    //----
    // D field divergence clean
    Real cd_scale = MFP::Larmor / (MFP::Debye * MFP::Debye * MFP::lightspeed);

    // assume zero charge density
    Amr& amr = *(mfp->get_parent());
    const Vector<BoxArray>& amr_grids = amr.boxArray(0, amr.finestLevel());
    const Vector<DistributionMapping>& amr_dmap = amr.DistributionMap(0, amr.finestLevel());
    const int nlevels = amr_grids.size();

    Vector<MultiFab> charge_density(nlevels);

    const int ngrow = 1;

    // go through level-by-level and calculate the charge and current density
    for (int ilev = 0; ilev < nlevels; ++ilev) {
    #ifdef AMREX_USE_EB
        EBFArrayBoxFactory const* amr_factory =
          mfp->get_eb_data(ilev, field->global_idx).ebfactory.get();
    #endif

        charge_density[ilev].define(amr_grids[ilev],
                                    amr_dmap[ilev],
                                    1,
                                    ngrow,
                                    MFInfo()
    #ifdef AMREX_USE_EB
                                      ,
                                    *amr_factory
    #endif
        );
        charge_density[ilev].setVal(0.0);

        AmrLevel& ilevel = amr.getLevel(ilev);
    #ifdef AMREX_USE_EB
        Vector<EBData>& eb_data = mfp->eb_data[ilev];
    #endif

        // get all of the source contributions
        for (const auto& src_idx : field->associated_actions) {
            if (mfp->actions[src_idx]->get_tag() != Plasma5::tag) continue;

            Plasma5& plasma = static_cast<Plasma5&>(*(mfp->actions[src_idx]));

            // now calculate the charge and current density
            const size_t num_src = plasma.species.size();

            // iterate over data
            for (MFIter mfi(ilevel.get_new_data(MFP::Cost_Idx)); mfi.isValid(); ++mfi) {
                const Box& box = mfi.tilebox();

                FArrayBox& local_cd = charge_density[ilev][mfi];
                local_cd.setVal(0.0);

                // get lists of the data for each state included in the src
                for (size_t src_idx = 0; src_idx < num_src; ++src_idx) {
                    HydroState& hydro = *plasma.species[src_idx];
                    const FArrayBox& cons = ilevel.get_new_data(hydro.data_idx)[mfi];
    #ifdef AMREX_USE_EB
                    const FArrayBox& vfrac = eb_data[hydro.global_idx].volfrac[mfi];
    #endif

                    hydro.calc_current_and_charge(box,
                                                  cons,
                                                  &local_cd,
                                                  nullptr

    #ifdef AMREX_USE_EB
                                                  ,
                                                  vfrac
    #endif
                    );
                }

                // scale by the relative permittivity and permeability

                FArrayBox& field_data = ilevel.get_new_data(field->data_idx)[mfi];

                local_cd.divide(field_data, box, +FieldDef::ConsIdx::ep, 0, 1);
            }
        }

        // scale the charge and current density appropriately
        charge_density[ilev].mult(-cd_scale, 0);

        //        plot_FAB_2d(charge_density[ilev],0, ngrow, "charge density "+num2str(ilev), false,
        //        true);
    }

    ParallelDescriptor::Barrier();

    // now that we have the charge density over all levels do the actual solve
    solve_divergence(mfp,
                     field->data_idx,
                     +FieldDef::ConsIdx::Dx,
                     +FieldDef::ConsIdx::phi,
                     field->boundary_conditions.phys_fill_bc[+FieldDef::ConsIdx::Dx],
                     GetVecOfPtrs(charge_density));

#endif
}

void Elliptic::solve_divergence(MFP* mfp,
                                const int state_idx,
                                const int vector_idx,
                                const int phi_idx,
                                const BCRec& bc,
                                Vector<MultiFab*> S_cc_ptr)
{
    BL_PROFILE("Elliptic::solve_divergence");

    Amr& amr = *(mfp->get_parent());
    const Geometry& geom = mfp->Geom();

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

    for (int ilev = 0; ilev < nlevels; ++ilev) {
        amr_factory[ilev] = mfp->get_eb_data(ilev, state_idx).ebfactory.get();
    }
#endif

    // Linear Solver
    LPInfo lp_info;

    MLNodeLaplacian matrix(amr_geom,
                           amr_grids,
                           amr_dmap,
                           lp_info
#ifdef AMREX_USE_EB
                           ,
                           amr_factory
#endif
    );

    // Set boundary conditions.
    // Note that Dirichlet boundary conditions are assumed to be homogeneous (i.e. phi = 0)

    // we assume that the BCRec given for the first component of the field vector
    // is representative and thus use it for determining the overall BC for phi

    Array<LinOpBCType, AMREX_SPACEDIM> lobc;
    Array<LinOpBCType, AMREX_SPACEDIM> hibc;

    set_field_bcs(lobc, hibc, bc, geom);

    matrix.setDomainBC(lobc, hibc);

    // Set matrix attributes to be used by MLMG solver
    matrix.setGaussSeidel(false);
    matrix.setHarmonicAverage(false);

    Vector<MultiFab> field_data(nlevels);
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
        S_cc_ptr.resize(nlevels, nullptr);
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            S_cc[ilev].define(amr_grids[ilev], amr_dmap[ilev], 1, 1, MFInfo());
            S_cc_ptr[ilev] = &S_cc[ilev];
            S_cc[ilev].setVal(0.0);  // Set it to zero for this example
        }
    }

    const Real time = amr.cumTime();

#ifdef AMREX_USE_EB
    State& istate = MFP::get_state(state_idx);
    EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate.eb2_index));
#endif

    for (int ilev = 0; ilev < nlevels; ++ilev) {
        AmrLevel& ilevel = amr.getLevel(ilev);

        //
        //  Create the cell-centered field we want to project
        //
        field_data[ilev]
          .define(amr_grids[ilev], amr_dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), ilevel.Factory());
        field_ptr[ilev] = &field_data[ilev];

        // Set field including ghost cells

        AmrLevel::FillPatch(ilevel,
                            field_data[ilev],
                            1,
                            time,
                            state_idx,
                            vector_idx,
                            AMREX_SPACEDIM,
                            0);

        // RHS is nodal
        const BoxArray& nd_grids =
          amrex::convert(amr_grids[ilev], IntVect::TheNodeVector());  // nodal grids

        // Multifab to host RHS
        rhs[ilev].define(nd_grids, amr_dmap[ilev], 1, 1, MFInfo());
        rhs[ilev].setVal(0.0);
        rhs_ptr[ilev] = &rhs[ilev];
        rhs_ptr2[ilev] = &rhs[ilev];

        // Node-centered contributions to RHS
        S_nd[ilev].define(nd_grids, amr_dmap[ilev], 1, 1, MFInfo());
        S_nd_ptr[ilev] = &S_nd[ilev];
        S_nd[ilev].setVal(0.0);  // Set it to zero for this example

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

    int v = amr.Verbose();
    nodal_solver.setVerbose(v);

    //
    // Solve div( sigma * grad(phi) ) = RHS
    //
    nodal_solver.solve(phi_ptr, rhs_ptr2, reltol, abstol);

    //
    // Create cell-centered multifab to hold value of -sigma*grad(phi) at cell-centers
    //
    Vector<MultiFab> fluxes(nlevels);
    Vector<MultiFab*> fluxes_ptr(nlevels, nullptr);
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        fluxes[ilev].define(amr_grids[ilev], amr_dmap[ilev], AMREX_SPACEDIM, 1, MFInfo());
        fluxes_ptr[ilev] = &fluxes[ilev];
        fluxes[ilev].setVal(0.0);
    }

    // Get fluxes from solver
    nodal_solver.getFluxes(fluxes_ptr);

    //
    // Apply projection explicitly --  vel = vel - sigma * grad(phi)
    //
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        AmrLevel& ilevel = amr.getLevel(ilev);
        MultiFab& ifield = ilevel.get_new_data(state_idx);
        MultiFab::Add(ifield, fluxes[ilev], 0, vector_idx, AMREX_SPACEDIM, 0);

        average_node_to_cellcenter(ifield, phi_idx, phi[ilev], 0, 1);
    }

    return;
}
