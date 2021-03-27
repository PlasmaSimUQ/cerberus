
#include "MFP.H"

#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#endif

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include "MFP_diagnostics.H"

#include <climits>
#include <fstream>

using namespace amrex;

constexpr int MFP::level_mask_interior;
constexpr int MFP::level_mask_covered;
constexpr int MFP::level_mask_notcovered;
constexpr int MFP::level_mask_physbnd;

bool MFP::archive_checkpoint = true;

bool MFP::first_step = true;

#ifdef AMREX_PARTICLES
int MFP::particle_verbosity = 0;
#endif

GlobalData MFP::gd;

MFP::MFP() {}

MFP::MFP(Amr& papa, int lev, const Geometry& level_geom, const BoxArray& bl,
         const DistributionMapping& dm, Real time)
    : AmrLevel(papa, lev, level_geom, bl, dm, time) {
    if (level > 0) {

        flux_reg.resize(gd.num_solve_state);
        for (int idx = 0; idx < gd.num_solve_state; ++idx) {
            State &istate = gd.get_state(idx);
            if (istate.reflux) {
                flux_reg[idx].define(
                            bl, papa.boxArray(level - 1),
                            dm, papa.DistributionMap(level - 1),
                            level_geom, papa.Geom(level - 1),
                            papa.refRatio(level - 1), level,
                            desc_lst[idx].nComp());
            }
        }
    }

    buildMetrics();
}

MFP::~MFP(){}

void MFP::init(AmrLevel& old) {
    auto& oldlev = dynamic_cast<MFP&>(old);

    Real dt_new = parent->dtLevel(level);
    Real cur_time = oldlev.state[0].curTime();
    Real prev_time = oldlev.state[0].prevTime();
    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    for (int idx = 0; idx<gd.num_solve_state; ++idx) {
        MultiFab& S_new = get_new_data(idx);
        int nc = gd.states[idx]->n_cons();
        FillPatch(old, S_new, 0, cur_time, idx, 0, nc);
    }

    MultiFab& C_new = get_new_data(gd.Cost_Idx);
    FillPatch(old, C_new, 0, cur_time, gd.Cost_Idx, 0, 1);

    if (gd.Shock_Idx > 0) {
        MultiFab& S_new = get_new_data(gd.Shock_Idx);
        FillPatch(old, S_new, 0, cur_time, gd.Shock_Idx, 0, gd.num_shock_detector);
    }

}

void MFP::init() {
    Real dt = parent->dtLevel(level);
    Real cur_time = getLevel(level - 1).state[0].curTime();
    Real prev_time = getLevel(level - 1).state[0].prevTime();
    Real dt_old = (cur_time - prev_time) /
                  static_cast<Real>(parent->MaxRefRatio(level - 1));
    setTimeLevel(cur_time, dt_old, dt);

    for (int idx = 0; idx<gd.num_solve_state; ++idx) {
        MultiFab& S_new = get_new_data(idx);
        int nc = gd.states[idx]->n_cons();
        FillCoarsePatch(S_new, 0, cur_time, idx, 0, nc);
    }

    MultiFab& C_new = get_new_data(gd.Cost_Idx);
    FillCoarsePatch(C_new, 0, cur_time, gd.Cost_Idx, 0, 1);

    if (gd.Shock_Idx > 0) {
        MultiFab& S_new = get_new_data(gd.Shock_Idx);
        FillCoarsePatch(S_new, 0, cur_time, gd.Shock_Idx, 0, gd.num_shock_detector);
    }

}

void MFP::initData() {
    BL_PROFILE("MFP::initData()");

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    // initialise states
    for (int idx = 0; idx < gd.num_solve_state; ++idx) {

        MultiFab& S_new = get_new_data(idx);
        //    Real cur_time = state[si].curTime();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();

#ifdef AMREX_USE_EB
            const EBCellFlagFab& flag = getEBData(idx).flags[mfi];
#endif
            init_data(box,
                      EB_OPTIONAL(flag,)
                      S_new[mfi],
                      dx,
                      prob_lo,
                      idx);
        }
    }

    // initialise cost state
    MultiFab& C_new = get_new_data(gd.Cost_Idx);
    C_new.setVal(1.0);

    // initialise particles
#ifdef AMREX_PARTICLES

    ParmParse pp("particles");
    pp.query("v", particle_verbosity);

    if (level == 0) {
        init_particles();
    }
#endif

    if (gd.Shock_Idx > 0) {
        MultiFab& S_new = get_new_data(gd.Shock_Idx);
        S_new.setVal(0.0);
    }

}

void MFP::computeInitialDt(int finest_level, int sub_cycle,
                           Vector<int>& n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>& dt_level, Real stop_time) {
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0) {
        return;
    }

    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0, n_factor * dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
    Real cur_time = state[0].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps)) dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

void MFP::computeNewDt(int finest_level, int sub_cycle, Vector<int>& n_cycle,
                       const Vector<IntVect>& ref_ratio, Vector<Real>& dt_min,
                       Vector<Real>& dt_level, Real stop_time,
                       int post_regrid_flag) {
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0) {
        return;
    }

    for (int i = 0; i <= finest_level; i++) {
        dt_min[i] = getLevel(i).estTimeStep();
    }

    if (post_regrid_flag == 1) {
        //
        // Limit dt's by pre-regrid dt
        //
        for (int i = 0; i <= finest_level; i++) {
            dt_min[i] = std::min(dt_min[i], dt_level[i]);
        }
    } else {
        //
        // Limit dt's by change_max * old dt
        //
        static Real change_max = 1.1;
        for (int i = 0; i <= finest_level; i++) {
            dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
        }
    }

    //
    // Find the minimum over all levels
    //
    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0, n_factor * dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
    Real cur_time = state[0].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
        }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

void MFP::post_regrid(int lbase, int new_finest) {

#ifdef AMREX_PARTICLES
    for (AmrTracerParticleContainer* TracerPC : particles) {
        TracerPC->Redistribute();
        //    if (TracerPC && level == lbase) {
        //        TracerPC->Redistribute(lbase);
        //    }
    }
#endif

}

void MFP::post_timestep(int iteration) {
    // reflux
    if (level < parent->finestLevel()) {
        MFP& fine_level = getLevel(level + 1);

        for (int idx = 0; idx < gd.num_solve_state; ++idx) {
            State &istate = gd.get_state(idx);
            if (istate.reflux) {
                MultiFab& S_crse = get_new_data(idx);
#ifdef AMREX_USE_EB
                MultiFab& S_fine = fine_level.get_new_data(idx);
                fine_level.flux_reg[idx].Reflux(S_crse, getEBData(idx).volfrac, S_fine,
                                                getFineEBData(idx).volfrac);
#else
                fine_level.flux_reg[idx].Reflux(S_crse, 0);
#endif
            }
        }
    }

    if (level < parent->finestLevel()) {
        avgDown();
    }

#ifdef AMREX_PARTICLES
    if (gd.do_tracer_particles)
    {
        const int ncycle = parent->nCycle(level);
        //
        // Don't redistribute on the final subiteration except on the coarsest grid.
        //
        if (iteration < ncycle || level == 0)
        {
            int ngrow = (level == 0) ? 0 : iteration;

            for (AmrTracerParticleContainer* TracerPC : particles) {
                TracerPC->Redistribute(level, parent->finestLevel(), ngrow);
            }
        }
    }
#endif
}

void MFP::postCoarseTimeStep(Real time) {

    // perform projection of divergence for fields
//    correct_dynamic_fields(time);
    project_divergence(time);

    if (gd.verbose >= 1) {
        amrex::Print().SetPrecision(17)
                << "[MFP] : step = " << parent->levelSteps(0) << ", time = " << time
                << std::endl;
    }

    first_step = false;

    // This only computes sum on level 0
    if (gd.verbose >= 3) {
        printTotal();
    }
}

void MFP::printTotal()  {
    for (int idx = 0; idx < gd.num_solve_state; ++idx) {
        State &istate = gd.get_state(idx);
        const MultiFab& S_new = get_new_data(idx);
        MultiFab mf(grids, dmap, 1, 0);
        int nc = istate.n_cons();
        Vector<Real> tot(nc);
        for (int comp = 0; comp < nc; ++comp) {
            MultiFab::Copy(mf, S_new, comp, 0, 1, 0);
#ifdef AMREX_USE_EB
            MultiFab::Multiply(mf, getEBData(idx).volfrac, 0, 0, 1, 0);
#endif
            tot[comp] = mf.sum(0, true) * geom.ProbSize();
        }
        ParallelDescriptor::ReduceRealSum(tot.data(), nc, ParallelDescriptor::IOProcessorNumber());

        std::stringstream account;

        account << "\n[MFP] " << istate.name << "\n";

        for (int comp = 0; comp < nc; ++comp) {
            std::string name = istate.get_cons_names()[comp];
            account << "      Total " << name << " is " << tot[comp];

            if (!first_step) {
                account << " (|\u0394|=" << std::abs(istate.sum_cons[comp] -  tot[comp]) << ")";
            }

            account << "\n";

            if (first_step) {
                istate.sum_cons[comp] = tot[comp];
            }
        }

        amrex::Print().SetPrecision(17) << account.str();

    }
}

void MFP::post_init(Real) {

    // check our reconstruction doesn't go out-of-bounds

    const Box& box = geom.Domain();
    const IntVect size = box.bigEnd() - box.smallEnd();
    for (int idx=0; idx<gd.num_solve_state; ++idx) {
        State& istate = gd.get_state(idx);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            if (istate.reconstruction->num_grow > size[d])
                Abort("Reconstruction stencil is larger than domain in dim="+num2str(d));
        }
    }

    //#ifdef AMREX_PARTICLES
    //    for (AmrTracerParticleContainer* TracerPC : particles) {
    //        TracerPC->Redistribute();
    //    }
    //#endif

    if (level > 0) return;
    for (int k = parent->finestLevel() - 1; k >= 0; --k) {
        getLevel(k).avgDown();
    }

    solve_static_fields(parent->cumTime());

    if (gd.verbose >= 3) {
        printTotal();
    }
}

void MFP::post_restart() {
    if (level > 0) {
        flux_reg.resize(gd.num_solve_state);
        for (int idx = 0; idx < gd.num_solve_state; ++idx) {
            State &istate = gd.get_state(idx);
            if (istate.reflux) {
                flux_reg[idx].define(
                            grids, parent->boxArray(level - 1), dmap,
                            parent->DistributionMap(level - 1), geom,
                            parent->Geom(level - 1),
                            parent->refRatio(level - 1), level,
                            desc_lst[idx].nComp());
            }
        }


    }


    buildMetrics();

#ifdef AMREX_PARTICLES
    std::string restart_chkfile;
    ParmParse pp("amr");
    pp.query("restart", restart_chkfile);
    ParticlePostRestart(restart_chkfile);
#endif

    solve_static_fields(parent->cumTime());

}

void MFP::get_lua_script() {
    gd.lua_script =
        #include "default.lua"
            ;

    ParmParse pp("mfp");
    std::string input;

    if (pp.contains("lua")) {
        pp.get("lua", input);
    } else {
        amrex::Abort("Exiting as we have not defined the problem (mfp.lua = '')");
    }

    gd.lua_script += input;

    // insert call to preprocessor defined in default.lua
    gd.lua_script += "\npreprocess()\n";

    return;
}

void MFP::save_lua_script() {

    if (ParallelDescriptor::IOProcessor()) {

        std::ofstream ofs;

        ofs.open("MFP_Lua_initialisation.lua", std::ofstream::out | std::ofstream::trunc);
        ofs << gd.lua_script;
        ofs.close();
        amrex::Print() << "'MFP_Lua_initialisation.lua' written to file\n";

    }

    return;
}

void MFP::avgDown() {
    BL_PROFILE("MFP::avgDown()");

    if (level == parent->finestLevel()) return;

    auto& fine_lev = getLevel(level + 1);

    for (int idx = 0; idx < gd.num_solve_state; ++idx) {
        MultiFab& S_crse = get_new_data(idx);
        MultiFab& S_fine = fine_lev.get_new_data(idx);

#ifdef AMREX_USE_EB
        MultiFab volume(S_fine.boxArray(), S_fine.DistributionMap(), 1, 0);
        volume.setVal(1.0);
        amrex::EB_average_down(S_fine, S_crse, volume, getFineEBData(idx).volfrac, 0,
                               S_fine.nComp(), fine_ratio);
#else
        amrex::average_down(S_fine, S_crse, 0, S_fine.nComp(), fine_ratio);
#endif

    }
}


void MFP::buildMetrics() {
#ifdef AMREX_USE_EB
    BL_PROFILE("MFP::buildMetrics()");

    // make sure dx == dy == dz
    const Real* dx = geom.CellSize();
    for (int i = 1; i < AMREX_SPACEDIM; ++i) {
        if (std::abs(dx[0] - dx[i]) > 1.e-12 * dx[0]) {
            amrex::Abort("MFP: must have dx == dy == dz\n");
        }
    }

    Vector<EBData>& eb_data = getEBData();

    // get the information for each embedded boundary
    eb_data.resize(gd.num_solve_state);
    const Vector<int> ngrow = {m_eb_basic_grow_cells,m_eb_volume_grow_cells,m_eb_full_grow_cells};
    for (int idx=0; idx<gd.num_solve_state; ++idx) {
        State &istate = gd.get_state(idx);

        eb_data[idx].ebfactory = makeEBFabFactory (istate.eb2_index,
                                                   geom,
                                                   grids,
                                                   dmap,
                                                   ngrow,
                                                   EBSupport::full);


        const auto& flags_orig = eb_data[idx].ebfactory->getMultiEBCellFlagFab();
        eb_data[idx].flags.define(grids, dmap, 1, flags_orig.n_grow);
        eb_data[idx].flags.copy(flags_orig,0,0,1,flags_orig.n_grow,flags_orig.n_grow,geom.periodicity());

        const auto& vfrac_original = eb_data[idx].ebfactory->getVolFrac();
        eb_data[idx].volfrac.define(grids, dmap, 1, vfrac_original.n_grow);
        eb_data[idx].volfrac.copy(vfrac_original,0,0,1,vfrac_original.n_grow,vfrac_original.n_grow,geom.periodicity());

        eb_data[idx].bndrycent = &(eb_data[idx].ebfactory->getBndryCent());
        eb_data[idx].bndrynorm = &(eb_data[idx].ebfactory->getBndryNormal());
        eb_data[idx].areafrac = eb_data[idx].ebfactory->getAreaFrac();
        eb_data[idx].facecent = eb_data[idx].ebfactory->getFaceCent();

        // different boundary conditions
        eb_data[idx].bndryidx.define(grids, dmap, 1, m_eb_basic_grow_cells, flags_orig);
        eb_data[idx].bndryidx.setVal(-1.0);

        // update ghost cells
        auto& flags = eb_data[idx].flags;
        auto& vfrac = eb_data[idx].volfrac;

        for (MFIter mfi(vfrac); mfi.isValid(); ++mfi){

//            plot_FAB_2d(flags[mfi], "flags before", false);
//            plot_FAB_2d(vfrac[mfi], 0, "vfrac before", false,false);

            istate.update_eb_vfrac(geom,vfrac[mfi]);
            istate.update_eb_flags(geom,flags[mfi]);

//            plot_FAB_2d(flags[mfi], "flags after", false);
//            plot_FAB_2d(vfrac[mfi], 0, "vfrac after", false,true);
        }
    }



    // loop over all of the individual EB definitions and put an id into bndryidx wherever there is
    // a boundary
    for (const auto &eb : gd.eb_def) {

        std::unique_ptr<EBFArrayBoxFactory> bc_ebfactory = makeEBFabFactory(eb.index_space,
                                                                            geom,
                                                                            grids,
                                                                            dmap,
                                                                            ngrow,
                                                                            EBSupport::full);

        const FabArray<EBCellFlagFab>& bc_flags = bc_ebfactory->getMultiEBCellFlagFab();

        for (const auto& si : eb.states) {
            MultiCutFab& bndryidx = eb_data[si.first].bndryidx;
            const FabArray<EBCellFlagFab>& state_flags = eb_data[si.first].flags;

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


Real MFP::estTimeStep() {
    BL_PROFILE("MFP::estTimeStep()");

    if (gd.force_dt > 0.0) {
        return gd.force_dt;
    }

    Real estdt = std::numeric_limits<Real>::max();

    const Real* dx = geom.CellSize();
    Real dt = estdt;

    // get a list of all the datasets
    Vector<const MultiFab*> grab_mf(gd.num_solve_state);
    for (int idx=0; idx<gd.num_solve_state; ++idx) {
        grab_mf[idx] = &get_new_data(idx);
    }

    // now iterate over datasets and grab the same box of data
    // from each dataset as we go
    Vector<const FArrayBox*> grab_fbox(gd.num_solve_state);


#ifdef AMREX_USE_EB
    Vector<const FArrayBox*> grab_vfrac(gd.num_solve_state);
#endif

    for (MFIter mfi(*grab_mf[0]); mfi.isValid(); ++mfi) {

        const Box& box = mfi.tilebox();


        for (int idx=0; idx<gd.num_solve_state; ++idx) {

            const MultiFab &S = *grab_mf[idx];

            grab_fbox[idx] = nullptr;

            FabType flag = FabType::regular;

#ifdef AMREX_USE_EB
            // get the EB data required for later calls
            const FArrayBox& vfrac = getEBData(idx).volfrac[mfi];

            flag = vfrac.getType();

            if (flag != FabType::covered) {
                grab_vfrac[idx] = &vfrac;
            } else {
                grab_vfrac[idx] = nullptr;
            }

#endif

            if (flag != FabType::covered) {
                if (is_vacuum(box, S[mfi], idx)) {
                    flag = FabType::covered;
                }
            }

            if (flag != FabType::covered) {
                grab_fbox[idx] = &S[mfi];
            } else {
                grab_fbox[idx] = nullptr;
            }
        }

        // now that we have a full set of data we can pass it off to the
        // function that does the hard work

        estimate_dt(box,
                    grab_fbox,
                    EB_OPTIONAL(grab_vfrac,)
                    dx,
                    dt);

        estdt = std::min(estdt, dt);


    }

    estdt *= gd.cfl;
    ParallelDescriptor::ReduceRealMin(estdt);


    // update any divergence cleaning factors
    for (auto &istate : gd.states) {
        istate->update_div_clean(dx);
    }

    return estdt;
}

Real MFP::initialTimeStep() { return estTimeStep(); }

#ifdef AMREX_USE_EB
// https://www.asawicki.info/news_1428_finding_polygon_of_plane-aabb_intersection

bool ray_to_plane(const RealVect &ray_origin, const RealVect &ray_direction, const Array<Real, AMREX_SPACEDIM+1>& plane, Real &T) {
    Real VD = 0.0;
    
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        VD += plane[i]*ray_direction[i];
    }

    if (VD == 0.0) {
        return false;
    }

    T = 0.0;
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        T += plane[i]*ray_origin[i];
    }

    T += plane[AMREX_SPACEDIM];

    T *= - 1.0 / VD;

    return true;
}

Vector<RealVect> plane_cell_intersections(const Array<Real, AMREX_SPACEDIM+1>& plane) {
    Vector<RealVect> intersect;
    intersect.reserve(6); // max of 6 intersection points in 3D

    
    Real t;
    RealVect orig;

    // Test edges along X axis, pointing right.
    RealVect dir(AMREX_D_DECL(1.0, 0.0, 0.0));

    orig = RealVect(AMREX_D_DECL(-0.5, -0.5, -0.5));
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);
    
    orig = RealVect(AMREX_D_DECL(-0.5, 0.5, -0.5));
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);

#if AMREX_SPACEDIM == 3
    orig = RealVect(-0.5, -0.5, 0.5);
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);

    orig = RealVect(-0.5, 0.5, 0.5);
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);
#endif

    // Test edges along Y axis, pointing up.
    dir = RealVect(AMREX_D_DECL(0.0, 1.0, 0.0));

    orig = RealVect(AMREX_D_DECL(-0.5, -0.5, -0.5));
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);

    orig = RealVect(AMREX_D_DECL(0.5, -0.5, -0.5));
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);

#if AMREX_SPACEDIM == 3
    orig = RealVect(-0.5, -0.5, 0.5);
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);

    orig = RealVect(0.5, -0.5, 0.5);
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);


    // Test edges along Z axis, pointing forward.
    dir = RealVect(0.0, 0.0, 1.0);

    orig = RealVect(-0.5, -0.5, -0.5);
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);

    orig = RealVect(0.5, -0.5, -0.5);
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);

    orig = RealVect(-0.5, 0.5, -0.5);
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);

    orig = RealVect(0.5, 0.5, -0.5);
    if (ray_to_plane(orig, dir, plane, t) && t >= 0.0 && t <= 1.0)
        intersect.push_back(orig + dir * t);
#endif



    // sort the intersection points
#if AMREX_SPACEDIM == 2
    std::sort(intersect.begin(), intersect.end(), [&](const RealVect &lhs, const RealVect &rhs) -> bool {
        return lhs[0]*rhs[1] - rhs[0]*lhs[1] < 0.0;
    } );
#elif AMREX_SPACEDIM == 3
    const RealVect plane_normal(AMREX_D_DECL(plane[0], plane[1], plane[2]));
    const RealVect& origin = intersect[0];
    std::sort(intersect.begin(), intersect.end(), [&](const RealVect &lhs, const RealVect &rhs) -> bool {
        RealVect v = lhs - origin;
        v.crossProduct(rhs - origin);
        return v.dotProduct(plane_normal) < 0.0;
    } );

#endif

    // make sure we don't have any duplicates
    const Real tol = 1e-12;
    intersect.erase( std::unique( intersect.begin(), intersect.end(),  [&](const RealVect &lhs, const RealVect &rhs) -> bool {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            if (std::abs(lhs[i] - rhs[i]) > tol) {
                return false;
            }
        }
        return true;
    }), intersect.end() );

    return intersect;
}

FArrayBox MFP::get_eb_faces(const MFIter &mfi, const int idx) {

    //    const Box& box = mfi.fabbox();

    Vector<EBData>& eb_data = getEBData();

    const FArrayBox &bcent = (*eb_data[idx].bndrycent)[mfi];
    const FArrayBox &bnorm = (*eb_data[idx].bndrynorm)[mfi];
    const EBCellFlagFab& flag = eb_data[idx].flags[mfi];

    const Box& box = bnorm.box();

    Array4<const Real> const& bcent4 = bcent.array();
    Array4<const Real> const& bnorm4 = bnorm.array();
    Array4<const EBCellFlag> const& flag4 = flag.array();

    FArrayBox polys(box, AMREX_D_PICK(0, 5, 19));
    Array4<Real> const& poly4 = polys.array();

    polys.setVal(-1.0);

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    // iterate over faces
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

                const EBCellFlag &flag_eb = flag4(i,j,k);

                if (flag_eb.isSingleValued()) {

                    Vector<RealVect> pts;
                    Array<Real, AMREX_SPACEDIM+1> plane;

                    // get the coefficients to define a plane
                    plane.fill(0.0);
                    for (int di=0; di<AMREX_SPACEDIM; ++di) {
                        plane[di] = bnorm4(i,j,k,di);
                        plane[AMREX_SPACEDIM] -= plane[di]*bcent4(i,j,k,di);
                    }

                    // calculate all of the intersections between the plane and the cell edges
                    pts = plane_cell_intersections(plane);

                    int cnt = 1;
                    poly4(i,j,k,0) = pts.size(); // first load in how many points there are
                    for (const RealVect &v : pts) {
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            poly4(i,j,k,cnt) = v[d];
                            cnt++;
                        }
                    }

                }
            }
        }
    }

    return polys;
}
#endif

bool MFP::is_vacuum(const Box& bx, const FArrayBox& fab, const int idx)
{
    // check if density is effectively zero and thus no fluxes
    // need to be calculated

    int di = gd.states[idx]->get_cons_density_idx();

    // check that the state type corresponds to something that has density
    if (di < 0) {
        return false;
    }

    if (fab.max(bx, di) <= gd.effective_zero) {
        return true;
    } else {
        return false;
    }
}
