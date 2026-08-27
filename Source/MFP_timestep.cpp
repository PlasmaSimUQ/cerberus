#include "MFP.H"
#include "MFP_state.H"

Real MFP::estTimeStep()
{
    BL_PROFILE("MFP::estTimeStep");

    if (force_dt > 0.0) { return force_dt; }

    Real estdt = std::numeric_limits<Real>::max();

    // handle all of the intrinsic time step limitations of a state
    for (auto& state : states) { estdt = std::min(estdt, state->get_allowed_time_step(this)); }

    for (auto& action : actions) { estdt = std::min(estdt, action->get_allowed_time_step(this)); }

    estdt *= cfl;
    ParallelDescriptor::ReduceRealMin(estdt);

    return estdt;
}

void MFP::computeInitialDt(int finest_level,
                           int sub_cycle,
                           Vector<int>& n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>& dt_level,
                           Real stop_time)
{
    BL_PROFILE("MFP::computeInitialDt");

    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0) { return; }

    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        dt_level[i] = getLevel(i).estTimeStep();
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0, n_factor * dt_level[i]);
    }

    //
    // Optional gentle start (Lua `dt_init_shrink`, 1 = off). This is applied
    // ONCE, here: AMReX calls computeInitialDt only when levelSteps(0) == 0
    // (AMReX_Amr.cpp:2101), so it never fires again during the run and never
    // on restart -- a restarted run resumes from the checkpointed dt_level.
    // computeNewDt then grows dt by at most dt_change_max per coarse step.
    // Skipped when the step is forced, since estTimeStep returns force_dt
    // before cfl is applied and a forced dt must not be scaled.
    //
    if ((force_dt <= 0.0) && (dt_init_shrink > 0.0)) { dt_0 *= dt_init_shrink; }

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

void MFP::computeNewDt(int finest_level,
                       int sub_cycle,
                       Vector<int>& n_cycle,
                       const Vector<IntVect>& ref_ratio,
                       Vector<Real>& dt_min,
                       Vector<Real>& dt_level,
                       Real stop_time,
                       int post_regrid_flag)
{
    BL_PROFILE("MFP::computeNewDt");

    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0) { return; }

    for (int i = 0; i <= finest_level; i++) { dt_min[i] = getLevel(i).estTimeStep(); }

    if (post_regrid_flag == 1) {
        //
        // Limit dt's by pre-regrid dt
        //
        for (int i = 0; i <= finest_level; i++) { dt_min[i] = std::min(dt_min[i], dt_level[i]); }
    } else {
        //
        // Limit dt's by dt_change_max * old dt (Lua `dt_change_max`)
        //
        for (int i = 0; i <= finest_level; i++) {
            dt_min[i] = std::min(dt_min[i], dt_change_max * dt_level[i]);
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
        if ((cur_time + dt_0) > (stop_time - eps)) { dt_0 = stop_time - cur_time; }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}
