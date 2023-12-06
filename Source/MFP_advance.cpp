#include "MFP.H"
#include "MFP_action.H"
#include "MFP_diagnostics.H"
#include "MFP_eulerian.H"

Real MFP::advance(Real time, Real dt, int iteration, int ncycle)
{
    BL_PROFILE("MFP::advance");

    for (int i = 0; i < state.size(); ++i) {
        state[i].allocOldData();
        state[i].swapTimeLevels(dt);
    }

    MultiFab& C_new = get_new_data(Cost_Idx);
    C_new.setVal(0.0);

    // copy old into new
    for (const int& global_idx : eulerian_states) {
        EulerianState& istate = EulerianState::get_state_global(global_idx);
        const int data_idx = istate.data_idx;
        MultiFab& old_data = get_old_data(data_idx);
        MultiFab& new_data = get_new_data(data_idx);
        MultiFab::Copy(new_data, old_data, 0, 0, old_data.nComp(), 0);
    }

    // reset any flux registers
    const int finest_level = parent->finestLevel();
    for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
        EulerianState& istate = EulerianState::get_state(data_idx);
        if (istate.reflux && level < finest_level) {
            MFP& fine_level = getLevel(level + 1);
            fine_level.flux_reg[data_idx].reset();
        }
    }

    switch (time_integration_scheme) {
    case TimeIntegrator::RungeKutta: advance_RK(time, dt, iteration, ncycle); break;
    case TimeIntegrator::StrangSplitting: advance_strang(time, dt, iteration, ncycle); break;
#ifdef SYMPLECTIC
    case TimeIntegrator::Symplectic:
        // do nothing, handled by 'apply_change'
        break;
#endif
    default: advance_RK(time, dt, iteration, ncycle);
    }

    // apply any post timestep corrections
    for (const auto& act : actions) { act->apply_change(this, time, dt); }

    // do some lua garbage collection
    lua.script("collectgarbage('collect')");
    lua.script("collectgarbage('collect')");

    return dt;
}

void MFP::advance_RK(Real time, Real dt, int iteration, int ncycle)
{
    BL_PROFILE("MFP::advance_RK");

    Vector<UpdateData> RK_step(eulerian_states.size());

    // set the proportion of contributions for refluxing
    Real reflux_scaling;
    if (time_integration_nsteps >= 2) {
        reflux_scaling = 0.5 * dt;
    } else {
        reflux_scaling = dt;
    }

    for (const auto& act : actions) {
        // calculate any contributions to dU (dU += f(old))
        act->get_data(this, RK_step, time);
        act->calc_time_derivative(this, RK_step, time, dt);
        act->calc_spatial_derivative(this, RK_step, time, dt, reflux_scaling);
    }

    // add dU to new data (new = old + dU)
    for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
        if (RK_step[data_idx].dU_status != UpdateData::Status::Changed) continue;

        const int nc = RK_step[data_idx].dU.nComp();

        // U^1 = U^0 + dt*dU^0
        MultiFab& new_data = get_new_data(data_idx);
        MultiFab::Add(new_data, RK_step[data_idx].dU, 0, 0, nc, 0);

        RK_step[data_idx].U_status = UpdateData::Status::Expired;
        RK_step[data_idx].dU_status = UpdateData::Status::Expired;
    }

    if (time_integration_nsteps >= 2) {
        // Euler update (RK1) finished, now do RK2
        // now holds the Euler update based on time t0

        for (const auto& act : actions) {
            // load in the new data from the previous RK step
            // this is done again so that ghost cells are filled in
            act->get_data(this, RK_step, time + dt);

            act->calc_time_derivative(this, RK_step, time, dt);
            act->calc_spatial_derivative(this, RK_step, time, dt, reflux_scaling);
        }

        // final update
        for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
            if (RK_step[data_idx].dU_status != UpdateData::Status::Changed) continue;

            const int nc = RK_step[data_idx].dU.nComp();

            //            plot_FAB_2d(RK_step[data_idx].U, 1, 2, "U1", false,true);
            //            plot_FAB_2d(RK_step[data_idx].dU, 1, 2, "dU1", false,true);

            // U^2 = U^1 + dt*dU^1
            MultiFab::Add(RK_step[data_idx].U, RK_step[data_idx].dU, 0, 0, nc, 0);

            MultiFab& old_data = get_old_data(data_idx);
            MultiFab& new_data = get_new_data(data_idx);

            // the updated solution
            // U^n+1 = 0.5*U^0 + 0.5*U^2
            //       = 0.5*U^0 + 0.5*(U^1 + dt*dU^1)
            //       = 0.5*U^0 + 0.5*(U^0 + dt*dU^0 + dt*dU^1)
            //       = U^0 + 0.5*dt*dU^0 + 0.5*dt*dU^1
            MultiFab::LinComb(new_data, 0.5, old_data, 0, 0.5, RK_step[data_idx].U, 0, 0, nc, 0);
        }
    }

    // update new data with any one-shot updates based on old data (new = f(old))
    for (const auto& act : actions) {
        act->apply_time_derivative(this, time, dt);
        act->apply_spatial_derivative(this, time, dt);
    }
}

void MFP::advance_strang(Real time, Real dt, int iteration, int ncycle)
{
    BL_PROFILE("MFP::advance_strang");

    Vector<UpdateData> RK_step(eulerian_states.size());

    ///////////////////////////////////////////////////////////////////////////
    /// SOURCE TERMS

    for (const auto& act : actions) {
        // calculate any contributions to dU (dU += f(old))
        act->get_data(this, RK_step, time);
        act->calc_time_derivative(this, RK_step, time, 0.5 * dt);
    }

    // add dU to new data (new = old + dU)
    for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
        if (RK_step[data_idx].dU_status != UpdateData::Status::Changed) continue;

        const int nc = RK_step[data_idx].dU.nComp();

        // U^1 = U^0 + dt*dU^0
        if (time_integration_nsteps >= 2) {
            MultiFab::Add(RK_step[data_idx].U, RK_step[data_idx].dU, 0, 0, nc, 0);
            RK_step[data_idx].dU.setVal(0.0);
            RK_step[data_idx].dU_status = UpdateData::Status::Expired;
        } else {
            MultiFab& new_data = get_new_data(data_idx);
            MultiFab::Add(new_data, RK_step[data_idx].dU, 0, 0, nc, 0);
        }
    }

    if (time_integration_nsteps >= 2) {
        // Euler update (RK1) finished, now do RK2
        // RK_step now holds the Euler update based on time t0

        for (const auto& act : actions) {
            act->calc_time_derivative(this, RK_step, time, 0.5 * dt);
        }

        for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
            if (RK_step[data_idx].dU_status != UpdateData::Status::Changed) continue;

            const int nc = RK_step[data_idx].dU.nComp();

            // U^2 = U^1 + dt*dU^1
            MultiFab::Add(RK_step[data_idx].U, RK_step[data_idx].dU, 0, 0, nc, 0);
            RK_step[data_idx].dU.setVal(0.0);

            MultiFab& old_data = get_old_data(data_idx);
            MultiFab& new_data = get_new_data(data_idx);

            // the updated solution
            // U^n+1 = 0.5*U^0 + 0.5*U^2
            //       = 0.5*U^0 + 0.5*(U^1 + dt*dU^1)
            //       = 0.5*U^0 + 0.5*(U^0 + dt*dU^0 + dt*dU^1)
            //       = U^0 + 0.5*dt*dU^0 + 0.5*dt*dU^1
            MultiFab::LinComb(new_data, 0.5, old_data, 0, 0.5, RK_step[data_idx].U, 0, 0, nc, 0);
            RK_step[data_idx].U_status = UpdateData::Status::Expired;
            RK_step[data_idx].dU_status = UpdateData::Status::Expired;
        }
    }

    for (const auto& act : actions) {
        // update new data with any one-shot updates based on old data (new* = f(old))
        act->apply_time_derivative(this, time, 0.5 * dt);
    }

    ///////////////////////////////////////////////////////////////////////////
    /// FLUXES

    for (const auto& act : actions) {
        act->get_data(this, RK_step, time + dt);

        // calculate any contributions to dU (dU += f(new*))

        act->calc_spatial_derivative(this, RK_step, time + dt, dt, dt);

        // update new data with any one-shot updates based on new* data (new** = f(new*))
        act->apply_spatial_derivative(this, time + dt, dt);
    }

    // add dU to new data (new** = new* + dU)
    for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
        if (RK_step[data_idx].dU_status != UpdateData::Status::Changed) continue;

        const int nc = RK_step[data_idx].dU.nComp();

        MultiFab& new_data = get_new_data(data_idx);
        MultiFab::Add(new_data, RK_step[data_idx].dU, 0, 0, nc, 0);

        RK_step[data_idx].U_status = UpdateData::Status::Expired;
        RK_step[data_idx].dU_status = UpdateData::Status::Expired;
    }

    ///////////////////////////////////////////////////////////////////////////
    /// SOURCE TERMS

    for (const auto& act : actions) {
        // calculate any contributions to dU (dU += f(old))
        act->get_data(this, RK_step, time + dt);
        act->calc_time_derivative(this, RK_step, time + dt, 0.5 * dt);
    }

    // add dU to new data (new = old + dU)
    for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
        if (RK_step[data_idx].dU_status != UpdateData::Status::Changed) continue;

        const int nc = RK_step[data_idx].dU.nComp();

        // U^1 = U^0 + dt*dU^0
        if (time_integration_nsteps >= 2) {
            MultiFab::Add(RK_step[data_idx].U, RK_step[data_idx].dU, 0, 0, nc, 0);
            RK_step[data_idx].dU.setVal(0.0);
            RK_step[data_idx].dU_status = UpdateData::Status::Inactive;
        } else {
            MultiFab& new_data = get_new_data(data_idx);
            MultiFab::Add(new_data, RK_step[data_idx].dU, 0, 0, nc, 0);
        }
    }

    // Euler update (RK1) finished, now do RK2
    // RK_step now holds the Euler update based on time t1

    if (time_integration_nsteps >= 2) {
        for (const auto& act : actions) {
            act->calc_time_derivative(this, RK_step, time + dt, 0.5 * dt);
        }

        for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {
            if (RK_step[data_idx].dU_status != UpdateData::Status::Changed) continue;

            const int nc = RK_step[data_idx].dU.nComp();

            // U^2 = U^1 + dt*dU^1
            MultiFab::Add(RK_step[data_idx].U, RK_step[data_idx].dU, 0, 0, nc, 0);
            RK_step[data_idx].dU.setVal(0.0);

            MultiFab& new_data = get_new_data(data_idx);

            // the updated solution
            // U^n+1 = 0.5*U^0 + 0.5*U^2
            MultiFab::LinComb(new_data, 0.5, new_data, 0, 0.5, RK_step[data_idx].U, 0, 0, nc, 0);
        }
    }

    for (const auto& act : actions) {
        // update new data with any one-shot updates based on old data (new* = f(old))
        act->apply_time_derivative(this, time + dt, 0.5 * dt);
    }

    ParallelDescriptor::Barrier();
}
