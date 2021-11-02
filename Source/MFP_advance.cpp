#include "MFP.H"
#include "MFP_eulerian.H"
#include "MFP_action.H"
#include "MFP_diagnostics.H"



Real MFP::advance(Real time, Real dt, int iteration, int ncycle)
{
    BL_PROFILE("MFP::advance()");

    for (int i = 0; i < state.size(); ++i) {
        state[i].allocOldData();
        state[i].swapTimeLevels(dt);
    }

    MultiFab& C_new = get_new_data(Cost_Idx);
    C_new.setVal(0.0);

    // sety up the scratch space for derivatives of the eulerian data
    if (eulerian_du.empty() && need_scratch_space) {
        eulerian_du.resize(eulerian_states.size());

#ifdef AMREX_USE_EB
        const int num_grow_eb = 2;
#else
        const int num_grow_eb = 0;
#endif

        for (int data_idx=0; data_idx<eulerian_states.size(); ++data_idx) {
            int nc = state[data_idx].descriptor()->nComp();
            eulerian_du[data_idx].first = 0;
            eulerian_du[data_idx].second.define(grids, dmap, nc, num_grow_eb, MFInfo(), Factory());
        }
    }

    // copy old into new
    for (const int& global_idx : eulerian_states) {
        EulerianState& istate = EulerianState::get_state_global(global_idx);
        const int data_idx = istate.data_idx;
        MultiFab& old_data = get_old_data(data_idx);
        MultiFab& new_data = get_new_data(data_idx);
        MultiFab::Copy(new_data, old_data,0,0,old_data.nComp(),0);
    }

    switch (time_integration_scheme) {
    case TimeIntegrator::OneStep:
        advance_euler(time, dt, iteration, ncycle);
        break;
    case TimeIntegrator::StrangSplitting:
        advance_strang(time, dt, iteration, ncycle);
        break;
#ifdef SYMPLECTIC
    case TimeIntegrator::Symplectic:
        // do nothing, handled by 'apply_change'
        break;
#endif
    default:
        advance_euler(time, dt, iteration, ncycle);
    }

    // apply any post timestep corrections
    for (const auto& act : actions) {
        act->apply_change(this,time, dt);
    }

    // do some lua garbage collection
    lua.script("collectgarbage('collect')");
    lua.script("collectgarbage('collect')");

    return dt;
}

void MFP::apply_derivative(Vector<std::pair<int,MultiFab>>& dU)
{
    for (int data_idx = 0; data_idx < eulerian_states.size(); ++data_idx) {

        if (!dU[data_idx].first) continue;

//        plot_FAB_2d(dU[data_idx].second,1, 0, "dU[1] "+EulerianState::get_state(data_idx).name, false, true);

        MultiFab& new_data = get_new_data(data_idx);

        const int nc = dU[data_idx].second.nComp();

        MultiFab::Add(new_data, dU[data_idx].second, 0, 0, nc, 0);

        // zero it out now that it has been applied
        dU[data_idx].first = 0;
        dU[data_idx].second.setVal(0.0);
    }
}

void MFP::advance_euler(Real time, Real dt, int iteration, int ncycle)
{
    BL_PROFILE("MFP::advance_euler");

    // zero out the derivative accumulator
    for (auto& dU : eulerian_du) {
        dU.second.setVal(0.0);
    }

    for (const auto& act : actions) {
        // calculate any contributions to dU (dU += f(old))
        act->calc_time_derivative(this, eulerian_du, time, dt);
        act->calc_spatial_derivative(this, eulerian_du, time, dt);
    }

    // add dU to new data (new = old + dU)
    apply_derivative(eulerian_du);

    //    plot_FAB_2d(eulerian_du[0].second, 0, 0, "eulerian_du", false, true);

    // update new data with any one-shot updates based on old data (new = f(old))
    for (const auto& act : actions) {
        act->apply_time_derivative(this, time, dt);
        act->apply_spatial_derivative(this, time, dt);
    }

}

void MFP::advance_strang(Real time, Real dt, int iteration, int ncycle)
{
    // zero out the derivative accumulator
    for (auto& dU : eulerian_du) {
        dU.second.setVal(0.0);
    }

    for (const auto& act : actions) {
        // calculate any contributions to dU (dU += f(old))
        act->calc_time_derivative(this, eulerian_du, time, dt/2.0);

        // update new data with any one-shot updates based on old data (new* = f(old))
        act->apply_time_derivative(this, time, dt/2.0);
    }

    // add dU to new data (new* = old + dU)
    apply_derivative(eulerian_du);

    for (const auto& act : actions) {
        // calculate any contributions to dU (dU += f(new*))
        act->calc_spatial_derivative(this, eulerian_du, time+dt, dt);

        // update new data with any one-shot updates based on new* data (new** = f(new*))
        act->apply_spatial_derivative(this, time+dt, dt);
    }

    // add dU to new data (new** = new* + dU)
    apply_derivative(eulerian_du);

    for (const auto& act : actions) {
        // calculate any contributions to dU (dU += f(new**))
        act->calc_time_derivative(this, eulerian_du, time+dt, dt/2.0);

        // update new data with any one-shot updates based on old data (new = f(old))
        act->apply_time_derivative(this, time+dt, dt/2.0);

    }

    // add dU to new data (new = new** + dU)
    apply_derivative(eulerian_du);
}
