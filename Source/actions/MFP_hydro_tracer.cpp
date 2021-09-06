#ifdef AMREX_PARTICLES
#include "MFP_hydro_tracer.H"
#include "sol.hpp"
#include "MFP_diagnostics.H"

std::string HydroTracer::tag = "hydro_tracer";
bool HydroTracer::registered = GetActionFactory().Register(HydroTracer::tag, ActionBuilder<HydroTracer>);


HydroTracer::HydroTracer(){}
HydroTracer::~HydroTracer(){}

HydroTracer::HydroTracer(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    std::string hydro_name = def["fluid"];
    std::string particle_name = def["particles"];

    for (const auto& i : MFP::eulerian_states) {
        EulerianState& istate = EulerianState::get_state_global(i);
        if (istate.get_type() == State::StateType::Hydro) {
            if (istate.name == hydro_name) {
                hydro_state = static_cast<HydroState*>(&istate);
            }
        }
    }

    if (!hydro_state) Abort("Source '"+name+"' unable to attach to hydro state '"+hydro_name+"'");

    for (const auto& i : MFP::lagrangian_states) {
        LagrangianState& istate = LagrangianState::get_state_global(i);
        if (istate.get_type() == State::StateType::TracerParticle) {
            if (istate.name == particle_name) {
                tracer_state = static_cast<TracerParticle*>(&istate);
            }
        }
    }

    if (!tracer_state) Abort("Source '"+name+"' unable to attach to tracer particles '"+particle_name+"'");
}

void HydroTracer::apply_spatial_derivative(MFP* mfp, const Real time, const Real dt)
{
    BL_PROFILE("HydroTracer::solve");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    // grab the density and as many momentum components as necessary
    int num_grow = 1;
    MultiFab hydro_cons(mfp->boxArray(), mfp->DistributionMap(), AMREX_SPACEDIM+1, num_grow, MFInfo(),mfp->Factory());
    mfp->FillPatch(*mfp, hydro_cons, num_grow, time, hydro_state->data_idx, +HydroDef::ConsIdx::Density, AMREX_SPACEDIM+1);

    FArrayBox vel;

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();

#ifdef AMREX_USE_EB
        EBData& eb = mfp->get_eb_data(hydro_state->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        const EBCellFlagFab& flag = eb.flags[mfi];
        if (vfrac.getType() == FabType::covered) continue;
#endif

        hydro_state->calc_velocity(grow(box,1),
                                   hydro_cons[mfi],
                                   vel
                           #ifdef AMREX_USE_EB
                                   ,vfrac
                           #endif
                                   );

//        plot_FAB_2d(vel, 0, "x-vel",false,true);

        tracer_state->push_particles(mfp->get_level(),
                                     mfi,
                                     vel,
                                     mfp->Geom(),
                                     dt
                             #ifdef AMREX_USE_EB
                                     ,flag
                             #endif
                                     );




        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}

#endif
