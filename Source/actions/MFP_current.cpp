#include "MFP_current.H"

std::string Current::tag = "current";
bool Current::registered = GetActionFactory().Register(Current::tag, ActionBuilder<Current>);


Current::Current(){}
Current::~Current(){}

Current::Current(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    const std::string state_name = def.get_or<std::string>("state","");

    if (!state_name.empty()) {
        State& istate = MFP::get_state(state_name);
        switch (istate.get_type()) {
        case State::StateType::Field:
            field = static_cast<FieldState*>(&istate);
            state_indexes.push_back(istate.global_idx);
            break;
        default:
            Abort("An invalid state has been defined for the Current source "+name);
        }
    } else {
        Abort("Action '"+name+"' requires 'state' to be defined.");
    }

    scale_factor = -MFP::Larmor/(MFP::lightspeed*MFP::Debye*MFP::Debye);

    // get our current source

    get_udf(def["x"], current[0], 0.0);
    get_udf(def["y"], current[1], 0.0);
    get_udf(def["z"], current[2], 0.0);


    return;
}

void Current::calc_time_derivative(MFP* mfp, Vector<std::pair<int,MultiFab>>& dU, const Real time, const Real dt)
{
    BL_PROFILE("Current::calc_time_derivative");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    // mark dU components that have been touched
    dU[field->data_idx].first = 1;

    std::map<std::string, Real> Q{{"x",0.0}, {"y",0.0}, {"z",0.0}, {"t",time}};
    Real x, y, z;

    const Real* dx = mfp->Geom().CellSize();
    const Real* prob_lo = mfp->Geom().ProbLo();

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(field->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif

        Array4<Real> const& field_dU4 = dU[field->data_idx].second.array(mfi);


        for     (int k = lo.z; k <= hi.z; ++k) {
            z = prob_lo[2] + (k + 0.5)*dx[2];
            Q["z"] = z;
            for   (int j = lo.y; j <= hi.y; ++j) {
                y = prob_lo[1] + (j + 0.5)*dx[1];
                Q["y"] = y;
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {
                    x = prob_lo[0] + (i + 0.5)*dx[0];
                    Q["x"] = x;

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) {
                        continue;
                    }
#endif

                    for (int d = 0; d<3; ++d) {
                        field_dU4(i,j,k,+FieldDef::ConsIdx::Dx + d) += dt*scale_factor*current[d](Q);
                    }
                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
