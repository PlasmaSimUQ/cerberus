#include "MFP_udf.H"

std::string UserDefined::tag = "user";
bool UserDefined::registered = GetActionFactory().Register(UserDefined::tag, ActionBuilder<UserDefined>);


UserDefined::UserDefined(){}
UserDefined::~UserDefined(){}

UserDefined::UserDefined(const int idx, const sol::table &def)
{

    action_idx = idx;
    name = def["name"];

    const std::string state_name = def.get_or<std::string>("state","");

    if (!state_name.empty()) {
        state = &EulerianState::get_state(state_name);
        state_indexes.push_back(state->global_idx);
    } else {
        Abort("Action '"+name+"' requires 'state' to be defined.");
    }

    // get all of our user defined functions
    terms.resize(state->get_num_cons());
    const Vector<std::string>& cons_names = state->get_cons_names();
    for (int j = 0; j<cons_names.size(); ++j) {
        get_udf(def["value"][cons_names[j]], terms[j], 0.0);
    }

    return;
}

void UserDefined::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("UserDefined::get_data");

    Vector<Array<int,2>> options = {{state->global_idx, 0}};

    Action::get_data(mfp, options, update, time);

}

void UserDefined::calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt)
{
    BL_PROFILE("UserDefined::calc_time_derivative");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    update[state->data_idx].dU_status = UpdateData::Status::Changed;

    std::map<std::string, Real> Q{{"x",0.0}, {"y",0.0}, {"z",0.0}, {"t",time}};
    Real x, y, z;

    const Real* dx = mfp->Geom().CellSize();
    const Real* prob_lo = mfp->Geom().ProbLo();

    const size_t n_cons = terms.size();

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(state->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif

        Array4<Real> const& dU4 = update[state->data_idx].dU.array(mfi);


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

                    for (int d = 0; d<n_cons; ++d) {
                        const Real udf = dt*terms[d](Q);
                        dU4(i,j,k,d) += udf;
                    }
                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
