#include "MFP_damp_div_mhd.H"

std::string DampDivergenceMHD::tag = "damp_divergence_mhd";
bool DampDivergenceMHD::registered = GetActionFactory().Register(DampDivergenceMHD::tag, ActionBuilder<DampDivergenceMHD>);

DampDivergenceMHD::DampDivergenceMHD(){}
DampDivergenceMHD::~DampDivergenceMHD(){}

DampDivergenceMHD::DampDivergenceMHD(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    div_damp = def.get_or("div_damp",0.0);

    if (div_damp <= 0.0) {
        Abort("Why does action '"+name+"' of type 'damp_divergence_mhd' action have damp <= 0.0?");
    }

    const std::string state_name = def.get_or<std::string>("state","");

    if (!state_name.empty()) {
        State& istate = MFP::get_state(state_name);
        switch (istate.get_type()) {
        case State::StateType::MHD:
            mhd = static_cast<MHDState*>(&istate);
            state_indexes.push_back(istate.global_idx);
            break;
        default:
            Abort("An invalid state has been defined for the DampDivergenceMHD source "+name);
        }
    } else {
        Abort("Action '"+name+"' requires 'state' to be defined.");
    }

    div_transport = mhd->div_transport;

    if (div_transport <= 0.0) {
        Abort("action '"+name+"' is associated with a state that has div_transport <= 0.0");
    }

    return;
}

void DampDivergenceMHD::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("DampDivergenceMHD::get_data");

    Vector<Array<int,2>> options = {{mhd->global_idx, 0}};

    Action::get_data(mfp, options, update, time);

}

void DampDivergenceMHD::calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt)
{
    BL_PROFILE("DampDivergenceMHD::calc_time_derivative");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    // mark dU components that have been touched
    update[mhd->data_idx].dU_status = UpdateData::Status::Changed;

    Vector<Real> U(+MHDDef::ConsIdx::NUM);

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(mhd->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif

        Array4<Real> const& mhd4 = update[mhd->data_idx].U.array(mfi);
        Array4<Real> const& mhd_dU4 = update[mhd->data_idx].dU.array(mfi);


        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) {
                        continue;
                    }
#endif
                    // use the local speed to calculate the damping

                    for (size_t n=0; n<+MHDDef::ConsIdx::NUM; ++n) {
                        U[n] = mhd4(i,j,k,n);
                    }

                    const RealArray s = mhd->get_speed_from_cons(U);

                    Real ch = 0.0;

                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        ch = std::max(ch, std::abs(s[d]));
                    }


                    ch *= div_transport;

                    if (ch > 0.0) {
                        const Real cd = ch/div_damp;
                        mhd_dU4(i,j,k,+MHDDef::ConsIdx::psi) -= dt*cd*mhd4(i,j,k,+MHDDef::ConsIdx::psi);

                    }
                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
