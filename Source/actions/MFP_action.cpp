#include "MFP_action.H"
#include "MFP.H"
#include "MFP_eulerian.H"

Action::Action(){}
Action::~Action(){}

ClassFactory<Action>& GetActionFactory() {
    static ClassFactory<Action> F;
    return F;
}


void Action::get_data(MFP* mfp, Vector<Array<int,2>>& options, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("Action::get_U0");

#ifdef AMREX_USE_EB
    constexpr int num_grow_eb = 2;
#else
    constexpr int num_grow_eb = 0;
#endif

    // copy the species data
    for (const auto& index : options) {
        const EulerianState& state = EulerianState::get_state_global(index[0]);
        const int data_idx = state.data_idx;
        bool fill = index[1];

        UpdateData::Status U_status = update[data_idx].U_status;
        UpdateData::Status dU_status = update[data_idx].dU_status;

        int ns = state.n_cons();

        ////////////
        /// U

        if (fill) {

            if (U_status != UpdateData::Status::Filled) {

                int ng = state.get_num_grow() + num_grow_eb;

                if ((U_status == UpdateData::Status::Inactive) || (update[data_idx].U.nGrow() < ng)) {
                    update[data_idx].U.define(mfp->boxArray(), mfp->DistributionMap(),
                                              ns,ng,MFInfo(),mfp->Factory());
                }

#ifdef AMREX_USE_EB
                EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(state.eb2_index));
#endif
                mfp->FillPatch(*mfp, update[data_idx].U, ng, time, data_idx, 0, ns);
                update[data_idx].U_status = UpdateData::Status::Filled;



            }

        } else {

            if ((U_status != UpdateData::Status::Local) && (U_status != UpdateData::Status::Filled)) {
                MultiFab& species_data_ref = mfp->get_data(data_idx,time);

                int ng = 0;

                if (U_status == UpdateData::Status::Inactive) {
                    update[data_idx].U.define(mfp->boxArray(), mfp->DistributionMap(),
                                              ns,ng,MFInfo(),mfp->Factory());
                }

                MultiFab::Copy(update[data_idx].U, species_data_ref, 0, 0, ns, ng);
                update[data_idx].U_status = UpdateData::Status::Local;
            }

        }

        /////////
        /// dU


        if (dU_status == UpdateData::Status::Inactive) {
            update[data_idx].dU.define(mfp->boxArray(), mfp->DistributionMap(),
                                       ns,num_grow_eb,
                                       MFInfo(),mfp->Factory());
        }

        if ((dU_status == UpdateData::Status::Expired) || (dU_status == UpdateData::Status::Inactive)) {
            update[data_idx].dU.setVal(0.0);
            update[data_idx].dU_status = UpdateData::Status::Zero;
        }
    }
}
