#include "MFP_acceleration.H"

#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"

std::string Acceleration::tag = "acceleration";
bool Acceleration::registered =
  GetActionFactory().Register(Acceleration::tag, ActionBuilder<Acceleration>);

Acceleration::Acceleration() {}
Acceleration::~Acceleration() {}

Acceleration::Acceleration(const int idx, const sol::table& def)
{
    BL_PROFILE("Acceleration::Acceleration");

    action_idx = idx;
    name = def["name"];

    sol::table acc_vector = def["vector"];
    for (const auto& key_value_pair : acc_vector) {
        int idx = key_value_pair.first.as<int>();
        Real val = key_value_pair.second.as<Real>();
        acc[idx - 1] = val;
    }

    const sol::table state_names = def["states"];

    for (const auto& key_value_pair : state_names) {
        std::string state_name = key_value_pair.second.as<std::string>();
        State& istate = MFP::get_state(state_name);

        switch (istate.get_type()) {
        case State::StateType::Hydro:
            species.push_back(static_cast<HydroState*>(&istate));
            state_indexes.push_back(istate.global_idx);
            break;
        default: Abort("An invalid state has been defined for the Acceleration source " + name);
        }
    }
    return;
}

void Acceleration::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("Acceleration::get_data");

    Vector<Array<int, 2>> options(species.size());

    for (size_t i = 0; i < species.size(); ++i) { options[i] = {species[i]->global_idx, 0}; }

    Action::get_data(mfp, options, update, time);
}

void Acceleration::calc_time_derivative(MFP* mfp,
                                        Vector<UpdateData>& update,
                                        const Real time,
                                        const Real dt)
{
    BL_PROFILE("Acceleration::calc_time_derivative");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    size_t n_species = species.size();

    for (const HydroState* hstate : species) {
        // mark dU components that have been touched
        update[hstate->data_idx].dU_status = UpdateData::Status::Changed;
    }

    Vector<Array4<Real>> species4(n_species);
    Vector<Array4<Real>> species_dU4(n_species);

    // define some 'registers'

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {
        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(species[0]->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif

        for (int n = 0; n < n_species; ++n) {
            species4[n] = update[species[n]->data_idx].U.array(mfi);
            species_dU4[n] = update[species[n]->data_idx].dU.array(mfi);
        }

        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
#ifdef AMREX_USE_EB
                    if (vf4(i, j, k) == 0.0) { continue; }
#endif

                    for (size_t n = 0; n < n_species; ++n) {
                        const Real rho = species4[n](i, j, k, +HydroDef::ConsIdx::Density);

                        // momentum and energy
                        for (int d = 0; d < 3; ++d) {
                            const Real m = species4[n](i, j, k, +HydroDef::ConsIdx::Xmom + d);
                            species_dU4[n](i, j, k, +HydroDef::ConsIdx::Xmom + d) +=
                              dt * acc[d] * rho;  // g*rho
                            species_dU4[n](i, j, k, +HydroDef::ConsIdx::Eden) +=
                              dt * acc[d] * m;  // g*rho*u
                        }
                    }
                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}
