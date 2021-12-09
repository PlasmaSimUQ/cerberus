#ifdef EILMER_GAS
#include "MFP_gas_kinetics.H"
#include "MFP.H"
#include "MFP_state.H"
#include "MFP_wrap_gas.H"
#include "sol.hpp"
#include <algorithm>

std::string GasKinetics::tag = "gas_kinetics";
bool GasKinetics::registered = GetActionFactory().Register(GasKinetics::tag, ActionBuilder<GasKinetics>);

GasKinetics::GasKinetics(){}
GasKinetics::~GasKinetics(){}

GasKinetics::GasKinetics(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    const sol::table state_names = def["states"];

    for (const auto& key_value_pair : state_names) {
        std::string state_name = key_value_pair.second.as<std::string>();
        State& istate = MFP::get_state(state_name);

        switch (istate.get_type()) {
        case State::StateType::Hydro:
            states.push_back(static_cast<HydroState*>(&istate));
            state_indexes.push_back(istate.global_idx);
            istate.associated_actions.push_back(action_idx);
            break;
        default:
            Abort("An invalid state has been defined for the 'gas' action '"+name+"'");
        }
    }

    num_states = states.size();

    // gas construction
    EilmerGas::initialise();

    std::string gas_model_file = def.get_or<std::string>("gas_model","");
    if (gas_model_file.empty()) Abort("Action '"+name+"' requires 'gas_model' to be defined (a lua file)");

    gas_model_id = gas_model_new(gas_model_file.data());

    if (gas_model_id < 0) Abort("Action '"+name+"' has failed when trying to create a new gas model");

    n_species = gas_model_n_species(gas_model_id);

    species_info.resize(n_species);

    std::pair<bool, int > found;
    int n;
    for (int i=0; i<n_species; ++i) {
        std::string& sname = species_info[i].name;
        sname.resize(10);
        gas_model_species_name_and_length(gas_model_id, i, sname.data(), &n);
        sname.resize(n);

        // search through the neutrals and ions to get where it is located (neutral/ion, idx)
        for (int s=0; s<num_states; ++s) {
            found = findInVector(states[s]->gas->comp_names, sname);
            if (found.first) {
                species_info[i].state_idx = s;
                species_info[i].alpha_idx = found.second;
                break;
            }
        }

        if (!found.first)
            Abort("Couldn't find sub-component for '"+sname+"'");
    }

    gas_state_id = gas_state_new(gas_model_id);

    if (gas_state_id < 0) Abort("Action '"+name+"' has failed when trying to create a new gas state");


    std::string chemistry_update_file = def.get_or<std::string>("chemistry_update","");
    if (chemistry_update_file.empty()) Abort("Action '"+name+"' requires 'chemistry_update' to be defined (a lua file)");

    std::string chemistry_update_file2 = def.get_or<std::string>("chemistry_update2","");

    thermochemical_reactor_id = thermochemical_reactor_new(gas_model_id, chemistry_update_file.data(), chemistry_update_file2.data());

    if (thermochemical_reactor_id < 0) Abort("Action '"+name+"' has failed when trying to create a thermochemical reactor");


    return;
}

void GasKinetics::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("Collisions::get_data");

    Vector<Array<int,2>> options(states.size());

    for (size_t s=0; s<num_states; ++s) {
        options.push_back({states[s]->global_idx, 0});
    }

    Action::get_data(mfp, options, update, time);

}


void GasKinetics::calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt)
{
    BL_PROFILE("Collisions::calc_time_derivative");

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    Vector<Array4<const Real>> U4(num_states);
    Vector<Array4<Real>> dU4(num_states);
    Vector<Vector<Real>> U_old(num_states), U_new(num_states), Q(num_states), alpha(num_states);
    Vector<Vector<Array<Real,+HydroDef::ConsIdx::NUM>>> accounting(num_states);

    for (int si=0; si<num_states; ++si) {
        update[states[si]->data_idx].dU_status = UpdateData::Status::Changed;
        U_old[si].resize(states[si]->n_cons());
        U_new[si].resize(states[si]->n_cons());
        Q[si].resize(states[si]->n_prim());
        alpha[si].resize(states[si]->gas->n_species());
        accounting[si].resize(states[si]->gas->n_species());
    }

    // overall gas state
    Vector<Real> massf(n_species);

    // somewhere to put the contributions to the gas from each state
    Array<Vector<Real>,+HydroDef::ConsIdx::NUM> s_conserved;
    for (int i=0; i<+HydroDef::ConsIdx::NUM; ++i) {
        s_conserved[i].resize(num_states);
    }

    Real dt_suggest;

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);


#ifdef AMREX_USE_EB
        // get the EB data required for later calls and check if we can skip this FAB entirely

        EBData& eb = mfp->get_eb_data(states[0]->global_idx);
        const FArrayBox& vfrac = eb.volfrac[mfi];
        if (vfrac.getType() == FabType::covered) continue;

        Array4<const Real> const& vf4 = vfrac.array();

#endif
        for (int si=0; si<num_states; ++si) {
            U4[si] = update[states[si]->data_idx].U.array(mfi);
            dU4[si] = update[states[si]->data_idx].dU.array(mfi);
        }

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (vf4(i,j,k) == 0.0) continue;
#endif

                    // get the conserved quantities and the alpha values
                    for (int si=0; si<num_states; ++si) {
                        Vector<Real>& UU = U_old[si];
                        const Array4<const Real>& UU4 = U4[si];
                        for (size_t l=0; l<UU.size(); ++l) {
                            UU[l] = UU4(i,j,k,l);
                        }

                        if (UU[+HydroDef::ConsIdx::Density] < states[si]->effective_zero) {
                            std::fill(alpha[si].begin(), alpha[si].end(), 0.0);
                            std::fill(Q[si].begin(), Q[si].end(), 0.0);
                        } else {
                            states[si]->gas->cons2prim(UU, Q[si]);
                            states[si]->gas->get_alpha_fractions_from_prim(Q[si], alpha[si]);
                        }
                    }

                    // get the overall sum of mass, momentum and energy held by the gas state
                    for (int ci=0; ci<+HydroDef::ConsIdx::NUM; ++ci) {
                        std::fill(s_conserved[ci].begin(), s_conserved[ci].end(), 0.0);
                    }
                    for (int isp=0; isp<n_species; ++isp) {
                        const int sid = species_info[isp].state_idx;
                        const int aid = species_info[isp].alpha_idx;
                        for (int ci=0; ci<+HydroDef::ConsIdx::NUM; ++ci) {
                            s_conserved[ci][sid] += alpha[sid][aid]*U_old[sid][ci];
                        }
                    }

                    Array<Real, +HydroDef::ConsIdx::NUM> cons_sum;
                    for (int ci=0; ci<+HydroDef::ConsIdx::NUM; ++ci) {
                        cons_sum[ci] = sum_vec(s_conserved[ci]);
                    }

                    const Real rho_sum = cons_sum[+HydroDef::ConsIdx::Density];
                    const Real mx_sum  = cons_sum[+HydroDef::ConsIdx::Xmom];
                    const Real my_sum  = cons_sum[+HydroDef::ConsIdx::Ymom];
                    const Real mz_sum  = cons_sum[+HydroDef::ConsIdx::Zmom];
                    const Real nrg_sum = cons_sum[+HydroDef::ConsIdx::Eden];

                    const Real int_nrg_sum = nrg_sum - 0.5*(mx_sum*mx_sum + my_sum*my_sum + mz_sum*mz_sum)/rho_sum;
                    const Real specific_internal_nrg = int_nrg_sum/rho_sum;

                    // get the mass fractions
                    for (int isp=0; isp<n_species; ++isp) {
                        const int sid = species_info[isp].state_idx;
                        const int aid = species_info[isp].alpha_idx;
                        massf[isp] = alpha[sid][aid]*U_old[sid][+HydroDef::ConsIdx::Density]/rho_sum;
                    }

                    // calculate the mass weighted temperature
                    Real guess_T = 0.0;
                    for (int si=0; si<num_states; ++si) {
                        guess_T += Q[si][+HydroDef::PrimIdx::Temp]*s_conserved[+HydroDef::ConsIdx::Density][si]/rho_sum;
                    }

                    // set properties - note the conversion to dimensional values
                    gas_state_set_scalar_field(gas_model_id, "rho", rho_sum*MFP::rho_ref);
                    gas_state_set_scalar_field(gas_model_id, "u", specific_internal_nrg*MFP::prs_ref/MFP::rho_ref);
                    gas_state_set_scalar_field(gas_model_id, "T", guess_T*MFP::T_ref); // need an initial guess at the temperature
                    gas_state_set_array_field(gas_state_id, "massf", massf.data(), n_species);

                    // update model and solve for update
                    gas_model_gas_state_update_thermo_from_rhou(gas_model_id, gas_state_id);
                    thermochemical_reactor_gas_state_update(thermochemical_reactor_id, gas_state_id, dt*MFP::t_ref, &dt_suggest);

                    // retrieve mass fraction data - this is all we can get out
                    gas_state_get_array_field(gas_state_id, "massf", massf.data(), n_species);

                    // get the original amount of stuff across all states
                    for (int si=0; si<num_states; ++si) {
                        for (int ai=0; ai<alpha[si].size(); ++ai) {

                            const int ci = +HydroDef::ConsIdx::Density;
                            accounting[si][ai][ci] = U_old[si][ci]*alpha[si][ai]; // original amount
                        }
                    }

                    // get the updated amounts of stuff for all sub-components that have been updated by the reactor
                    // 1. get the change in mass for each species
                    // 2. assume constant temperature and velocity for each grouping and recalculate momentum and energy
                    // 3. shuffle momentum and energy to ensure conservation

                    for (int isp=0; isp<n_species; ++isp) {
                        SpeciesInfo& si = species_info[isp];
                        const int ci = +HydroDef::ConsIdx::Density;
                        accounting[si.state_idx][si.alpha_idx][ci] = massf[isp]*cons_sum[ci]; // updated amount
                    }

                    Array<Real,+HydroDef::ConsIdx::NUM> cons_sums_new, cons_sums_old;
                    cons_sums_new.fill(0.0);
                    cons_sums_old.fill(0.0);
                    for (int si=0; si<num_states; ++si) {

                        const int ci = +HydroDef::ConsIdx::Density;
                        // accumulate how much stuff is held by all the sub-components in this state
                        Real total_density = 0.0;
                        for (const auto& acc : accounting[si]) {
                            total_density += acc[ci];
                        }

                        // apply update for density
                        dU4[si](i,j,k,ci) += total_density - U4[si](i,j,k,ci);

                        // set the density in the primitives vector
                        Q[si][+HydroDef::PrimIdx::Density] = total_density;

                        // apply the update for tracers
                        // note that we can simply use the density of the sub-component which we have already calculated
                        for (int ti=0; ti<states[si]->n_tracers; ++ti) {
                            dU4[si](i,j,k,+HydroDef::ConsIdx::NUM + ti) += accounting[si][ti][ci]  - U4[si](i,j,k,+HydroDef::ConsIdx::NUM + ti);

                            // update the alpha values in the primitives vector
                            Q[si][+HydroDef::PrimIdx::NUM + ti] = accounting[si][ti][ci]/total_density;
                        }

                        // now recalculate the conserved properties
                        states[si]->gas->prim2cons(Q[si],U_new[si]);

                        // calculate the sums of momentum and energy before and after
                        for (int csi=0; csi<+HydroDef::ConsIdx::NUM; ++csi) {
                            cons_sums_new[csi] += U_new[si][csi];
                            cons_sums_old[csi] += U_old[si][csi];
                        }
                    }

                    // calculate and apply the conservative updates
                    for (int csi=+HydroDef::ConsIdx::Xmom; csi<+HydroDef::ConsIdx::NUM; ++csi) {
                        if (cons_sums_new[csi] > 0) {
                            const Real factor = cons_sums_old[csi]/cons_sums_new[csi] ;
                            if (factor != 1.0) {
                                for (int si=0; si<num_states; ++si) {
                                    U_new[si][csi] *= factor;
                                }
                            }
                            for (int si=0; si<num_states; ++si) {
                                dU4[si](i,j,k,csi) += U_new[si][csi] - U4[si](i,j,k,csi);
                            }
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
#endif
