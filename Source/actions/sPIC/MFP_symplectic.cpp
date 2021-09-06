#ifdef SYMPLECTIC
#include "MFP_symplectic.H"
#include "MFP.H"

std::string Symplectic::tag = "symplectic";
bool Symplectic::registered = GetActionFactory().Register(Symplectic::tag, ActionBuilder<Symplectic>);

Symplectic::Symplectic(){}
Symplectic::~Symplectic(){}

Symplectic::Symplectic(const int idx, const sol::table &def)
{
    action_idx = idx;
    name = def["name"];

    time_order = def.get_or("time_order",2);
    allow_field_generation = def.get_or("allow_field_generation",1);

    const sol::table state_names = def["states"];

    bool has_field = false;

    for (const auto& key_value_pair : state_names) {
        std::string state_name = key_value_pair.second.as<std::string>();

        State& istate = MFP::get_state(state_name);

        switch (istate.get_type()) {
        case State::StateType::Field: {
            if (has_field) Abort("Only one field state can be set for the symplectic action '"+name+"'");
            field = static_cast<FieldState*>(&istate);
            state_indexes.push_back(istate.global_idx);
            field->associated_actions.push_back(action_idx);
            has_field = true;
            break;
        }
        case State::StateType::ChargedParticle: {
            ChargedParticle* cpart = static_cast<ChargedParticle*>(&istate);
            species.push_back(cpart);
            state_indexes.push_back(istate.global_idx);
            cpart->associated_actions.push_back(action_idx);
            break;
        }
        default:
            Abort("An invalid state has been defined for the Symplectic action '"+name+"'");
        }
    }

    return;
}

void Symplectic::apply_change(MFP* mfp, const Real time, const Real dt)
{
    switch(time_order) {
    case 1:
        Theta_B(mfp,time,dt);
        Theta_E(mfp,time,dt);
        Theta_R(mfp,time,dt,Z);
        Theta_R(mfp,time,dt,Y);
        Theta_R(mfp,time,dt,X);
        break;
    case 2:
        Theta_E(mfp,time,dt/2);
        Theta_R(mfp,time,dt/2,X);
        Theta_R(mfp,time,dt/2,Y);
        Theta_R(mfp,time,dt/2,Z);
        Theta_B(mfp,time,dt);
        Theta_R(mfp,time,dt/2,Z);
        Theta_R(mfp,time,dt/2,Y);
        Theta_R(mfp,time,dt/2,X);
        Theta_E(mfp,time,dt/2);
        break;
    default:
        Abort("Time order of "+num2str(time_order)+" in action '"+name+"' is invalid");
    }
}
#endif
