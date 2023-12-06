#include "MFP_lagrangian.H"

#include "MFP.H"
#ifdef AMREX_PARTICLES

LagrangianState::LagrangianState() {}

LagrangianState::~LagrangianState() {}

void LagrangianState::init_from_lua() {}

void LagrangianState::write_info(nlohmann::json& js) const
{
    BL_PROFILE("LagrangianState::write_info");

    State::write_info(js);
}

#endif
