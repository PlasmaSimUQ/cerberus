#include "MFP_lagrangian.H"
#include "MFP.H"
#ifdef AMREX_PARTICLES

LagrangianState::LagrangianState(){}

LagrangianState::~LagrangianState(){}

void LagrangianState::init_from_lua()
{
    BL_PROFILE("LagrangianState::init_from_lua");


}




void LagrangianState::write_info(nlohmann::json &js) const
{

    State::write_info(js);

}

#endif
