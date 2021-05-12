#include "MFP_shockdetector.H"
#include"MFP_global.H"

using GD = GlobalData;

ShockDetector::ShockDetector(){}

ShockDetector::~ShockDetector(){}

Real ShockDetector::solve(Vector<Real> &L,Vector<Real> &R) const {return 0.0;}

PhysicsFactory<ShockDetector>& GetShockDetectorFactory()
{
    static PhysicsFactory<ShockDetector> F;
    return F;
}

//================================================================================

#include "MFP_hydro.H"

std::string HydroShockDetector::tag = "pressure_jump_detector";
bool HydroShockDetector::registered = GetShockDetectorFactory().Register(HydroShockDetector::tag, ShockDetectorBuilder<HydroShockDetector>);

HydroShockDetector::HydroShockDetector(){}
HydroShockDetector::HydroShockDetector(const sol::table& def)
{
    idx = def["global_idx"];
    istate = GD::get_state_ptr(idx);

    shock_threshold = def["threshold"].get_or(-1.0);

    if (shock_threshold < 0.0) {
        Abort("Error: 'pressure_jump_detector' requires 'threshold' to be set (>0.0)");
    }

}

Real HydroShockDetector::solve(Vector<Real> &L, Vector<Real> &R) const
{
    BL_PROFILE("HydroShockDetector::solve");

    Real pL = L[+HydroState::PrimIdx::Prs];
    Real pR = R[+HydroState::PrimIdx::Prs];
    Real varphi = std::abs(pR - pL)/(pR + pL);

    return 0.5 + 0.5*tanh_approx(5*(varphi - 0.75*shock_threshold)/shock_threshold);
}

bool HydroShockDetector::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s == +StateType::isHydro ) return true;

    return false;
}

void HydroShockDetector::write_info(nlohmann::json& js) const
{
    nlohmann::json& sd = js["shock_detector"];

    sd["name"] = tag;
    sd["threshold"] = shock_threshold;
}
