#include "MFP_hydro_shockdetector.H"

HydroShockDetector::HydroShockDetector() {}

HydroShockDetector::~HydroShockDetector() {}

ClassFactory<HydroShockDetector>& GetHydroShockDetectorFactory()
{
    static ClassFactory<HydroShockDetector> F;
    return F;
}

//================================================================================

#include "MFP_hydro.H"

std::string PressureJumpShockDetector::tag = "pressure_jump_detector";
bool PressureJumpShockDetector::registered = GetHydroShockDetectorFactory().Register(
  PressureJumpShockDetector::tag, HydroShockDetectorBuilder<PressureJumpShockDetector>);

PressureJumpShockDetector::PressureJumpShockDetector() {}
PressureJumpShockDetector::PressureJumpShockDetector(const sol::table& def)
{
    BL_PROFILE("PressureJumpShockDetector::PressureJumpShockDetector");

    idx = def["global_idx"];

    shock_threshold = def["threshold"].get_or(-1.0);

    if (shock_threshold < 0.0) {
        Abort("Error: 'pressure_jump_detector' requires 'threshold' to be set (>0.0)");
    }
}

Real PressureJumpShockDetector::solve(Vector<Real>& L, Vector<Real>& R) const
{
    BL_PROFILE("PressureJumpShockDetector::solve");

    Real pL = L[+HydroDef::PrimIdx::Prs];
    Real pR = R[+HydroDef::PrimIdx::Prs];
    Real varphi = std::abs(pR - pL) / (pR + pL);

    return 0.5 + 0.5 * tanh_approx(5 * (varphi - 0.75 * shock_threshold) / shock_threshold);
}

void PressureJumpShockDetector::write_info(nlohmann::json& js) const
{
    BL_PROFILE("PressureJumpShockDetector::write_info");

    nlohmann::json& sd = js["shock_detector"];

    sd["name"] = tag;
    sd["threshold"] = shock_threshold;
}
