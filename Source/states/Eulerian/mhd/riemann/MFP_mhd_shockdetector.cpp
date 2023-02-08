#include "MFP_mhd_shockdetector.H"

MHDShockDetector::MHDShockDetector(){}

MHDShockDetector::~MHDShockDetector(){}

ClassFactory<MHDShockDetector>& GetMHDShockDetectorFactory()
{
    static ClassFactory<MHDShockDetector> F;
    return F;
}

//================================================================================

#include "MFP_mhd.H"

std::string PressureJumpShockDetectorMHD::tag = "pressure_jump_detector";
bool PressureJumpShockDetectorMHD::registered = GetMHDShockDetectorFactory().Register(PressureJumpShockDetectorMHD::tag, MHDShockDetectorBuilder<PressureJumpShockDetectorMHD>);

PressureJumpShockDetectorMHD::PressureJumpShockDetectorMHD(){}
PressureJumpShockDetectorMHD::PressureJumpShockDetectorMHD(const sol::table& def)
{
    idx = def["global_idx"];

    shock_threshold = def["threshold"].get_or(-1.0);

    if (shock_threshold < 0.0) {
        Abort("Error: 'pressure_jump_detector' requires 'threshold' to be set (>0.0)");
    }

}

Real PressureJumpShockDetectorMHD::solve(Vector<Real> &L,
                                      Vector<Real> &R) const
{
    BL_PROFILE("PressureJumpShockDetector::solve");

    Real pL = L[+MHDDef::PrimIdx::Prs];
    Real pR = R[+MHDDef::PrimIdx::Prs];
    Real varphi = std::abs(pR - pL)/(pR + pL);

    return 0.5 + 0.5*tanh_approx(5*(varphi - 0.75*shock_threshold)/shock_threshold);
}

void PressureJumpShockDetectorMHD::write_info(nlohmann::json& js) const
{
    nlohmann::json& sd = js["shock_detector"];

    sd["name"] = tag;
    sd["threshold"] = shock_threshold;
}
