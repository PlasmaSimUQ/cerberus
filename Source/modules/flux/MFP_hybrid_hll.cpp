#include "MFP_hybrid_hll.H"
#include "MFP_utility.H"
#include "MFP_global.H"

using GD = GlobalData;

//================================================================================

std::string HydroHybridHLL::tag = "HLLE/HLLC";
bool HydroHybridHLL::registered = GetRiemannSolverFactory().Register(HydroHybridHLL::tag, RiemannSolverBuilder<HydroHybridHLL>);

HydroHybridHLL::HydroHybridHLL() {}
HydroHybridHLL::HydroHybridHLL(const int i)
{
    idx = i;
    hllc = HydroHLLC(i);
    hlle = HydroHLLE(i);
    istate = GD::get_state_ptr(i);
    F_hlle.resize(+HydroState::ConsIdx::NUM);
    F_hllc.resize(+HydroState::ConsIdx::NUM);
}

void HydroHybridHLL::solve(Vector<Real> &L,
                           Vector<Real> &R,
                           Vector<Real> &F,
                           Real* shk)
{
    BL_PROFILE("HydroHybridHLL::solve");
    const Real eps = 1e-14;

    if (*shk < eps) {
        hllc.solve(L, R, F, shk);
    } else if ((1.0 - *shk) < eps) {
        hlle.solve(L, R, F, shk);
    } else {
        hllc.solve(L, R, F_hllc, shk);
        hlle.solve(L, R, F_hlle, shk);
        for (int i=0; i<+HydroState::ConsIdx::NUM; ++i) {
            F[i] = (1.0 - *shk)*F_hllc[i] + *shk*F_hlle[i];
        }
    }
}

bool HydroHybridHLL::valid_state(const int idx)
{
    int s = GD::get_state(idx).get_type();

    if (s != +StateType::isHydro) {
        return false;
    }
    return true;
}
