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
    hlle.idx = i;
    hllc.idx = i;
}

void HydroHybridHLL::solve(Vector<Real> &L,
                           Vector<Real> &R,
                           Vector<Real> &F,
                           Real* shk) const
{
    BL_PROFILE("HydroHybridHLL::solve");
    const Real eps = 1e-14;

    if (*shk < eps) {
        hllc.solve(L, R, F, shk);
    } else if ((1.0 - *shk) < eps) {
        hlle.solve(L, R, F, shk);
    } else {
        Vector<Real> F_hlle(+HydroState::ConsIdx::NUM), F_hllc(+HydroState::ConsIdx::NUM);
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
