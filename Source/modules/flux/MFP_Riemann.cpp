#include "MFP_Riemann.H"
#include <math.h>
#include <algorithm>
#include <iostream>

#include "MFP_utility.H"
#include "MFP_global.H"

using GD = GlobalData;

RiemannSolver::RiemannSolver()
{
}

RiemannSolver::~RiemannSolver()
{
    // do nothing
}

void RiemannSolver::solve(Vector<Real> &L,
                               Vector<Real> &R,
                               Vector<Real> &F,
                               Real* shk) const
{
    // do nothing
}

PhysicsFactory<RiemannSolver>& GetRiemannSolverFactory()
{
    static PhysicsFactory<RiemannSolver> F;
    return F;
}
