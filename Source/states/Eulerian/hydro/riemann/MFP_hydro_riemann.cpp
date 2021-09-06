#include "MFP_hydro_riemann.H"
#include <math.h>
#include <algorithm>
#include <iostream>

#include "MFP_utility.H"

HydroRiemannSolver::HydroRiemannSolver()
{
}

HydroRiemannSolver::~HydroRiemannSolver()
{
    // do nothing
}

ClassFactory<HydroRiemannSolver>& GetHydroRiemannSolverFactory()
{
    static ClassFactory<HydroRiemannSolver> F;
    return F;
}
