#include "MFP_hydro_riemann.H"

#include "MFP_utility.H"

#include <algorithm>
#include <iostream>
#include <math.h>

HydroRiemannSolver::HydroRiemannSolver() {}

HydroRiemannSolver::~HydroRiemannSolver()
{
    // do nothing
}

ClassFactory<HydroRiemannSolver>& GetHydroRiemannSolverFactory()
{
    static ClassFactory<HydroRiemannSolver> F;
    return F;
}
