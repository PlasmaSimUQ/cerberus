#include "MFP_mhd_riemann.H"

#include "MFP_utility.H"

#include <algorithm>
#include <iostream>
#include <math.h>

MHDRiemannSolver::MHDRiemannSolver() {}

MHDRiemannSolver::~MHDRiemannSolver()
{
    // do nothing
}

ClassFactory<MHDRiemannSolver>& GetMHDRiemannSolverFactory()
{
    static ClassFactory<MHDRiemannSolver> F;
    return F;
}
