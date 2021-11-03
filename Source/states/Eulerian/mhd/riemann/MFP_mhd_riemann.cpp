#include "MFP_mhd_riemann.H"
#include <math.h>
#include <algorithm>
#include <iostream>

#include "MFP_utility.H"

MHDRiemannSolver::MHDRiemannSolver()
{
}

MHDRiemannSolver::~MHDRiemannSolver()
{
    // do nothing
}

ClassFactory<MHDRiemannSolver>& GetMHDRiemannSolverFactory()
{
    static ClassFactory<MHDRiemannSolver> F;
    return F;
}
