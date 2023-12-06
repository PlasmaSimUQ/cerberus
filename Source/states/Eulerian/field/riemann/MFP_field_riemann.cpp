#include "MFP_field_riemann.H"

#include "MFP_utility.H"

#include <algorithm>
#include <iostream>
#include <math.h>

FieldRiemannSolver::FieldRiemannSolver() {}

FieldRiemannSolver::~FieldRiemannSolver()
{
    // do nothing
}

ClassFactory<FieldRiemannSolver>& GetFieldRiemannSolverFactory()
{
    static ClassFactory<FieldRiemannSolver> F;
    return F;
}
