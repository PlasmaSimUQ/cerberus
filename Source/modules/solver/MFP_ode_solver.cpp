#include "MFP_ode_solver.h"

#include<AMReX_Print.H>
#include <Dense>

#include "MFP_utility.H"
#include "MFP_global.H"
#include "MFP_source.H"

using GD = GlobalData;

// ====================================================================================

SolveODE::SolveODE(){}
SolveODE::~SolveODE(){}

int SolveODE::init_data(Vector<Real> *y,
                        Vector<Real> *yout)
{
    y0 = y;
    y1 = yout;

    return 0;
}

bool SolveODE::valid_solution()
{
    // check for validity of state vector y1
    bool valid = true;
    int of = 0;
    int global_idx, nc;
    for (int solver_idx = 0; solver_idx < parent->n_states; ++solver_idx) {
        global_idx = parent->local2global_index[solver_idx];
        State &istate = *GD::states[global_idx];
        nc = istate.n_cons();

        if (parent->get_solver_state_valid(solver_idx)) {

            // get a copy of the species state
            Vector<Real> state(&y1->at(of), &y1->at(of)+nc);


            // check for validity
            valid = istate.cons_valid(state);
            if (!valid) return valid;
        }

        of += nc;
    }
    return valid;
}

int SolveODE::solve(Real x, Real y, Real z, Real t0, Real t1, int depth=0){return 0;}
void SolveODE::clear(){}

PhysicsFactory<SolveODE>& GetSolveODEFactory()
{
    static PhysicsFactory<SolveODE> F;
    return F;
}
