#ifdef USE_CVODE
#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.      */
#include <nvector/nvector_serial.h>         /* access to serial N_Vector                */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix                */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver          */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_types.h>         /* definition of realtype                   */
#include <sundials/sundials_math.h>          /* contains the macros ABS, SUNSQR, and EXP */

#include "MFP_cvode.H"

#include "MFP_utility.H"
#include "MFP_global.H"
#include "MFP_source.H"
#include "MFP_ode_system.h"

using GD = GlobalData;

std::string SolveCVODE::tag = "CVODE";
bool SolveCVODE::registered = GetSolveODEFactory().Register(SolveCVODE::tag, SolveODEBuilder<SolveCVODE>);


SolveCVODE::SolveCVODE(){}
SolveCVODE::SolveCVODE(const sol::table& def)
{
    sol::optional<sol::table> options = def["options"];
    if (options) {
        sol::table &opt = options.value();
        n_sub = opt["refine_factor"].get_or(2);
        max_depth = opt["refine_limit"].get_or(4);
        atol = opt["abstol"].get_or(1e-14);
        rtol = opt["reltol"].get_or(1e-6);
        max_step = opt["nsteps"].get_or(200);
        jac = opt["jacobian"].get_or(1);
        verbosity = opt["verbosity"].get_or(1);

        apply_constraints = opt["CVODE_apply_constraints"].get_or(0);
        cv_lmm = opt["CVODE_method"].get_or(2);
        cv_nls_fixed_point = opt["CVODE_nls_fixed_point"].get_or(0);
  }

  depth = 0;

}

SolveCVODE::~SolveCVODE(){}

int SolveCVODE::init_data(Vector<Real> *y, Vector<Real> *yout)
{
    BL_PROFILE_VAR("CVODE::init_data()", cvid);

    int N = y->size();

    // CVODE set up
    y0 = N_VMake_Serial(N, y->data());
    y1 = N_VMake_Serial(N, yout->data());

    if (cv_lmm == 1) {
        cvode_mem = CVodeCreate(CV_ADAMS);
    } else if (cv_lmm == 2) {
        cvode_mem = CVodeCreate(CV_BDF);
    }

    int flag;

    flag = CVodeInit(cvode_mem, parent->fun_rhs_cvode, 0.0, y0);
    flag = CVodeSStolerances(cvode_mem, rtol, atol);

    if (cv_nls_fixed_point == 1) {
        NLS = SUNNonlinSol_FixedPoint(y0, 0);
        flag = CVodeSetNonlinearSolver(cvode_mem, NLS);
    }

    // TODO: selection of the linear system could be made an option
    Amat = SUNDenseMatrix(N, N);
    Linsol = SUNDenseLinearSolver(y0, Amat);
    flag = CVodeSetLinearSolver(cvode_mem, Linsol, Amat);

    // use the user_defined fun_jac if needed
    if (jac) {
        flag = CVodeSetJacFn(cvode_mem, parent->fun_jac_cvode);
    }

    flag = CVodeSetMaxNumSteps(cvode_mem, max_step);

    // apply physical constraints if specified.
    // This should only be used if setting atol to
    // a small value doesn't work.
    // See CVODE user guide section 4.5.2, "Advice on controlling unphysical negative values"
    if (apply_constraints == 1) {
        N_Vector constraints = NULL;
        constraints = N_VNew_Serial(N);
        N_VConst_Serial(0.0, constraints);

        int global_idx, state_pos_in_y;
        for (int i = 0; i < parent->n_states; ++i) {

            // get the global index of the i-th state in the solver
            global_idx = parent->local2global_index[i];
            State &istate = *GD::states[global_idx];

            if (istate.get_enforce_positivity()) {
                state_pos_in_y = parent->global2solver_index[global_idx];

                //  0.0 -> no constraint
                //  1.0 -> Density or Eden >= 0.0
                // -1.0 -> Density or Eden <= 0.0
                //  2.0 -> Density or Eden >  0.0
                // -2.0 -> Density or Eden <  0.0
                for (const int cidx : istate.get_enforce_positivity_idx()) {

                    NV_Ith_S(constraints, state_pos_in_y + cidx) = 2.0;

                }
            }
        }

        flag = CVodeSetConstraints(cvode_mem, constraints);
    }

    BL_PROFILE_VAR_STOP(cvid);

    return flag;
}

int SolveCVODE::solve(Real x, Real y, Real z, Real t0, Real t1, int depth)
{
    BL_PROFILE_VAR("CVODE::solve()", cvs);

    int flag;
    realtype tout;

    user_data.x = x;
    user_data.y = y;
    user_data.z = z;
    user_data.sys = parent;

    flag = CVodeSetUserData(cvode_mem, &user_data);
    flag = CVodeReInit(cvode_mem, t0, y0);

    flag = CVode(cvode_mem, t1, y1, &tout, CV_NORMAL);

    if ((flag != CV_SUCCESS) || !valid_solution()) {

        if (depth == max_depth) {
            return 1;
        }

        depth += 1;

        if (verbosity >= 2) {
            amrex::Print() << "CVODE time refinement level : /\\ " << depth;
            if (flag != CV_SUCCESS) {
                amrex::Print() << " (cvode error " << flag << ")" << std::endl;
            } else {
                amrex::Print() << " (conservation error)" << std::endl;
            }
        }

        realtype dt = (t1 - t0)/n_sub;

        realtype t0_sub = t0;
        realtype t1_sub;
        for (int i=0; i<n_sub; ++i) {
            t0_sub += i*dt;
            t1_sub = t0_sub + dt;

            if (solve(x, y, z, t0_sub, t1_sub, depth) != CV_SUCCESS) return 1;

            N_VScale_Serial(1.0, y1, y0);

        }

        depth -= 1;

        if (verbosity >= 2)
            amrex::Print() << "CVODE time refinement level : \\/ " << depth << std::endl;

    }

    BL_PROFILE_VAR_STOP(cvs);
    return 0;
}

void SolveCVODE::clear()
{
    N_VDestroy(y0);
    N_VDestroy(y1);
    CVodeFree(&cvode_mem);
    SUNMatDestroy(Amat);
    SUNLinSolFree(Linsol);
}

bool SolveCVODE::valid_solution()
{
    BL_PROFILE("SolveCVODE::valid_solution");
    // check for validity of state vector y1
    bool valid = true;
    int of = 0;
    int global_idx, nc;
    int i = 0;

    for (int solver_idx = 0; solver_idx < parent->n_states; ++solver_idx) {

        global_idx = parent->local2global_index[solver_idx];
        State &istate = *GD::states[global_idx];

        nc = istate.n_cons();

        // get a copy of the species state
        Vector<Real> state(&NV_Ith_S(y1, of), &NV_Ith_S(y1, of)+nc);
        of += nc;

        // check for validity
        valid = istate.cons_valid(state);
        if (!valid) return valid;
    }
    return valid;
}

#endif
