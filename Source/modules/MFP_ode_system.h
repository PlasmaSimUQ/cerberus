#ifndef ODESYSTEM_H
#define ODESYSTEM_H

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#ifdef USE_CVODE
#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.      */
#include <nvector/nvector_serial.h>         /* access to serial N_Vector                */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix                */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver          */
#include <cvode/cvode_direct.h>           /* access to CVDls interface                */
#include <sundials/sundials_types.h>         /* definition of realtype                   */
#include <sundials/sundials_math.h>          /* contains the macros ABS, SUNSQR, and EXP */
#endif

#include "MFP_ode_solver.h"

using namespace amrex;

// forward declaration
class SourceTerm;
class ODESystem;

typedef struct UserData {Real x; Real y; Real z; ODESystem* sys;} UserData;

class ODESystem
{
public:
    ODESystem();
    ~ODESystem();

    void add_source(std::unique_ptr<SourceTerm> src);
    SourceTerm* get_source(int i) const {return sources.at(i).get();}
    void set_solver(std::unique_ptr<SolveODE> solve);
    bool empty();

    int fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt=0);
    int fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J);

    #ifdef USE_CVODE
    static int fun_rhs_cvode(realtype t, N_Vector y, N_Vector ydot, void *user_data);
    static int fun_jac_cvode(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    #endif

    const int get_solver_state_valid(const int local_idx) const;
    void set_solver_state_valid(const int local_idx, const int v);

    std::string print();

    Vector<std::unique_ptr<SourceTerm> > sources;

    // Each solver has n states, with each state having m components,
    // The y vector is serialized:
    //     y = [S00, S01, ..., S0m, ..., Sn0, Sn1, ..., Snm]
    // where the first index is the i-th state {i: 0..n} and
    // the second index is the j-th component {j: 0..m}.
    //
    // global_idx = solver2global_index[i] will give us the global index
    // for the i-th state of the solver.
    //
    // pos = global2solver_index[global_idx] will give us the position of the 0th
    // component for that state in the y-vector. In this case,
    // global_idx corresponds with the i-th state in the solver,
    //     y[pos] = Si0
    //     y[pos + 1] = Si1,
    // and so on.

    // note that we may have index vectors of size greater than there are
    // solvers in the system. It is up to the developer to check for a valid
    // index when using these index lists
    Vector<int> global2solver_index;
    Vector<int> local2global_index;
    Vector<int> has_face_src;
    Vector<int> is_valid;

    std::unique_ptr<SolveODE> solver;

    int n_sources = 0;
    int n_states = 0;
    int n_components = 0;
    std::string name;
};

#endif // ODESYSTEM_H
