#include "MFP_ode_system.h"

#ifdef USE_CVODE
#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.      */
#include <nvector/nvector_serial.h>         /* access to serial N_Vector                */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix                */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver          */
#include <cvode/cvode_direct.h>           /* access to CVDls interface                */
#include <sundials/sundials_types.h>         /* definition of realtype                   */
#include <sundials/sundials_math.h>          /* contains the macros ABS, SUNSQR, and EXP */
#endif

#include "MFP_source.H"

#include "MFP_global.H"
#include "MFP_ode_solver.h"

using GD = GlobalData;

ODESystem::ODESystem()
{
    global2solver_index.resize(GD::num_solve_state, -1);
    local2global_index.resize(GD::num_solve_state, -1);
}


ODESystem::~ODESystem()
{
    // do nothing
}

void ODESystem::add_source(std::unique_ptr<SourceTerm> src)
{

    src->parent = this;
    // now assign local index values
    // each source that is added is checked for what states it has
    // and a local index assigned in the order of their insertion
    int global_idx;
    const int src_size = src->offsets.size();
    for (int src_idx = 0; src_idx < src_size; ++src_idx) {
        global_idx = src->offsets[src_idx].global;
        State &istate = *GD::states[global_idx];
        if (global2solver_index[global_idx] == -1) { // hasn't been added yet
            global2solver_index[global_idx] = n_components;
            local2global_index[n_states] = global_idx;
            n_components += istate.n_cons();
            n_states += 1;
        }

        src->offsets[src_idx].solver = global2solver_index[global_idx];

    }

    if (src->has_face_src()) {
        has_face_src.push_back(n_sources);
    }

    sources.emplace_back(std::move(src));
    n_sources += 1;

    is_valid.resize(n_states, 1);
}

void ODESystem::set_solver(std::unique_ptr<SolveODE> solve)
{
    solver = std::move(solve);
    solver->parent = this;
}

bool ODESystem::empty()
{
    return sources.empty();
}

int ODESystem::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt)
{
    // initialise to zero
    ydot.resize(n_components);
    std::fill(ydot.begin(), ydot.end(), 0);


    // run through all the source terms
    Vector<Real> yd(n_components);
    for (const auto& src : sources) {

        // make sure every function call gets a clean vector
        std::fill(yd.begin(), yd.end(), 0);

        src->fun_rhs(x, y, z, t, y0, yd, dt);

        for (int i=0; i<ydot.size(); ++i) {
            ydot[i] += yd[i];
        }
    }

    return 0;
}

int ODESystem::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J)
{

    // initialise
    J.resize(n_components*n_components);
    std::fill(J.begin(), J.end(), 0);

    // run through all the source terms
    Vector<Real> JJ(n_components*n_components);
    for (const auto& src : sources) {

        // make sure each function call gets a clean matrix
        std::fill(JJ.begin(), JJ.end(), 0);

        src->fun_jac(x, y, z, t, y0, JJ);

        for (int i=0; i<JJ.size(); ++i) {
            J[i] += JJ[i];
        }

    }

    return 0;
}

#ifdef USE_CVODE
int ODESystem::fun_rhs_cvode(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    UserData* data = (UserData*) user_data;
    ODESystem* s = data->sys;

    // zero out ydot
    for (int i=0; i < s->n_components; ++i) {
        NV_Ith_S(ydot, i) = 0.0;
    }

    // copy of y
    Vector<Real> yy(NV_DATA_S(y), NV_DATA_S(y)+NV_LENGTH_S(y));

    // run through all the source terms
    Vector<Real> yd(s->n_components);
    for (const auto& src : s->sources) {

        // make sure every function call gets a clean vector
        std::fill(yd.begin(), yd.end(), 0);
        src->fun_rhs(data->x, data->y, data->z, t, yy, yd);

        for (int i=0; i < s->n_components; ++i) {
            NV_Ith_S(ydot, i) += yd[i];
        }
    }

    return 0;
}

int ODESystem::fun_jac_cvode(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
                               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{

    UserData* data = (UserData*) user_data;
    ODESystem* s = data ->sys;

    int n_terms = s->n_components;

    //zero out
    realtype* Jdata = SUNDenseMatrix_Data(J);

    // copy of y
    Vector<Real> yy(NV_DATA_S(y), NV_DATA_S(y)+NV_LENGTH_S(y));

    // run through all the source terms
    Vector<Real> JJ(n_terms*n_terms);
    for (const auto& src : s->sources) {

        std::fill(JJ.begin(), JJ.end(), 0);

        src->fun_jac(data->x, data->y, data->z, tn, yy, JJ);

        for (int i=0; i<JJ.size(); ++i) {
            Jdata[i] += JJ[i];
        }

    }

    return 0;
}
#endif


const int ODESystem::get_solver_state_valid(const int local_idx) const
{
    return is_valid.at(local_idx);
}

void ODESystem::set_solver_state_valid(const int local_idx, const int v)
{
    is_valid[local_idx] = v;

    for (auto &src : sources) {
        // need to check if this solver index is actually held by the source!!!
        src->set_local_idx_validity(local_idx, (bool)v);
    }
}

std::string ODESystem::print()
{
    std::string msg;

    msg =  "* ODE system '" + name + "' active:\n";
    msg += "    solver type : "+solver->get_tag()+"\n";
    msg += "    source terms: (type : name)\n";

    for (const auto& src : sources) {
        msg += "      " + src->print() + "\n";
    }

    msg += "    states: (type : name)\n";
    for (const int &global_idx : local2global_index) {
        if (global_idx < 0) break;
        State &istate = GD::get_state(global_idx);
        int tp = istate.get_type();

        msg += "      " + istate.get_tag() + " : " + istate.name + "\n";
    }

    return msg;
}
