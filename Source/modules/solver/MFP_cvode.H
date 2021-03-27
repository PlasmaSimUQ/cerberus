#ifdef USE_CVODE
#ifndef CVODE_H
#define CVODE_H

#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.      */
#include <nvector/nvector_serial.h>         /* access to serial N_Vector                */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix                */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver          */
#include <cvode/cvode_direct.h>           /* access to CVDls interface                */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sundials/sundials_types.h>         /* definition of realtype                   */
#include <sundials/sundials_math.h>          /* contains the macros ABS, SUNSQR, and EXP */

#include "MFP_ode_solver.h"
#include "MFP_ode_system.h"

class SolveCVODE : public SolveODE
{
public:
  SolveCVODE();
  SolveCVODE(const sol::table& def);
  ~SolveCVODE();
  void clear() override;
  virtual int init_data(Vector<Real> *y, Vector<Real> *yout) override;
  virtual int solve(Real x, Real y, Real z, Real t0, Real t1, int depth=0) override;
  virtual SolveType get_type() const override {return SolveType::CVODE;}
  virtual std::string get_tag() const override {return tag;}

  // custom definition of valid_solution() for cvode data types
  virtual bool valid_solution() override;
  static std::string tag;
  static bool registered;

  void* cvode_mem;
  N_Vector y0;
  N_Vector y1;
  SUNMatrix Amat;
  SUNLinearSolver Linsol;
  SUNNonlinearSolver NLS;
  UserData user_data;

  Real atol = 1e-14;
  Real rtol = 1e-6;
  int max_step = 200;
  int jac = 1;
  int apply_constraints = 0;

  int cv_lmm = 2;
  int cv_nls_fixed_point = 0;

  int n_sub = 2;
  int depth;
  int max_depth = 4;
  int verbosity = 1;


};

#endif // CVODE_H
#endif
