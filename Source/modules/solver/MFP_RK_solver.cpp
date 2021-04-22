#include "MFP_RK_solver.H"

#include "MFP_utility.H"
#include "MFP_global.H"
#include "MFP_source.H"

#include "Eigen"

using GD = GlobalData;

// ====================================================================================

SolveRK::SolveRK(){}
SolveRK::SolveRK(const sol::table& def)
{
    // get any solver options

    sol::optional<sol::table> options = def["options"];
    if (options) {
        sol::table &opt = options.value();
        n_sub = opt["refine_factor"].get_or(2);
        max_depth = opt["refine_limit"].get_or(4);
        verbosity = opt["verbosity"].get_or(0);
        order = opt["order"].get_or(0);
    }

    if ((order < 1) || (order > 3))
        Abort("Order of calculation in explicit solver must be in the range 1 to 3");
}
SolveRK::~SolveRK(){}


int SolveRK::solve(Real x, Real y, Real z, Real t0, Real t1, int depth)
{
    BL_PROFILE("SolveRK::solve");
    Real dt = t1 - t0;
    int n_terms = parent->n_components;

    // first order
    euler_step_solve(x, y, z, t0, dt, *y0, *y1);

    // RK2 - 2 stage - CFL<=1
    if (order == 2) {

        Vector<Real> f1(n_terms);
        euler_step_solve(x, y, z, t0+dt, dt, *y1, f1);

        for (int i=0; i<n_terms; ++i) {
            (*y1)[i] = 0.5*(*y0)[i] + 0.5*f1[i];
        }
    }

    // RK3 - 4 stage - CFL<=2
    if (order == 3) {
        Vector<Real> f1(n_terms);
        for (int i=0; i<n_terms; ++i) {
            f1[i] = 0.5*(*y0)[i] + 0.5*(*y1)[i];
        }

        Vector<Real> f2(n_terms);
        euler_step_solve(x, y, z, t0+0.5*dt, dt, f1, f2);
        for (int i=0; i<n_terms; ++i) {
            f2[i] = 0.5*f1[i] + 0.5*f2[i];
        }

        Vector<Real> f3(n_terms);
        euler_step_solve(x, y, z, t0+dt, dt, f2, f3);
        for (int i=0; i<n_terms; ++i) {
            f3[i] = (2.0/3.0)*(*y0)[i] + (1.0/6.0)*f2[i] + (1.0/6.0)*f3[i];
        }

        euler_step_solve(x, y, z, t0+0.5*dt, dt, f3, *y1);
        for (int i=0; i<n_terms; ++i) {
            (*y1)[i] = 0.5*f3[i] + 0.5*(*y1)[i];
        }
    }


    // check for validity of state vector and recursively reduce the time step
    // according to the passed in parameters
    if (!valid_solution()) {

        if (depth == max_depth) {
            amrex::Abort("Source integration time reduction went too far");
            return 1;
        }

        depth += 1;

        if (verbosity >= 2) {
            amrex::Print() << "ODE time refinement level : /\\ " << depth;
            amrex::Print() << " (conservation error)" << std::endl;
        }


        dt = (t1 - t0)/n_sub;

        Real t0_sub = t0;
        Real t1_sub;
        for (int i=0; i<n_sub; ++i) {
            t0_sub += i*dt;
            t1_sub = t0_sub + dt;

            if (solve(x, y, z, t0_sub, t1_sub, depth)) return 1;

            std::copy(y1->begin(), y1->end(), y0->begin());

        }

        depth -= 1;

        if (verbosity >= 2)
            amrex::Print() << "ODE time refinement level : \\/ " << depth << std::endl;

    }

    return 0;

}


// ====================================================================================

std::string SolveExplicit::tag = "explicit";
bool SolveExplicit::registered = GetSolveODEFactory().Register(SolveExplicit::tag, SolveODEBuilder<SolveExplicit>);


SolveExplicit::SolveExplicit(){}
SolveExplicit::SolveExplicit(const sol::table& def) : SolveRK(def){}
SolveExplicit::~SolveExplicit(){}

void SolveExplicit::euler_step_solve(Real x, Real y, Real z, Real t0, Real dt, Vector<Real>& y0, Vector<Real>& y1)
{
    BL_PROFILE("SolveExplicit::euler_step_solve");
    int n_terms = parent->n_components;

    parent->fun_rhs(x, y, z, t0, y0, y1, dt);
    for (int i=0; i<n_terms; ++i) {
        y1[i] = y0[i] + dt*y1[i];
    }

    return;

}

// ====================================================================================

Eigen::ColPivHouseholderQR<Eigen::MatrixXd> SolveImplicit::linear_solve;

std::string SolveImplicit::tag = "implicit";
bool SolveImplicit::registered = GetSolveODEFactory().Register(SolveImplicit::tag, SolveODEBuilder<SolveImplicit>);


SolveImplicit::SolveImplicit(){}
SolveImplicit::SolveImplicit(const sol::table& def) : SolveRK(def){}
SolveImplicit::~SolveImplicit(){}

void SolveImplicit::euler_step_solve(Real x, Real y, Real z, Real t0, Real dt, Vector<Real>& y0, Vector<Real>& y1)
{
    BL_PROFILE("SolveImplicit::euler_step_solve");

    int n_terms = parent->n_components;

    // calculate the jacobian
    Vector<Real> J;
    parent->fun_jac(x, y, z, t0, y0, J);

    // get an alias to the Jacobian
    Eigen::Map<Eigen::MatrixXd> A(J.data(), n_terms, n_terms);

    // I - dt*J
    A *= -dt;
    A.diagonal().array() += 1;

    // get an alias to the source terms vector and output vector
    Eigen::Map<Eigen::VectorXd> B(y0.data(), n_terms);
    Eigen::Map<Eigen::VectorXd> X(y1.data(), n_terms);

    X.setZero();

    //    if (GD::msg) {
    //        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
    //        std::cout << "A = " << A.format(HeavyFmt) << std::endl;
    //        std::cout << "b = " << B.format(HeavyFmt) << std::endl;
    //        std::cout << "#---------------------------------" << std::endl;
    //    }


    /* In order to reduce the dimensionality of the linear system as
     * much as possible we only solve a linear system where the
     * corresponding components of the jacobian (J) have non-zero
     * columns (within some tolerance).
     *
     * J[:,i] == 0 : nothing else depends on the i'th component
     *  if J[i,:] != 0 then compute change based on the updated values,
     *  if J[i,:] == 0 then copy it from the b vector
     *
     * J[i,:] != 0 and J[:,i] != 0 : It needs to be in the linear system
     *      as it is both updated and is required for other quantities
     *      to be updated.
     */

    // get a list of the entries to keep, remove, or post-calculate
    Vector<int> remove, keep, post;

    // get the maximum value of any entry in the rows and columns
    Eigen::VectorXd max_col = Eigen::VectorXd::Zero(n_terms);
    Eigen::VectorXd max_row = Eigen::VectorXd::Zero(n_terms);
    Real abs_val;
    for (int i=0; i<n_terms; ++i) {
        for (int j=0; j<n_terms; ++j) {
            if (i == j)
                continue;

            abs_val = std::abs(A(i,j));

            if (abs_val > max_row(i))
                max_row(i) = std::abs(A(i,j));

            if (abs_val > max_col(j))
                max_col(j) = std::abs(A(i,j));
        }
    }

    const Real tol = 1e-15;
    bool rm, ps;
    for (int i=0; i<n_terms; ++i) {
        rm = false;
        ps = false;

        // nothing elese depends on this element, remove it
        if (max_col(i) <= tol) {
            rm = true;
            // does it need updating post-solve?
            if (max_row(i) > tol) {
                ps = true;
            }
        }

        // check the diagonal
        if ((std::abs(A(i,i) - 1.0) > tol) && rm) {
            ps = true;
        }

        if (rm) {
            remove.push_back(i);
        } else {
            keep.push_back(i);
        }

        if (ps) {
            post.push_back(i);
        }
    }

    //    if (GD::msg) {
    //        std::cout << "remove = " << vec2str(remove) << std::endl;
    //        std::cout << "keep   = " << vec2str(keep) << std::endl;
    //        std::cout << "post   = " << vec2str(post) << std::endl;
    //        std::cout << "#---------------------------------" << std::endl;
    //    }

    int n = keep.size();

    if (n > 0) {

        // make a new linear system of reduced rank
        Eigen::MatrixXd A1(n, n);
        Eigen::VectorXd b1(n);
        Eigen::VectorXd x1(n);

        // get the subset
        A1 = A(keep, keep);
        b1 = B(keep);


        //        if (GD::msg) {
        //            Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
        //            std::cout << "A1 = " << A1.format(HeavyFmt) << std::endl;
        //            std::cout << "b1 = " << b1.format(HeavyFmt) << std::endl;
        //            std::cout << "dt = " << dt << std::endl;
        //            std::cout << "#---------------------------------" << std::endl;
        //        }

        // solve the linearised system for the updated state vector at time t1
        linear_solve.compute(A1);
        x1 = linear_solve.solve(b1);

        //        if (GD::msg) {
        //            Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
        //            std::cout << "x1 = " << x1.format(HeavyFmt) << std::endl;
        //            std::cout << "#---------------------------------" << std::endl;
        //        }

        // insert the reduced rank solution back into the overall solution
        for (int i=0; i<n; ++i) {
            X[keep[i]] = x1[i];
        }
    }

    // copy the items that would have been untouched by the source terms
    for (const int &i : remove) {
        X[i] = B[i];
    }

    if (!post.empty()) {
        // calculate the source terms from the updated state
        Vector<Real> ydot(n_terms);
        parent->fun_rhs(x, y, z, t0, y1, ydot, dt);

        // apply post-solve updates
        for (const int &i : post) {
            X[i] += dt*ydot[i];
        }
    }

    //    if (GD::msg) {
    //        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
    //        std::cout << "X = " << X.format(HeavyFmt) << std::endl;
    //        std::cout << "#---------------------------------" << std::endl;
    //    }

    return;


}
