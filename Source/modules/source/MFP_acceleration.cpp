#include "MFP_acceleration.H"

#include "eigen.hpp"
#include "forward.hpp"

using GD = GlobalData;

//---------------------------------------------------------------------------------------------

std::string Acceleration::tag = "acceleration";
bool Acceleration::registered = GetSourceTermFactory().Register(Acceleration::tag, SourceTermBuilder<Acceleration>);

Acceleration::Acceleration(){}

Acceleration::Acceleration(const sol::table &def)
{

    name = def.get<std::string>("name");

    if (!Acceleration::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &Acceleration::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }

    // get all of our manufactured solutions

    AMREX_D_TERM(get_udf(def["x"], acc[0], 0.0);,
                 get_udf(def["y"], acc[1], 0.0);,
                 get_udf(def["z"], acc[2], 0.0);)

    return;
}

Acceleration::~Acceleration()
{
    // do nothing
}

Vector<Real> Acceleration::get_acc(Real x, Real y, Real z, Real t) const
{
    BL_PROFILE("Acceleration::get_acc");
    // calculat the acceleration at this point in time and space
    Vector<Real> a(AMREX_SPACEDIM);

    std::map<std::string, Real> Q{{"x",x}, {"y",y}, {"z",z}, {"t",t}};

    for (int i = 0; i<AMREX_SPACEDIM; ++i) {
        a[i] = acc[i](Q);
    }

    return a;
}

Vector<dual> Acceleration::accelerate(const Vector<dual> &y0,
                                    const Vector<OffsetIndex> &offsets,
                                    const Vector<Real> &a)
{
    BL_PROFILE("Acceleration::accelerate");
    Vector<dual> ydot(y0.size());

    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        State &istate = GD::get_state(idx.global);

        const int di = istate.get_cons_density_idx();
        const int mxi = istate.get_cons_vector_idx()[0];
        const int ei = istate.get_nrg_idx()[0];

        // momentum and energy
        for (int i = 0; i<AMREX_SPACEDIM; ++i) {
            ydot[idx.solver + mxi + i] = a[i]*y0[idx.solver + di]; // g*rho
            ydot[idx.solver + ei]     += a[i]*y0[idx.solver + mxi + i]; // g*rho*u
        }

    }


    return ydot;
}

int Acceleration::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{

    BL_PROFILE("Acceleration::fun_rhs");
    Vector<Real> a = get_acc(x,y,z,t);

    const int n_terms = y0.size();

    // copy to eigen vector
    Vector<dual> yy(n_terms);
    for (int i=0; i<n_terms; ++i) {
        yy[i] = y0[i];
    }

    // call source function
    Vector<dual> yd = accelerate(yy, offsets, a);

    // copy to ydot
    for (int i=0; i<n_terms; ++i) {
        ydot[i] = yd[i].val;
    }

    return 0;
}

int Acceleration::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{
    BL_PROFILE("Acceleration::fun_jac");

//            Vector<Real> ydot;
//            num_jac(x, y, z, t, y0, ydot, J);
//            return 0;

    const int n_terms = y0.size();

    // make an alias
    J.resize(n_terms*n_terms);
    Eigen::Map<Eigen::MatrixXd> JJ(J.data(), n_terms, n_terms);

    // copy to dual vector
    Vector<dual> yy(n_terms);
    for (int i=0; i<n_terms; ++i) {
        yy[i] = y0[i];
    }

    Vector<Real> a = get_acc(x,y,z,t);

    // calculate the Jacobian matrix dF/dy
    JJ = jacobian(Acceleration::accelerate, autodiff::forward::wrt(yy), autodiff::forward::at(yy, offsets, a));

    // The estimate of the slope is calculated from the initial state plus the fluxes so
    // it is more J(t+dt/2) than J(t+dt). We thus apply a factor of a half to account for
    // the fact we are applying more of a mid-point rule rathern than backwards Euler.
    JJ *= 0.5;

    return 0;

}


bool Acceleration::valid_state(const int global_idx)
{
    State &istate = GD::get_state(global_idx);

    switch (istate.get_type()) {
    case +StateType::isMHD:
        return true;
    case +StateType::isHydro:
        return true;
    case +StateType::isHydro2P:
        return true;
    default:
        return false;
    }
}

bool Acceleration::valid_solver(const int solve_idx)
{
    if (solve_idx != +SolveType::Explicit) {
        return false;
    }
    return true;
}
