#include "MFP_elastic.H"
#include "MFP_global.H"

using GD = GlobalData;

std::string Elastic::tag = "elastic";
bool Elastic::registered = GetSourceTermFactory().Register(Elastic::tag, SourceTermBuilder<Elastic>);


Elastic::Elastic(){}

Elastic::Elastic(const sol::table& def)
{
    name = def.get<std::string>("name");

    if (!Elastic::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &Elastic::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }

    if (index.size() < 2) {
        Abort("Error: For elastic source `"+name+"`, you must define only two states: a + b -> a + b");
    }

    std::string el_mom_fit_file = def["el_mom_fit_file"];
    if (el_mom_fit_file.length() == 0) {
        Abort("Error: el_mom_fit_file must be defined for the elastic source term: `"+tag+"`");
    }

    std::string el_nrg_fit_file = def["el_nrg_fit_file"];
    if (el_mom_fit_file.length() == 0) {
        Abort("Error: el_nrg_fit_file must be defined for the elastic source term: `"+tag+"`");
    }
    mom_interp = Interp2D(el_mom_fit_file);
    nrg_interp = Interp2D(el_nrg_fit_file);

    return;
}

Elastic::~Elastic()
{
    // do nothing
}

Vector<Real> Elastic::elastic(const Vector<Real> &y0, const Vector<OffsetIndex> &offsets, const Interp2D &mom_interp, const Interp2D &nrg_interp)
{
    BL_PROFILE("Elastic::elastic");
    Vector<Real> ydot(y0.size());

    const Real T_ref = GD::T_ref;
    const Real x_ref = GD::x_ref;
    const Real n_ref = GD::n_ref;
    const Real u_ref = GD::u_ref;

    int num_hydro = offsets.size();

    // vector for hydro primitive values
    Vector<Vector<Real>> hydro_prim(offsets.size());

    int nc;
    for (const auto &idx : offsets) {
        State &istate = GD::get_state(idx.global);

        // get a copy of the conserved variables
        Vector<Real> &U = hydro_prim[idx.local];
        U.resize(istate.n_cons());
        for (int i = 0; i < U.size(); ++i) {
        U[i] = y0[idx.solver + i];
        }

        // convert to primitive
        istate.cons2prim(U);
    }
    // grab properties for state a and b
    const OffsetIndex &offset_a = offsets[0];
    State &state_a = GD::get_state(offset_a.global);
    if (!offset_a.valid) Abort("State '"+state_a.name+"' is unavailable for source of type '"+tag);
    Vector<Real> &Q_a = hydro_prim[offset_a.local];

    const OffsetIndex &offset_b = offsets[1];
    State &state_b = GD::get_state(offset_b.global);
    if (!offset_b.valid) Abort("State '"+state_b.name+"' is unavailable for source of type '"+tag);
    Vector<Real> &Q_b = hydro_prim[offset_b.local];

    Real rho_a = Q_a[+HydroState::PrimIdx::Density];
    Real alpha_a = Q_a[+HydroState::PrimIdx::Alpha];
    Real m_a = state_a.get_mass(alpha_a);
    Real p_a = Q_a[+HydroState::PrimIdx::Prs];
    Real n_a = rho_a / m_a;
    Real T_a = p_a / n_a;

    Real rho_b = Q_b[+HydroState::PrimIdx::Density];
    Real alpha_b = Q_b[+HydroState::PrimIdx::Alpha];
    Real m_b = state_b.get_mass(alpha_b);
    Real p_b = Q_b[+HydroState::PrimIdx::Prs];
    Real n_b = rho_b / m_b;
    Real T_b = p_b / n_b;

    // Total mass
    Real M = m_a + m_b;

    Array<Real, 3> du;
    Array<Real, 3> U_eq;
    Real du2;
    Real du_dot_U_eq;

    // temp
    Real ua;
    Real ub;
    Real Ueq;
    Real delu;
    for (int i = 0; i < 3; ++i) {
        ua = Q_a[+HydroState::PrimIdx::Xvel + i];
        ub = Q_b[+HydroState::PrimIdx::Xvel + i];

        delu = ua - ub;
        du[i] = delu;
        du2 += delu * delu;

        // COM velocity
        Ueq = (m_a / M) * ua + (m_b / M) * ub;
        U_eq[i] = Ueq;

        du_dot_U_eq += Ueq * delu;
    }


    // COM mass
    Real mu = (m_a * m_b) / M;

    // COM temperature
    Real T_eq = (m_a*T_b + m_b*T_a) / M;
    Real T_eq_real = T_eq * T_ref;

    // multi-fluid factor
    Real lam = 0.5 * mu * du2 / T_eq;

    // non-dimensional momentum rate
    Real mom_rate = mu * n_a * n_b * mom_interp.interpolate(T_eq_real, lam) * x_ref / u_ref * n_ref;
    // Real mom_rate = 0;

    // non_dimensional energy rate
    Real nrg_rate = n_a * n_b * nrg_interp.interpolate(T_eq_real, lam) * x_ref / u_ref * n_ref;

    for (int i = 0; i < 3; ++i) {
        ydot[offset_a.solver + +HydroState::ConsIdx::Xmom + i] += -du[i] * mom_rate;
        ydot[offset_b.solver + +HydroState::ConsIdx::Xmom + i] += du[i] * mom_rate;
    }

    ydot[offset_a.solver + +HydroState::ConsIdx::Eden] += -du_dot_U_eq * mom_rate - 2 * m_a / (m_a + m_b) * (T_a - T_b) * nrg_rate;
    ydot[offset_b.solver + +HydroState::ConsIdx::Eden] += du_dot_U_eq * mom_rate + 2 * m_a / (m_a + m_b) * (T_a - T_b) * nrg_rate;

    return ydot;
}

int Elastic::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{
    BL_PROFILE("Elastic::fun_rhs");
    const int n_terms = y0.size();

    // copy to eigen vector
    Vector<Real> yy(n_terms);
    for (int i=0; i<n_terms; ++i) {
        yy[i] = y0[i];
    }

    // call source function
    ydot = elastic(yy, offsets, mom_interp, nrg_interp);

    return 0;
}

int Elastic::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{
    BL_PROFILE("Elastic::fun_jac");
    const int n_terms = y0.size();

    Vector<Real> ydot(n_terms);
    num_jac(x, y, z, t, y0, ydot, J);

    return 0;
}

bool Elastic::valid_state(const int global_idx)
{

    State &istate = GD::get_state(global_idx);

    switch (istate.get_type()) {
    case +StateType::isHydro:
        return true;
    case +StateType::isHydro2P:
        return true;
    default:
        return false;
    }
}

bool Elastic::valid_solver(const int solve_idx)
{
    return true;
}
