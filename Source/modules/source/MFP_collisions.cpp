#include "MFP_collisions.H"
#include "MFP_global.H"

using GD = GlobalData;

std::string Collisions::tag = "collisions";
bool Collisions::registered = GetSourceTermFactory().Register(Collisions::tag, SourceTermBuilder<Collisions>);

Collisions::Collisions(){}

Collisions::Collisions(const sol::table& def)
{
    name = def.get<std::string>("name");

    if (!Collisions::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &Collisions::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }

    if (index.size() < 2)
        Abort("Source '"+name+" needs at least 2 states to collide");

    // grab any collision cross sections
    sol::table cross_sections = def["cross_sections"].get_or(sol::table());

    if (cross_sections.empty()) {

        // only need ccs with a neutral species
        bool has_neutral = false;
        for (const int& idx : index) {
            State& istate = GD::get_state(idx);
            if ((istate.charge[0] == 0) || (istate.charge[1] == 0)) {
                has_neutral = true;
            }
        }
        if (has_neutral)
            Abort("Error: collisions enabled with a neutral species without hydro_ccs being defined");
    }

    ccs.resize(index.size());
    for (int a=0; a<index.size(); ++a) {
        ccs[a].resize(index.size(),0);
    }

    // load the hydro cross sections
    std::pair<bool, int> index_a, index_b;
    Real sigma;
    for (const auto& a : cross_sections) {
        index_a = findInVector(includes, a.first.as<std::string>());
        if (!index_a.first) {continue;}
        for (const auto& b : a.second.as<sol::table>()) {
            index_b = findInVector(includes, b.first.as<std::string>());
            if (!index_b.first) {continue;}
            sigma = b.second.as<Real>()/(GD::x_ref*GD::x_ref); // non-dimensionalise
            ccs[index_a.second][index_b.second] = sigma;
            ccs[index_b.second][index_a.second] = sigma; // cross load as sigma_ab = sigma_ba
        }
    }


    return;
}

Collisions::~Collisions()
{
    // do nothing
}

Vector<dual> Collisions::collisions(const Vector<dual> &y0,
                                    const Vector<OffsetIndex> &offsets,
                                    const Vector<Vector<Real>> &ccs)
{

    Real Debye = GD::Debye;
    Real n0 = GD::n0;

    // coefficients
    const Real c1 = 0.02116454531141366/(n0*pow(Debye,4)); // sqrt(2)/(12*pi^(3/2))*...
    const Real c2 = 0.4135669939329334; // (2/(9*pi))**(1/3.0)
    const Real c3 = 2.127692162140974*n0; // (4/3)*n0*sqrt(8/pi)
    const Real lnC = 10.0; // Coulomb logarithm

    int num_hydro = offsets.size();

    // vector for hydro primitive values
    Vector<Vector<dual>> hydro_prim(offsets.size());
    for (const auto &idx : offsets) {

        State &istate = GD::get_state(idx.global);

        // get a copy of the conserved variables
        Vector<dual> &U = hydro_prim[idx.local];
        U.resize(istate.n_cons());
        for (int i=0; i<U.size(); ++i) {
            U[i] = y0[idx.solver+i];
        }

        // convert to primitive
        istate.cons2prim(U);
    }

    // define our output and set it to zero
    Vector<dual> ydot(y0.size());

    for (int a = 0; a < num_hydro; ++a) {

        const OffsetIndex &offset_a = offsets[a];

        State &state_a = GD::get_state(offset_a.global);

        Vector<dual> &Q_a = hydro_prim[offset_a.local];

        dual rho_a = state_a.get_density_from_prim(Q_a);

        if (rho_a.val <= GD::effective_zero || !offset_a.valid) {
            continue;
        }

        int vec_idx_a = state_a.get_cons_vector_idx()[0];
        int nrg_idx_a = state_a.get_nrg_idx()[0];
        Real alpha_a = state_a.get_alpha_from_prim(Q_a).val;

        Real m_a = state_a.get_mass(alpha_a);
        Real q_a = state_a.get_charge(alpha_a);
        Real q_a2 = q_a*q_a;

        dual T_a = state_a.get_temperature_from_prim(Q_a);

        for (int b = a+1; b < num_hydro; ++b) {

            const OffsetIndex &offset_b = offsets[b];

            State &state_b = GD::get_state(offset_b.global);

            Vector<dual> &Q_b = hydro_prim[offset_b.local];

            dual rho_b = state_b.get_density_from_prim(Q_b);

            if (rho_b.val <= GD::effective_zero || !offset_a.valid) {
                continue;
            }

            int vec_idx_b = state_b.get_cons_vector_idx()[0];
            int nrg_idx_b = state_b.get_nrg_idx()[0];
            Real alpha_b = state_b.get_alpha_from_prim(Q_b).val;

            Real m_b = state_b.get_mass(alpha_b);
            Real q_b = state_b.get_charge(alpha_b);
            Real q_b2 = q_b*q_b;

            dual T_b = state_b.get_temperature_from_prim(Q_b);
            dual n_b = state_b.get_density_from_prim(Q_b)/m_b;

            Real m_ab = (m_a*m_b)/(m_a + m_b);


            dual du2 = 0.0;
            Array<dual,3> du;
            for (int i=0; i<3; ++i) {
                du[i] = Q_b[vec_idx_b + i] - Q_a[vec_idx_a + i];
                du2 += du[i]*du[i];
            }

            //collision frequency
            dual nu;
            if ((q_a != 0) && (q_b != 0)) {
                Real coeff_1 = c1*((q_a2*q_b2*lnC)/(m_ab*m_a));
                nu = n_b*coeff_1*pow(c2*du2 + T_a/m_a + T_b/m_b, -1.5);
            } else {
                Real coeff_2 = (m_b*ccs[a][b])/(m_a + m_b);
                nu = n_b*coeff_2*c3*sqrt(T_a/m_a + T_b/m_b);
            }

            // effect of collisions

            dual Q = m_ab*rho_a/m_a*nu*du2 + 3*rho_a*nu/(m_a + m_b)*(T_b - T_a);

            Array<dual,3> R;
            for (int i=0; i<3; ++i) {
                R[i] = rho_a*nu*du[i];
                Q += R[i]*Q_a[vec_idx_a + i];
            }

            for (int i=0; i<3; ++i) {
                ydot[offset_a.solver + vec_idx_a + i] += R[i];
                ydot[offset_b.solver + vec_idx_b + i] -= R[i];
            }

            ydot[offset_a.solver + nrg_idx_a] += Q;
            ydot[offset_b.solver + nrg_idx_b] -= Q;

        }
    }

    return ydot;
}


int Collisions::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{

    const int n_terms = y0.size();

    // copy to eigen vector
    Vector<dual> yy(n_terms);
    for (int i=0; i<n_terms; ++i) {
        yy[i] = y0[i];
    }

    // call source function
    Vector<dual> yd = collisions(yy, offsets, ccs);

    // copy to ydot
    for (int i=0; i<n_terms; ++i) {
        ydot[i] = yd[i].val;
    }

    return 0;
}

int Collisions::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{

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

    // calculate the Jacobian matrix dF/dy
    JJ = jacobian(Collisions::collisions, autodiff::forward::wrt(yy), autodiff::forward::at(yy, offsets, ccs));

    // The estimate of the slope is calculated from the initial state plus the fluxes so
    // it is more J(t+dt/2) than J(t+dt). We thus apply a factor of a half to account for
    // the fact we are applying more of a mid-point rule rathern than backwards Euler.
    JJ *= 0.5;

    return 0;

}

bool Collisions::valid_state(const int global_idx)
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

bool Collisions::valid_solver(const int solve_idx)
{
    return true;
}

