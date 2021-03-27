#include "MFP_damp_div.H"

#include "MFP_global.H"

using GD = GlobalData;

std::string DampDivergenceCorrection::tag = "damp_divergence";
bool DampDivergenceCorrection::registered = GetSourceTermFactory().Register(DampDivergenceCorrection::tag, SourceTermBuilder<DampDivergenceCorrection>);

DampDivergenceCorrection::DampDivergenceCorrection(){}

DampDivergenceCorrection::DampDivergenceCorrection(const sol::table& def)
{

    name = def.get<std::string>("name");

    if (!DampDivergenceCorrection::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &DampDivergenceCorrection::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }

    // get the damping ratio
    damping = 0;
    damping = def["damping_ratio"].get_or(0.0);

    // if the damping won't achieve anything then don't apply it
    bool cancel = false;
    for (const int& global_idx : index) {
        State &istate = GD::get_state(global_idx);

        if ((istate.div_speed) <= 0.0) {
            cancel = true;
            break;
        }
    }

    return;
}

DampDivergenceCorrection::~DampDivergenceCorrection()
{
    // do nothing
}

int DampDivergenceCorrection::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{

    int i;
    Real ch, cd;

    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        State &istate = GD::get_state(idx.global);

        ch = istate.div_speed;
        if (damping > 0.0 && ch > 0.0) {
            cd = ch/damping;

            i = istate.get_cons_psi_idx();
            if (i > -1) {
                ydot[idx.solver + i] -= cd*y0[idx.solver + i];
            }

            i = istate.get_cons_phi_idx();
            if (i > -1) {
                ydot[idx.solver + i] -= cd*y0[idx.solver + i];
            }
        }
    }

    return 0;
}

int DampDivergenceCorrection::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{

    //  num_jac(t, y, ydot, J);

    const int n_terms = y0.size();

    // make an alias
    J.resize(n_terms*n_terms);
    Eigen::Map<Eigen::MatrixXd> JJ(J.data(), n_terms, n_terms);

    Real ch, cd;

    // analytical jacobian
    int i;
    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        State &istate = GD::get_state(idx.global);

        ch = istate.div_speed;
        if (damping > 0.0 && ch > 0.0) {
            cd = ch/damping;

            i = istate.get_cons_psi_idx();
            if (i > -1) {
                i += idx.solver;
                JJ(i,i) = -cd;
            }

            i = istate.get_cons_phi_idx();
            if (i > -1) {
                i += idx.solver;
                JJ(i,i) = -cd;
            }
        }

    }

    return 0;
}

bool DampDivergenceCorrection::valid_state(const int global_idx)
{

    State &istate = GD::get_state(global_idx);

    switch (istate.get_type()) {
    case +StateType::isField:
        return true;
    case +StateType::isMHD:
        return true;
    default:
        return false;
    }
}

bool DampDivergenceCorrection::valid_solver(const int solver_idx)
{
    return true;
}

void DampDivergenceCorrection::write_info(nlohmann::json &js) const
{

    js["damp"] = damping;

}

