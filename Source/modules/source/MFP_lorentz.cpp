#include "MFP_lorentz.H"
#include "MFP_global.H"

using GD = GlobalData;

std::string Lorentz::tag = "Lorentz";
bool Lorentz::registered = GetSourceTermFactory().Register(Lorentz::tag, SourceTermBuilder<Lorentz>);

Lorentz::Lorentz(){}

Lorentz::Lorentz(const sol::table &def)
{

    if (!GD::plasma_params_set) {
        Abort("Plasma parameters required for plasma5 source, define (>0) either skin_depth & beta or Larmor & Debye");
    }

    name = def.get<std::string>("name");

    if (!Lorentz::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &Lorentz::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }

    //------------------------Global state indexes

    linked_EM_idx = -1;

    bool found_em = false;

    for (const auto &idx : offsets) {

        State &istate = GD::get_state(idx.global);

        int tp = istate.get_type();
        switch(tp) {
            case +StateType::isMHD:
                continue;
            case +StateType::isField:
                if (found_em) Abort("Too many fields linked to plasma5 source");
                linked_EM_idx=idx.global;
                found_em = true;
                continue;
            case +StateType::isHydro:
            case +StateType::isHydro2P:
                if (istate.charge[0] < 0.) {
                    linked_electron_idx.push_back(idx.global);
                    continue;
                }
                else {
                    linked_ion_idx.push_back(idx.global);
                    continue;
                }
        }
    }

    if (linked_EM_idx < 0 )
        Abort("No electromagnetic field found for plasma5 source");

}

Lorentz::~Lorentz()
{
    // do nothing
}


int Lorentz::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{

    ydot.resize(y0.size());
    std::fill(ydot.begin(), ydot.end(),0.0);

    // get any magnetic field

    Real pD, pB;
    Real D_clean;
    Real Bx, By, Bz;
    Real Ex, Ey, Ez;
    Real ep;

    int field_offset;

    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        State &istate = GD::get_state(idx.global);
        int tp = istate.get_type();

        if (tp != +StateType::isField)
            continue;

        D_clean = istate.div_speed;

        field_offset = idx.solver;

        ep = y0[field_offset + +FieldState::ConsIdx::ep];

        // magnetic field
        Bx = y0[field_offset + +FieldState::ConsIdx::Bx];
        By = y0[field_offset + +FieldState::ConsIdx::By];
        Bz = y0[field_offset + +FieldState::ConsIdx::Bz];
        pB = y0[field_offset + +FieldState::ConsIdx::psi];

        // electric field
        Ex = y0[field_offset + +FieldState::ConsIdx::Dx]/ep;
        Ey = y0[field_offset + +FieldState::ConsIdx::Dy]/ep;
        Ez = y0[field_offset + +FieldState::ConsIdx::Dz]/ep;
        pD = y0[field_offset + +FieldState::ConsIdx::phi];
    }


    Real q, m, r;
    Real rho, mx, my, mz, alpha;
    Real u, v, w;

    Real Larmor = GD::Larmor;
    Real lightspeed = GD::lightspeed;

    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        State &istate = GD::get_state(idx.global);
        int tp = istate.get_type();

        if (tp == +StateType::isField)
            continue;

        rho =   y0[idx.solver + +HydroState::ConsIdx::Density];
        mx =    y0[idx.solver + +HydroState::ConsIdx::Xmom];
        my =    y0[idx.solver + +HydroState::ConsIdx::Ymom];
        mz =    y0[idx.solver + +HydroState::ConsIdx::Zmom];
        alpha = y0[idx.solver + +HydroState::ConsIdx::Tracer]/rho;

        m = istate.get_mass(alpha);
        q = istate.get_charge(alpha);

        r = q/m;

        u = mx/rho;
        v = my/rho;
        w = mz/rho;

        ydot[idx.solver + +HydroState::ConsIdx::Xmom] += (rho*r/Larmor)*(lightspeed*Ex + v*Bz - w*By);
        ydot[idx.solver + +HydroState::ConsIdx::Ymom] += (rho*r/Larmor)*(lightspeed*Ey + w*Bx - u*Bz);
        ydot[idx.solver + +HydroState::ConsIdx::Zmom] += (rho*r/Larmor)*(lightspeed*Ez + u*By - v*Bx);
    }

    return 0;
}

int Lorentz::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{

    //    Vector<Real> ydot;
    //    num_jac(x, y, z, t, y0, ydot, J);
    //    return 0;

    const int n_terms = y0.size();

    // make an alias for easier access
    J.resize(n_terms*n_terms);
    Eigen::Map<Eigen::MatrixXd> JJ(J.data(), n_terms, n_terms);

    JJ.setZero();

    //

    Real dL = GD::Larmor;

    // get all of the field data (should only be one field associated with this source
    Real ep;
    Real Bx, By, Bz;
    Real Ex, Ey, Ez;
    Real D_clean;

    int field_offset;

    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        State &istate = GD::get_state(idx.global);
        int tp = istate.get_type();

        if (tp != +StateType::isField)
            continue;

        field_offset = idx.solver;

        ep = y0[field_offset + +FieldState::ConsIdx::ep];

        // magnetic field
        Bx = y0[field_offset + +FieldState::ConsIdx::Bx];
        By = y0[field_offset + +FieldState::ConsIdx::By];
        Bz = y0[field_offset + +FieldState::ConsIdx::Bz];

        // electric field (convert from D to E)
        Ex = y0[field_offset + +FieldState::ConsIdx::Dx]/ep;
        Ey = y0[field_offset + +FieldState::ConsIdx::Dy]/ep;
        Ez = y0[field_offset + +FieldState::ConsIdx::Dz]/ep;

        D_clean = istate.div_speed;

        break;
    }

    // fill in the jacobian
    int off;
    Real rho, mx, my, mz, alpha;
    Real q, m, r, r_dL;

    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        off = idx.solver;

        State &istate = GD::get_state(idx.global);

        int tp = istate.get_type();

        if (tp == +StateType::isField) {
            continue;
        }



        rho =   y0[off + +HydroState::ConsIdx::Density];
        mx =    y0[off + +HydroState::ConsIdx::Xmom];
        my =    y0[off + +HydroState::ConsIdx::Ymom];
        mz =    y0[off + +HydroState::ConsIdx::Zmom];
        alpha = y0[off + +HydroState::ConsIdx::Tracer]/rho;

        m = istate.get_mass(alpha);
        q = istate.get_charge(alpha);

        r = q/m;
        r_dL = r/dL;

        // internal

        JJ(off + +HydroState::ConsIdx::Xmom, off + +HydroState::ConsIdx::Ymom) =  Bz*r_dL;
        JJ(off + +HydroState::ConsIdx::Xmom, off + +HydroState::ConsIdx::Zmom) = -By*r_dL;

        JJ(off + +HydroState::ConsIdx::Ymom, off + +HydroState::ConsIdx::Xmom) = -Bz*r_dL;
        JJ(off + +HydroState::ConsIdx::Ymom, off + +HydroState::ConsIdx::Zmom) =  Bx*r_dL;

        JJ(off + +HydroState::ConsIdx::Zmom, off + +HydroState::ConsIdx::Xmom) =  By*r_dL;
        JJ(off + +HydroState::ConsIdx::Zmom, off + +HydroState::ConsIdx::Ymom) = -Bx*r_dL;

    }

    return 0;
}


