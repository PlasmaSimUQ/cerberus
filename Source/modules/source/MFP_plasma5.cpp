#include "MFP_plasma5.H"
#include "MFP_global.H"

using GD = GlobalData;

std::string Plasma5::tag = "plasma5";
bool Plasma5::registered = GetSourceTermFactory().Register(Plasma5::tag, SourceTermBuilder<Plasma5>);

Plasma5::Plasma5(){}

Plasma5::Plasma5(const sol::table &def)
{

    if (!GD::plasma_params_set) {
        Abort("Plasma parameters required for plasma5 source, define (>0) either skin_depth & beta or Larmor & Debye");
    }

    name = def.get<std::string>("name");

    if (!Plasma5::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &Plasma5::valid_state, includes, index);

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
//    if (linked_electron_idx.size() < 1 )
//        Abort("No electron fluid found for plasma5 source");
//    if (linked_ion_idx.size() < 1 )
//        Abort("No ion fluid found for plasma5 source");

}

Plasma5::~Plasma5()
{
    // do nothing
}


int Plasma5::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{
    BL_PROFILE("Plasma5::fun_rhs");
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
    Real rho, mx, my, mz, nrg, alpha;
    Real u, v, w;

    Real Larmor = GD::Larmor;
    Real Debye = GD::Debye;
    Real lightspeed = GD::lightspeed;




    // get charge and current density
    Real charge_density = 0.0;
    Real current_x = 0.0;
    Real current_y = 0.0;
    Real current_z = 0.0;

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
        nrg =   y0[idx.solver + +HydroState::ConsIdx::Eden];
        alpha = y0[idx.solver + +HydroState::ConsIdx::Tracer]/rho;

        m = istate.get_mass(alpha);
        q = istate.get_charge(alpha);

        r = q/m;

        charge_density += rho*r;
        current_x += r*mx;
        current_y += r*my;
        current_z += r*mz;

        u = mx/rho;
        v = my/rho;
        w = mz/rho;

        ydot[idx.solver + +HydroState::ConsIdx::Xmom] += (rho*r/Larmor)*(lightspeed*Ex + v*Bz - w*By);
        ydot[idx.solver + +HydroState::ConsIdx::Ymom] += (rho*r/Larmor)*(lightspeed*Ey + w*Bx - u*Bz);
        ydot[idx.solver + +HydroState::ConsIdx::Zmom] += (rho*r/Larmor)*(lightspeed*Ez + u*By - v*Bx);
        ydot[idx.solver + +HydroState::ConsIdx::Eden]  += (rho*r*lightspeed)/Larmor*(u*Ex + v*Ey + w*Ez);
    }

    // electric field and divergence constraint sources

    Real f1 = Larmor/(Debye*Debye*lightspeed);
    Real f2 = D_clean*D_clean*f1/lightspeed;

    ydot[field_offset + +FieldState::ConsIdx::Dx]   -= f1*current_x;
    ydot[field_offset + +FieldState::ConsIdx::Dy]   -= f1*current_y;
    ydot[field_offset + +FieldState::ConsIdx::Dz]   -= f1*current_z;
    ydot[field_offset + +FieldState::ConsIdx::phi] +=  f2*charge_density;


    return 0;
}

int Plasma5::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{
    BL_PROFILE("Plasma5::fun_jac");
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
    Real dD = GD::Debye;
    Real c0 = GD::lightspeed;

    Real c2 = c0*c0;
    Real dD2 = dD*dD;


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
    Real rho, mx, my, mz, nrg, alpha;
    Real q, m, r, r_dL, c0r_dL, cf1;

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
        nrg =   y0[off + +HydroState::ConsIdx::Eden];
        alpha = y0[off + +HydroState::ConsIdx::Tracer]/rho;

        m = istate.get_mass(alpha);
        q = istate.get_charge(alpha);

        r = q/m;
        r_dL = r/dL;
        c0r_dL = c0*r_dL;

        // internal

        JJ(off + +HydroState::ConsIdx::Xmom, off + +HydroState::ConsIdx::Ymom) =  Bz*r_dL;
        JJ(off + +HydroState::ConsIdx::Xmom, off + +HydroState::ConsIdx::Zmom) = -By*r_dL;

        JJ(off + +HydroState::ConsIdx::Ymom, off + +HydroState::ConsIdx::Xmom) = -Bz*r_dL;
        JJ(off + +HydroState::ConsIdx::Ymom, off + +HydroState::ConsIdx::Zmom) =  Bx*r_dL;

        JJ(off + +HydroState::ConsIdx::Zmom, off + +HydroState::ConsIdx::Xmom) =  By*r_dL;
        JJ(off + +HydroState::ConsIdx::Zmom, off + +HydroState::ConsIdx::Ymom) = -Bx*r_dL;

        // state - D

        JJ(off + +HydroState::ConsIdx::Xmom, field_offset + +FieldState::ConsIdx::Dx) = rho*c0r_dL;
        JJ(off + +HydroState::ConsIdx::Ymom, field_offset + +FieldState::ConsIdx::Dy) = rho*c0r_dL;
        JJ(off + +HydroState::ConsIdx::Zmom, field_offset + +FieldState::ConsIdx::Dz) = rho*c0r_dL;

        // D - state
        cf1 = -r*dL/(c0*dD2);
        JJ(field_offset + +FieldState::ConsIdx::Dx, off + +HydroState::ConsIdx::Xmom) = cf1;
        JJ(field_offset + +FieldState::ConsIdx::Dy, off + +HydroState::ConsIdx::Ymom) = cf1;
        JJ(field_offset + +FieldState::ConsIdx::Dz, off + +HydroState::ConsIdx::Zmom) = cf1;

        // nothing else varies due to changes in the following rows, but they should be updated with values at t+dt (non-linear)
        // best solved in conjunction with the 'rank_reduce' option

        // divergence cleaning source (damping is done in another source term)
        JJ(field_offset + +FieldState::ConsIdx::phi, off + +HydroState::ConsIdx::Density) = D_clean*D_clean*dL*r/(c2*dD2);

        // energy
        JJ(off + +HydroState::ConsIdx::Eden, off + +HydroState::ConsIdx::Xmom) = Ex*c0r_dL;
        JJ(off + +HydroState::ConsIdx::Eden, off + +HydroState::ConsIdx::Ymom) = Ey*c0r_dL;
        JJ(off + +HydroState::ConsIdx::Eden, off + +HydroState::ConsIdx::Zmom) = Ez*c0r_dL;

        JJ(off + +HydroState::ConsIdx::Eden, field_offset + +FieldState::ConsIdx::Dx) = mx*c0r_dL;
        JJ(off + +HydroState::ConsIdx::Eden, field_offset + +FieldState::ConsIdx::Dy) = my*c0r_dL;
        JJ(off + +HydroState::ConsIdx::Eden, field_offset + +FieldState::ConsIdx::Dz) = mz*c0r_dL;

    }

    return 0;
}

int Plasma5::face_src(Real x, Real y, Real z, Real t, Vector<Real> &y0, Array<Vector<Real>, AMREX_SPACEDIM> &ydot_lo, Array<Vector<Real>, AMREX_SPACEDIM> &ydot_hi) const
{
    BL_PROFILE("Plasma5::face_src");
    // J. Moreno, E. Oliva, P. Velarde, J.C.P. 2020, In Press

    Array<Real, AMREX_SPACEDIM> current;

    current.fill(0.0);

    Real alpha, rho, r, q, m;

    // calculate the current

    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        State &istate = GD::get_state(idx.global);
        int tt = istate.get_type();

        if (tt == +StateType::isField)
            continue;

        rho =   y0[idx.solver + +HydroState::ConsIdx::Density];
        alpha = y0[idx.solver + +HydroState::ConsIdx::Tracer]/rho;

        m = istate.get_mass(alpha);
        q = istate.get_charge(alpha);

        r = q/m;

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            current[d] += r*y0[idx.solver + +HydroState::ConsIdx::Xmom + d];
        }
    }

    Real f1 = -GD::Larmor/(GD::lightspeed*GD::Debye*GD::Debye);


    // update the face delta values
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        current[d] *= f1;

        for (const auto &idx : offsets) {

            State &istate = GD::get_state(idx.global);
            int tt = istate.get_type();

            if (tt != +StateType::isField)
                continue;

            // current sources go into the corresponding face both hi and lo
            ydot_lo[d][idx.solver + +FieldState::ConsIdx::Dx + d] += current[d];
            ydot_hi[d][idx.solver + +FieldState::ConsIdx::Dx + d] += current[d];
        }
    }

    return 0;
}

Real Plasma5::get_max_freq(Vector<Real> &y) const
{
    BL_PROFILE("Plasma5::get_max_freq");
    // get any magnetic field

    Real Bx=0, By=0, Bz=0;

    int field_offset;

    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        State &istate = GD::get_state(idx.global);
        int t = istate.get_type();

        if (t != +StateType::isField)
            continue;

        field_offset = idx.solver;

        // magnetic field
        Bx = y[field_offset + +FieldState::ConsIdx::Bx];
        By = y[field_offset + +FieldState::ConsIdx::By];
        Bz = y[field_offset + +FieldState::ConsIdx::Bz];

        break;

    }

    Real B = std::sqrt(Bx*Bx + By*By + Bz*Bz);


    Real q, m, r;
    Real rho, alpha;
    Real omega_p, omega_c;

    Real D2 = GD::Debye*GD::Debye;
    Real L = GD::Larmor;

    Real f = 0;

    for (const auto &idx : offsets) {

        State &istate = GD::get_state(idx.global);
        int t = istate.get_type();

        if (t == +StateType::isField)
            continue;

        rho =   y[idx.solver + +HydroState::ConsIdx::Density];
        alpha = y[idx.solver + +HydroState::ConsIdx::Tracer]/rho;

        m = istate.get_mass(alpha);
        q = istate.get_charge(alpha);

        r = q/m;

        omega_p = 10*std::sqrt(rho*r*r/D2)/(2*PI);

        omega_c = 10*(std::abs(r)*B)/(L*2*PI);

        f = std::max(f, omega_p);
        f = std::max(f, omega_c);

    }

    return f;
}

void Plasma5::calc_charge_density(const Box& box,
                                  const Real* prob_lo,
                                  const Real* dx,
                                  Real time,
                                  const Vector<FArrayBox*>& src,
                                  FArrayBox& cd,
                                  FArrayBox& J
                                  EB_OPTIONAL(,const Vector<const EBCellFlagFab*>& flag)
                                  ) const
{
    BL_PROFILE("Plasma5::calc_charge_density");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Real q, m, r;
    Real rho, alpha;
    Real mx, my, mz;

    cd.setVal(0.0);
    J.setVal(0.0);

    Array4<Real> const& cd4 = cd.array();
    Array4<Real> const& J4 = J.array();

    OffsetIndex field_idx;

    for (const auto &idx : offsets) {

        State &istate = GD::get_state(idx.global);
        int t = istate.get_type();

        if (t == +StateType::isField) {
            field_idx = idx;
            continue;
        }

        if (!idx.valid) continue;

#ifdef AMREX_USE_EB
        Array4<const EBCellFlag> const& f4 = flag[idx.local]->array();
#endif

        const int rho_idx = istate.get_cons_density_idx();
        const int mom_idx = istate.get_cons_vector_idx()[0];
        const int trc_idx = istate.get_cons_tracer_idx();

        Array4<const Real> const& src4 = src[idx.local]->array();

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                    if (f4(i,j,k).isCovered()) {
                        continue;
                    }
#endif

                    rho   = src4(i,j,k,rho_idx);



                    if (rho <= 0.0)
                        continue;

                    alpha = src4(i,j,k,trc_idx)/rho;

                    q = istate.get_charge(alpha);
                    m = istate.get_mass(alpha);
                    r = q/m;

                    cd4(i,j,k) += rho*r;

                    mx   = src4(i,j,k,mom_idx+0);
                    my   = src4(i,j,k,mom_idx+1);
                    mz   = src4(i,j,k,mom_idx+2);
                    J4(i,j,k,0) += mx*r;
                    J4(i,j,k,1) += my*r;
                    J4(i,j,k,2) += mz*r;
                }
            }
        }
    }

    // factor by the relative permittivity and permeability

    Real mu, ep;

    State &istate = GD::get_state(field_idx.global);

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag[field_idx.local]->array();
#endif

    Array4<const Real> const& src4 = src[field_idx.local]->array();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    continue;
                }
#endif

                mu = src4(i,j,k,+FieldState::ConsIdx::mu);
                ep = src4(i,j,k,+FieldState::ConsIdx::ep);

                cd4(i,j,k) /= ep;
                J4(i,j,k,0) *= mu;
                J4(i,j,k,1) *= mu;
                J4(i,j,k,2) *= mu;
            }
        }
    }

    return;
}

bool Plasma5::valid_state(const int global_idx)
{

    State &istate = GD::get_state(global_idx);

    switch (istate.get_type()) {
        case +StateType::isField:
            return true;
        case +StateType::isHydro:
            return true;
        case +StateType::isHydro2P:
            return true;
        default:
            return false;
    }
}

bool Plasma5::valid_solver(const int solve_idx)
{
    return true;
}

std::string Plasma5::print() const
{
    std::stringstream msg;

    msg << tag << " : " << name;

    msg << " (ions = [";
    bool first = true;
    for (const int i : linked_ion_idx) {
        if (not first) msg << ", ";
        msg << GD::get_state(i).name;
        first = false;
    }

    msg <<  "], electrons = [";
    first = true;
    for (const int i : linked_electron_idx) {
        if (not first) msg << ", ";
        msg << GD::get_state(i).name;
        first = false;
    }

    msg << "], EM = " + GD::get_state(linked_EM_idx).name + ")";

    return msg.str();

}


