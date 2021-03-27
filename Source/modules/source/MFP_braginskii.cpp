#include "MFP_braginskii.H"
#include "MFP_global.H"

using GD = GlobalData;

// INTENTIONALLY NOT REGISTERED
std::string BraginskiiSource::tag = "braginskii";
//bool BraginskiiSource::registered = GetSourceTermFactory().Register(BraginskiiSource::tag, SourceTermBuilder<BraginskiiSource>);

BraginskiiSource::BraginskiiSource(){}

BraginskiiSource::BraginskiiSource(const sol::table& def) 
{
    name = def.get<std::string>("name"); 

    if (!BraginskiiSource::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &BraginskiiSource::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }

    linked_em = -1;
    linked_hydro = -1;

    Real q0, q1;
    for (const auto &idx : offsets) {
        State &istate = GD::get_state(idx.global);

        q0 = istate.charge[0];
        q1 = istate.charge[1];
        if ((std::abs(q0) == 0) || (std::abs(q1) == 0))
            Abort("Source '"+tag+"' applied to state '"+istate.name+"' with zero charge");

        // also check if we have an active field
        if (istate.get_type() == +StateType::isField) {
            linked_em = istate.global_idx;
        }

        if (istate.get_type() == +StateType::isHydro) {
            if (linked_hydro > -1)
                Abort("Source '"+tag+"' has more than one hydro state linked to it");
            linked_hydro = istate.global_idx;
        }
    }

    // get the reconstruction option
    set_reconstruction(def);
    if (!reconstruction)
        Abort("Source '"+tag+"' requires a reconstruction method to be defined (reconstruction=)");

}

BraginskiiSource::~BraginskiiSource()
{
    // do nothing
}

int BraginskiiSource::num_slopes() const
{
    return 6*AMREX_SPACEDIM; // rho, u, v, w, p, T
}

// calculate slopes and pack them serially into a vector
void BraginskiiSource::calc_slopes(const Box& box,
                                   Vector<FArrayBox*>& src_dat,
                                   Vector<FArrayBox>& slopes,
                                   EB_OPTIONAL(Vector<const EBCellFlagFab*> &flags,)
                                   const Real *dx) const
{

    slopes.resize(num_slopes());
    int cnt = 0;

    const Box slope_box = grow(box, 1);

    const Dim3 lo = amrex::lbound(slope_box);
    const Dim3 hi = amrex::ubound(slope_box);

    // First calculate the primitive variables magnetic field
    int nc = +HydroState::ConsIdx::NUM;

    FArrayBox buffer(slope_box, +HydroState::PrimIdx::Prs);
    Array4<Real> const& b4 = buffer.array();

    State& istate = GD::get_state(linked_hydro);
    Array4<const Real> const& h4 = src_dat[linked_hydro]->array();
    Vector<Real> U(nc);

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

                // grab the conserved quantities
                for (int n=0; n<nc; ++n) {
                    U[n] = h4(i,j,k,n);
                }

                // convert to primitive
                istate.cons2prim(U);

                // load into buffer
                for (int n=0; n<+HydroState::PrimIdx::Temp; ++n) {
                    b4(i,j,k,n) = U[n];
                }
            }
        }
    }

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int n=0; n<+HydroState::PrimIdx::Temp; ++n) {
            State::calc_slope(box, buffer, slopes[cnt], EB_OPTIONAL(*flags[linked_hydro],) dx, n, d, *reconstruction); cnt++;
        }
    }
}

void BraginskiiSource::retrieve_slopes(
        Vector<FArrayBox>& slopes,
        const int i,
        const int j,
        const int k)
{

    slope.resize(num_slopes());
    int cnt = 0;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int n=0; n<+HydroState::PrimIdx::Prs; ++n) {
            slope[cnt] = slopes[cnt].array()(i,j,k); cnt++;
        }
    }
}

Vector<Real> BraginskiiSource::source(const Vector<Real>& y, const Vector<OffsetIndex> &index) const {
    /*Would like to have kept this general like Daryls above but BRaginskii stuff is very 
    prescriptive and only really accounts for a fully ionised simple plasma or a three 
    component plasam of ion, neutral, and electrons. The current set up is to allow reversion
    to the for loop arrangment. `a' represents the electrons, `b' represents the ions. */

    Vector<Real> ydot(y.size());

    return ydot;
}

int BraginskiiSource::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{

    ydot = source(y0, offsets);

    return 0;
}

int BraginskiiSource::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{

    return 0;
}

Real BraginskiiSource::get_max_freq(Vector<Real> &y) const{
    // get any magnetic field

    Real Bx=0., By=0., Bz=0.;

    int field_offset;

    for (const auto &idx : offsets) {
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
    if (B<0.) {
    amrex::Abort("Negative B field in Braginskii source");
    }

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

        if (!idx.valid) Abort("State '"+istate.name+"' is unavailable for source of type '"+tag+"'");

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

bool BraginskiiSource::valid_state(const int global_idx)
{

    State &istate = GD::get_state(global_idx);

    switch (istate.get_type()) {
    case +StateType::isField:
        return true;
    case +StateType::isHydro:
        return true;
    default:
        return false;
    }
}

bool BraginskiiSource::valid_solver(const int solver_idx)
{
    return true;
}

std::string BraginskiiSource::print() const
{
    std::stringstream msg;

    msg << tag << " : " << name;
    return msg.str();

}




