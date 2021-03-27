#include "MFP_hydro_2p.H"
#include "MFP_global.H"
#include "MFP_state.H"

using GD = GlobalData;

std::string HydroTwoPressure::tag = "two_pressure";
bool HydroTwoPressure::registered = GetSourceTermFactory().Register(HydroTwoPressure::tag, SourceTermBuilder<HydroTwoPressure>);

HydroTwoPressure::HydroTwoPressure(){}

HydroTwoPressure::HydroTwoPressure(const sol::table& def)
{
    // only need slopes if we have non-zero charge

    name = def.get<std::string>("name");

    if (!HydroTwoPressure::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &HydroTwoPressure::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }


    linked_em = -1;

    Real q0, q1;
    for (const auto &idx : offsets) {
        State &istate = GD::get_state(idx.global);

        q0 = istate.charge[0];
        q1 = istate.charge[1];
        if (((std::abs(q0) == 0) || (std::abs(q1) == 0)) && (istate.get_type() != +StateType::isField))
            Abort("Source '"+tag+"' applied to state '"+istate.name+"' with zero charge");

        // also check if we have an active field
        if (istate.get_type() == +StateType::isField) {
            linked_em = istate.global_idx;
        }

        if (istate.get_type() == +StateType::isHydro2P) {
            linked_hydro.push_back(istate.global_idx);
        }
    }

    // get the reconstruction option
    set_reconstruction(def);
    if (!reconstruction)
        Abort("Source '"+tag+"' requires a reconstruction method to be defined (reconstruction=)");
}

HydroTwoPressure::~HydroTwoPressure()
{
    // do nothing
}

int HydroTwoPressure::num_slopes() const
{
    int cnt = linked_hydro.size()*(3*AMREX_SPACEDIM);

    return cnt;
}

// calculate slopes and pack them serially into a vector
void HydroTwoPressure::calc_slopes(const Box& box,
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

    FArrayBox buffer(slope_box, 3);
    Array4<Real> const& b4 = buffer.array();

    // calculate the slopes of the velocities

    for (const int& hydro_idx : linked_hydro) {

        Array4<const Real> const& h4 = src_dat[hydro_idx]->array();

        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {
                    const Real rho = h4(i,j,k,+Hydro2PState::ConsIdx::Density);
                    const Real mx = h4(i,j,k,+Hydro2PState::ConsIdx::Xmom);
                    b4(i,j,k,0) = mx/rho;

                    const Real my = h4(i,j,k,+Hydro2PState::ConsIdx::Ymom);
                    b4(i,j,k,1) = my/rho;

                    const Real mz = h4(i,j,k,+Hydro2PState::ConsIdx::Zmom);
                    b4(i,j,k,2) = mz/rho;

                }
            }
        }

        for (int d=0; d < AMREX_SPACEDIM; ++d) {
            for (int ui=0; ui<3; ++ui) {
                State::calc_slope(box, buffer, slopes[cnt], EB_OPTIONAL(*flags[hydro_idx],) dx, ui, d, *reconstruction); cnt++;
            }
        }

    }
}

void HydroTwoPressure::retrieve_slopes(
        Vector<FArrayBox>& slopes,
        const int i,
        const int j,
        const int k)
{

    slope.resize(num_slopes());
    int cnt = 0;

    for (const int& hydro_idx : linked_hydro) {
        for (int d=0; d < AMREX_SPACEDIM; ++d) {
            for (int ui=0; ui<3; ++ui) {
                slope[cnt] = slopes[cnt].array()(i,j,k); cnt++;

            }
        }
    }

}

int HydroTwoPressure::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{

    // get any magnetic field

    Real bx = 0.0;
    Real by = 0.0;
    Real bz = 0.0;

    int cnt = 0;
    for (const auto &idx : offsets) {

        State &istate = GD::get_state(idx.global);

        if (istate.get_type() != +StateType::isField)
            continue;

        if (!idx.valid) Abort("State '"+istate.name+"' is unavailable for source of type '"+tag+"'");

        // magnetic field
        Real xB = y0[idx.solver + +FieldState::ConsIdx::Bx];
        Real yB = y0[idx.solver + +FieldState::ConsIdx::By];
        Real zB = y0[idx.solver + +FieldState::ConsIdx::Bz];

        Real B = xB*xB + yB*yB + zB*zB;

        if (B <= 0.0) {
            B = 0.0;
        } else {
            B = 1/sqrt(B);
        }

        bx = xB*B;
        by = yB*B;
        bz = zB*B;
    }

    int hydro_cnt = cnt;
    for (const auto &idx : offsets) {

        if (!idx.valid) continue;

        State &istate = GD::get_state(idx.global);

        if (istate.get_type() == +StateType::isField)
            continue;

        Real tau = istate.pressure_relaxation_rate;

        Real rho   = y0[idx.solver + +Hydro2PState::ConsIdx::Density];
        Real mx    = y0[idx.solver + +Hydro2PState::ConsIdx::Xmom];
        Real my    = y0[idx.solver + +Hydro2PState::ConsIdx::Ymom];
        Real mz    = y0[idx.solver + +Hydro2PState::ConsIdx::Zmom];
        Real alpha = y0[idx.solver + +Hydro2PState::ConsIdx::Tracer]/rho;
        Real pp    = y0[idx.solver + +Hydro2PState::ConsIdx::PrsP];

        // calculate total pressure
        Real nrg   = y0[idx.solver + +Hydro2PState::ConsIdx::Eden];
        Real gam = istate.get_gamma(alpha);
        Real prs = (nrg - 0.5*(mx*mx + my*my + mz*mz)/rho)*(gam - 1);

        if (tau > 0.0) {
            // calculate relaxation
            ydot[idx.solver + +Hydro2PState::ConsIdx::PrsP] += (prs - pp)/tau;
        }

        Real q = istate.get_charge(alpha);
        if (q != 0.0) {

            Real u = mx/rho;
            Real v = my/rho;
            Real w = mz/rho;



            // grab the velocity gradients
            const Real du_dx = slope[hydro_cnt]; hydro_cnt++;
            const Real dv_dx = slope[hydro_cnt]; hydro_cnt++;
            const Real dw_dx = slope[hydro_cnt]; hydro_cnt++;

#if AMREX_SPACEDIM >= 2
            const Real du_dy = slope[hydro_cnt]; hydro_cnt++;
            const Real dv_dy = slope[hydro_cnt]; hydro_cnt++;
            const Real dw_dy = slope[hydro_cnt]; hydro_cnt++;
#else
            const Real du_dy = 0;
            const Real dv_dy = 0;
            const Real dw_dy = 0;
#endif

#if AMREX_SPACEDIM == 3
            const Real du_dz = slope[hydro_cnt]; hydro_cnt++;
            const Real dv_dz = slope[hydro_cnt]; hydro_cnt++;
            const Real dw_dz = slope[hydro_cnt]; hydro_cnt++;
#else
            const Real du_dz = 0;
            const Real dv_dz = 0;
            const Real dw_dz = 0;
#endif

//            const Real duvw = du_dx + dv_dy + dw_dz;

            // calculate the source terms
            const Real dpp_dt = 2*pp*(bx*(bx*du_dx + by*dv_dx + bz*dw_dx)
                                     +by*(bx*du_dy + by*dv_dy + bz*dw_dy)
                                     +bz*(bx*du_dz + by*dv_dz + bz*dw_dz));

            ydot[idx.solver + +Hydro2PState::ConsIdx::PrsP] -= dpp_dt;

        }

    }


    return 0;
}

int HydroTwoPressure::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{

    const int n_terms = y0.size();

    Vector<Real> ydot(n_terms);
    num_jac(x, y, z, t, y0, ydot, J);

    return 0;
}



bool HydroTwoPressure::valid_state(const int global_idx)
{
    State &istate = GD::get_state(global_idx);

    switch (istate.get_type()) {
    case +StateType::isField:
        return true;
    case +StateType::isHydro2P:
        return true;
    default:
        return false;
    }
}

bool HydroTwoPressure::valid_solver(const int solve_idx)
{
    if (solve_idx != +SolveType::Explicit) {
        return false;
    }
    return true;
}

