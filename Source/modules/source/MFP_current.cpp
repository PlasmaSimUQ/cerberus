#include "MFP_current.H"

#include "eigen.hpp"
#include "forward.hpp"

using GD = GlobalData;

//---------------------------------------------------------------------------------------------

std::string CurrentSource::tag = "current";
bool CurrentSource::registered = GetSourceTermFactory().Register(CurrentSource::tag, SourceTermBuilder<CurrentSource>);

CurrentSource::CurrentSource(){}

CurrentSource::CurrentSource(const sol::table &def)
{

    name = def.get<std::string>("name");

    if (!CurrentSource::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &CurrentSource::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }

    // get our current source

    get_udf(def["x"], current[0], 0.0);
    get_udf(def["y"], current[1], 0.0);
    get_udf(def["z"], current[2], 0.0);

    return;
}

CurrentSource::~CurrentSource()
{
    // do nothing
}

Array<Real,3> CurrentSource::get_current(Real x, Real y, Real z, Real t) const
{
    BL_PROFILE("CurrentSource::get_current");
    // calculate the current at this point in time and space
    Array<Real,3> j;

    std::map<std::string, Real> Q{{"x",x}, {"y",y}, {"z",z}, {"t",t}};

    Real f1 = -GD::Larmor/(GD::lightspeed*GD::Debye*GD::Debye);

    for (int i = 0; i<3; ++i) {
        j[i] = f1*current[i](Q);
    }

    return j;
}

Vector<dual> CurrentSource::apply_current(const Vector<dual> &y0,
                                          const Vector<OffsetIndex> &offsets,
                                          const Array<Real,3> &j)
{
    BL_PROFILE("CurrentSource::apply_current");
    Vector<dual> ydot(y0.size());

    for (const auto &idx : offsets) {
        if (!idx.valid) continue;
        // electric field
        for (int i = 0; i<3; ++i) {
            ydot[idx.solver + +FieldState::ConsIdx::Dx + i] += j[i];
        }
    }

    return ydot;
}


int CurrentSource::face_src(Real x, Real y, Real z, Real t, Vector<Real> &y0, Array<Vector<Real>, AMREX_SPACEDIM> &ydot_lo, Array<Vector<Real>, AMREX_SPACEDIM> &ydot_hi) const
{

    BL_PROFILE("CurrentSource::face_src");
    Array<Real,3> j = get_current(x,y,z,t);

    Real alpha, rho, r, q, m;

    // update the face delta values
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        for (const auto &idx : offsets) {
            if (!idx.valid) continue;

            State &istate = GD::get_state(idx.global);
            int tt = istate.get_type();

            if (tt != +StateType::isField)
                continue;

            // current sources go into the corresponding face both hi and lo
            ydot_lo[d][idx.solver + +FieldState::ConsIdx::Dx + d] += j[d];
            ydot_hi[d][idx.solver + +FieldState::ConsIdx::Dx + d] += j[d];
        }
    }

    return 0;
}

int CurrentSource::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{
    BL_PROFILE("CurrentSource::fun_rhs");

    Array<Real,3> j = get_current(x,y,z,t);

    const int n_terms = y0.size();

    // copy to eigen vector
    Vector<dual> yy(n_terms);
    for (int i=0; i<n_terms; ++i) {
        yy[i] = y0[i];
    }

    // call source function
    Vector<dual> yd = apply_current(yy, offsets, j);

    // copy to ydot
    for (int i=0; i<n_terms; ++i) {
        ydot[i] = yd[i].val;
    }

    return 0;
}

void CurrentSource::calc_charge_density(const Box& box,
                                        const Real* prob_lo,
                                        const Real* dx,
                                        Real time,
                                        const Vector<FArrayBox*>& src,
                                        FArrayBox& cd,
                                        FArrayBox& J
                                        EB_OPTIONAL(,const Vector<const EBCellFlagFab*>& flag)
                                        ) const
{
    BL_PROFILE("CurrentSource::calc_charge_density");

    cd.setVal(0.0);
    J.setVal(0.0);

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);


    Array4<Real> const& cd4 = cd.array();
    Array4<Real> const& J4 = J.array();
    Real x, y, z;


    for     (int k = lo.z; k <= hi.z; ++k) {
        z = prob_lo[2] + (k + 0.5)*dx[2];
        for   (int j = lo.y; j <= hi.y; ++j) {
            y = prob_lo[1] + (j + 0.5)*dx[1];
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                x = prob_lo[0] + (i + 0.5)*dx[0];

                const Array<Real,3> J = get_current(x, y, z, time);

                // setting the charge density doesn't seem to make any sense here
                // what is the velocity??
//                cd4(i,j,k) = std::sqrt(J[0]*J[0] + J[1]*J[1] + J[2]*J[2]);

                J4(i,j,k,0) = J[0];
                J4(i,j,k,1) = J[1];
                J4(i,j,k,2) = J[2];

            }
        }
    }


    OffsetIndex field_idx;

    for (const auto &idx : offsets) {
        State &istate = GD::get_state(idx.global);

        if (istate.get_type() == +StateType::isField) {

            if (!idx.valid) Abort("State '"+istate.name+"' is unavailable for source of type '"+tag+"'");

            field_idx = idx;
            break;
        }
    }


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
}



bool CurrentSource::valid_state(const int global_idx)
{
    State &istate = GD::get_state(global_idx);

    switch (istate.get_type()) {
        case +StateType::isField:
            return true;
        default:
            return false;
    }
}

bool CurrentSource::valid_solver(const int solve_idx)
{
    if (solve_idx != +SolveType::Explicit) {
        return false;
    }
    return true;
}
