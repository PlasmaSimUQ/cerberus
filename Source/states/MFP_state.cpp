#include "MFP_state.H"

#include "MFP_fillbc.H"
#include "MFP_fillcc_eb.H"

State::State()
{
    BL_PROFILE("State::State");

    name = "null";
    global_idx = -1;
    data_idx = -1;
}

State::~State() {}

Vector<BoundaryInfo> State::get_bc_limits(const Box& box, const Geometry& geom) const
{
    BL_PROFILE("State::get_bc_limits");

    Vector<BoundaryInfo> limits;

    if (AMREX_D_TERM(geom.isPeriodic(0), &&geom.isPeriodic(1), &&geom.isPeriodic(2))) return limits;

    const Box& domain = geom.Domain();

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    const Dim3 domlo = amrex::lbound(domain);
    const Dim3 domhi = amrex::ubound(domain);

    // define the limits of our operations
    BoundaryInfo bi;
    bi.imin = lo.x;
    bi.imax = hi.x;
    bi.jmin = lo.y;
    bi.jmax = hi.y;
    bi.kmin = lo.z;
    bi.kmax = hi.z;

    int ilo = domlo.x;
    int ihi = domhi.x;

    if (!geom.isPeriodic(0)) {
        if (lo.x < ilo) {
            limits.push_back(bi);
            BoundaryInfo& b = limits.back();
            b.imin = lo.x;
            b.imax = ilo - 1;
            b.dir = 0;
            b.lo_hi = 0;
        }

        if (hi.x > ihi) {
            limits.push_back(bi);
            BoundaryInfo& b = limits.back();
            b.imin = ihi + 1;
            b.imax = hi.x;
            b.dir = 0;
            b.lo_hi = 1;
        }
    }

#if AMREX_SPACEDIM >= 2
    if (!geom.isPeriodic(1)) {
        int jlo = domlo.y;
        int jhi = domhi.y;

        if (lo.y < jlo) {
            limits.push_back(bi);
            BoundaryInfo& b = limits.back();
            b.jmin = lo.y;
            b.jmax = jlo - 1;
            b.dir = 1;
            b.lo_hi = 0;
        }

        if (hi.y > jhi) {
            limits.push_back(bi);
            BoundaryInfo& b = limits.back();
            b.jmin = jhi + 1;
            b.jmax = hi.y;
            b.dir = 1;
            b.lo_hi = 1;
        }
    }
#endif

#if AMREX_SPACEDIM == 3
    if (!geom.isPeriodic(2)) {
        int klo = domlo.z;
        int khi = domhi.z;

        if (lo.z < klo) {
            limits.push_back(bi);
            BoundaryInfo& b = limits.back();
            b.kmin = lo.z;
            b.kmax = klo - 1;
            b.dir = 2;
            b.lo_hi = 0;
        }

        if (hi.z > khi) {
            limits.push_back(bi);
            BoundaryInfo& b = limits.back();
            b.kmin = khi + 1;
            b.kmax = hi.z;
            b.dir = 2;
            b.lo_hi = 1;
        }
    }
#endif

    return limits;
}

#ifdef AMREX_USE_EB

void State::update_eb_vfrac(const Geometry& geom, FArrayBox& vfrac) const
{
    BL_PROFILE("State::update_eb_vfrac");

    const BCRec& bc = boundary_conditions.eb_bc;

    // get supporting info for calling fillcc
    const Real* dx = geom.CellSize();
    const Real* problo = geom.ProbLo();
    const int* lo = vfrac.loVect();
    const Box& domain = geom.Domain();
    const int* dom_lo = domain.loVect();
    Real xlo[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) { xlo[i] = problo[i] + dx[i] * (lo[i] - dom_lo[i]); }

    //    plot_FAB_2d(vfrac, 0, "vfrac before", false, false);

    Vector<BoundaryInfo> limits = get_bc_limits(vfrac.box(), geom);
    if (!limits.empty()) {
        Array4<Real> const& vfrac4 = vfrac.array();

        // if the cell we are copying from has no solid in it, don't copy
        const Real check = 1.0;

        for (size_t i = 0; i < limits.size(); ++i) {
            const auto L = limits[i];

            const IntVect box_lo(AMREX_D_DECL(L.imin, L.jmin, L.kmin));
            const IntVect box_hi(AMREX_D_DECL(L.imax, L.jmax, L.kmax));
            const Box box(box_lo, box_hi);
            fab_filcc_eb(box, vfrac4, 1, geom.Domain(), dx, xlo, &bc, check);
        }
    }
    //        plot_FAB_2d(vfrac, 0, "vfrac after", false, true);
}

void State::update_eb_flags(const Geometry& geom, EBCellFlagFab& flags) const
{
    BL_PROFILE("State::update_eb_flags");

    const BCRec& bc = boundary_conditions.eb_bc;

    // get supporting info for calling fillcc
    const Real* dx = geom.CellSize();
    const Real* problo = geom.ProbLo();
    const int* lo = flags.loVect();
    const Box& domain = geom.Domain();
    const int* dom_lo = domain.loVect();
    Real xlo[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) { xlo[i] = problo[i] + dx[i] * (lo[i] - dom_lo[i]); }

    Vector<BoundaryInfo> limits = get_bc_limits(flags.box(), geom);
    if (!limits.empty()) {
        Array4<EBCellFlag> const& flags4 = flags.array();

        // if the cell we are copying is regular then don't copy
        EBCellFlag check;
        check.setRegular();

        for (size_t i = 0; i < limits.size(); ++i) {
            const auto L = limits[i];

            const IntVect box_lo(AMREX_D_DECL(L.imin, L.jmin, L.kmin));
            const IntVect box_hi(AMREX_D_DECL(L.imax, L.jmax, L.kmax));
            const Box box(box_lo, box_hi);
            fab_filcc_eb(box, flags4, 1, geom.Domain(), dx, xlo, &bc, check);
        }
    }
}

#endif

void State::get_refinement_tags(MFP* mfp, TagBoxArray& tags)
{
    BL_PROFILE("State::get_refinement_tags");

    if (refine) refine->get_tags(mfp, tags);

#ifdef AMREX_USE_EB
    if (MFP::refine_cutcells) Refinement::tag_cut_cells(mfp, tags, global_idx);
#endif
}

void State::write_info(nlohmann::json& js) const
{
    BL_PROFILE("State::write_info");

    // write out stuff that is common to all states

    js["name"] = name;
    js["global_idx"] = global_idx;
    js["data_idx"] = data_idx;
}

std::string State::str() const
{
    BL_PROFILE("State::str");

    return "";
}

ClassFactory<State>& GetStateFactory()
{
    static ClassFactory<State> F;
    return F;
}
