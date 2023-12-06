#ifdef SYMPLECTIC

    #include "MFP_chargedparticle.H"
    #include "MFP_field.H"
    #include "MFP_interpolation.H"
    #include "MFP_symplectic.H"

    #include <AMReX_AmrLevel.H>

using namespace amrex;

// Central difference first order
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Array1D<Real, 0, 2>
  curl_cdiff_1(Array4<Real const> const& a, int i, int j, int k, const double* ics)
{
    #if AMREX_SPACEDIM == 3
    return {((a(i, j + 1, k, Z) - a(i, j - 1, k, Z)) * ics[Y] * 0.5 -
             (a(i, j, k + 1, Y) - a(i, j, k - 1, Y)) * ics[Z] * 0.5),
            ((a(i, j, k + 1, X) - a(i, j, k - 1, X)) * ics[Z] * 0.5 -
             (a(i + 1, j, k, Z) - a(i - 1, j, k, Z)) * ics[X] * 0.5),
            ((a(i + 1, j, k, Y) - a(i - 1, j, k, Y)) * ics[X] * 0.5 -
             (a(i, j + 1, k, X) - a(i, j - 1, k, X)) * ics[Y] * 0.5)};
    #elif AMREX_SPACEDIM == 2
    return {(a(i, j + 1, k, Z) - a(i, j - 1, k, Z)) * ics[Y] * 0.5,
            -(a(i + 1, j, k, Z) - a(i - 1, j, k, Z)) * ics[X] * 0.5,
            ((a(i + 1, j, k, Y) - a(i - 1, j, k, Y)) * ics[X] * 0.5 -
             (a(i, j + 1, k, X) - a(i, j - 1, k, X)) * ics[Y] * 0.5)};
    #else
    return {0,
            -(a(i + 1, j, k, Z) - a(i - 1, j, k, Z)) * ics[X] * 0.5,
            (a(i + 1, j, k, Y) - a(i - 1, j, k, Y)) * ics[X] * 0.5};
    #endif
}

// forward difference, first order
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Array1D<Real, 0, 2>
  curl_fdiff_1(Array4<Real const> const& a, int i, int j, int k, const double* ics)
{
    #if AMREX_SPACEDIM == 3
    return {
      ((a(i, j + 1, k, Z) - a(i, j, k, Z)) * ics[Y] - (a(i, j, k + 1, Y) - a(i, j, k, Y)) * ics[Z]),
      ((a(i, j, k + 1, X) - a(i, j, k, X)) * ics[Z] - (a(i + 1, j, k, Z) - a(i, j, k, Z)) * ics[X]),
      ((a(i + 1, j, k, Y) - a(i, j, k, Y)) * ics[X] -
       (a(i, j + 1, k, X) - a(i, j, k, X)) * ics[Y])};
    #elif AMREX_SPACEDIM == 2
    return {(a(i, j + 1, k, Z) - a(i, j, k, Z)) * ics[Y],
            -(a(i + 1, j, k, Z) - a(i, j, k, Z)) * ics[X],
            ((a(i + 1, j, k, Y) - a(i, j, k, Y)) * ics[X] -
             (a(i, j + 1, k, X) - a(i, j, k, X)) * ics[Y])};
    #else
    return {0,
            -(a(i + 1, j, k, Z) - a(i, j, k, Z)) * ics[X],
            (a(i + 1, j, k, Y) - a(i, j, k, Y)) * ics[X]};
    #endif
}

// backward difference, first order
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Array1D<Real, 0, 2>
  curl_bdiff_1(Array4<Real const> const& a, int i, int j, int k, const double* ics)
{
    #if AMREX_SPACEDIM == 3
    return {
      ((a(i, j, k, Z) - a(i, j - 1, k, Z)) * ics[Y] - (a(i, j, k, Y) - a(i, j, k - 1, Y)) * ics[Z]),
      ((a(i, j, k, X) - a(i, j, k - 1, X)) * ics[Z] - (a(i, j, k, Z) - a(i - 1, j, k, Z)) * ics[X]),
      ((a(i, j, k, Y) - a(i - 1, j, k, Y)) * ics[X] -
       (a(i, j, k, X) - a(i, j - 1, k, X)) * ics[Y])};
    #elif AMREX_SPACEDIM == 2
    return {(a(i, j, k, Z) - a(i, j - 1, k, Z)) * ics[Y],
            -(a(i, j, k, Z) - a(i - 1, j, k, Z)) * ics[X],
            (a(i, j, k, Y) - a(i - 1, j, k, Y)) * ics[X] -
              (a(i, j, k, X) - a(i, j - 1, k, X)) * ics[Y]};
    #else
    return {0,
            -(a(i, j, k, Z) - a(i - 1, j, k, Z)) * ics[X],
            (a(i, j, k, Y) - a(i - 1, j, k, Y)) * ics[X]};
    #endif
}

void E_curl(
  const Geometry geom, Box const& bx, Array4<Real const> const& E, Array4<Real> const& B, double dt)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        auto res = curl_fdiff_1(E, i, j, k, geom.InvCellSize());
        B(i, j, k, X) = B(i, j, k, X) - dt * res(0);
        B(i, j, k, Y) = B(i, j, k, Y) - dt * res(1);
        B(i, j, k, Z) = B(i, j, k, Z) - dt * res(2);
    });
}

void B_curl(
  const Geometry geom, Box const& bx, Array4<Real const> const& B, Array4<Real> const& E, double dt)
{
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        auto res = curl_bdiff_1(B, i, j, k, geom.InvCellSize());
        E(i, j, k, X) = E(i, j, k, X) + dt * res(0);
        E(i, j, k, Y) = E(i, j, k, Y) + dt * res(1);
        E(i, j, k, Z) = E(i, j, k, Z) + dt * res(2);
    });
}

void push_B_E(
  const Geometry geom, Box const& bx, Array4<Real> const& B, Array4<const Real> const& E, double dt)
{
    push_ff<E_curl>(geom, bx, E, B, dt);
}

void push_E_B(
  const Geometry geom, Box const& bx, Array4<Real> const& E, Array4<const Real> const& B, double dt)
{
    push_ff<B_curl>(geom, bx, B, E, dt);
}

void push_V_E(
  CParticleType* particles, long np, const Geometry geom, Array4<Real const> const& E, double dt)
{
    // Basis functions are supported over two cells in each direction

    const auto _Ics = geom.InvCellSize();
    GpuArray<Real, AMREX_SPACEDIM> lb;
    GpuArray<Real, AMREX_SPACEDIM> Ics;

    AMREX_D_TERM(lb[X] = geom.ProbLo(X); Ics[X] = _Ics[X];, lb[Y] = geom.ProbLo(Y);
                 Ics[Y] = _Ics[Y];
                 , lb[Z] = geom.ProbLo(Z);
                 Ics[Z] = _Ics[Z];)

    ParallelFor(np, [=] AMREX_GPU_DEVICE(long i) {
        const double m = particles[i].rdata(+ChargedParticle::ParticleIdxR::Mass);
        const double q = particles[i].rdata(+ChargedParticle::ParticleIdxR::Charge);
        const double coef = dt * q / m;

        // Particle coordinate shifted to cell
        Real pcoord[3] = {0, 0, 0};
        AMREX_D_TERM(pcoord[X] = (particles[i].pos(X) - lb[X]) * Ics[X];
                     , pcoord[Y] = (particles[i].pos(Y) - lb[Y]) * Ics[Y];
                     , pcoord[Z] = (particles[i].pos(Z) - lb[Z]) * Ics[Z];)

        int coord[3] = {0, 0, 0};
        AMREX_D_TERM(coord[X] = floor(pcoord[0]);, coord[Y] = floor(pcoord[1]);
                     , coord[Z] = floor(pcoord[2]);)

        double dvx = 0;
        double dvy = 0;
        double dvz = 0;

        constexpr int W1_r = WRANGE * 2;
        constexpr int W1_li = -WRANGE + 1;
        constexpr int W1_hi = WRANGE;
        constexpr int Wp_li = W1_li;
        constexpr int Wp_hi = W1_hi - 1;

        Array1D<Real, 0, W1_r> W1X = {0};
        Array1D<Real, 0, W1_r> W1Y = {0};
        Array1D<Real, 0, W1_r> W1Z = {0};
        Array1D<Real, 0, W1_r> WpX = {0};
        Array1D<Real, 0, W1_r> WpY = {0};
        Array1D<Real, 0, W1_r> WpZ = {0};

        int idx = 0;
        for (int i = W1_li; i <= W1_hi; i++) {
            auto cx = coord[X] + i;
            W1X(idx) = W1(pcoord[0] - cx);
            auto cy = coord[Y] + i;
            W1Y(idx) = W1(pcoord[1] - cy);
            auto cz = coord[Z] + i;
            W1Z(idx) = W1(pcoord[2] - cz);
            idx++;
        }
        idx = 0;
        for (int i = Wp_li; i <= Wp_hi; i++) {
            auto cx = coord[X] + i;
            WpX(idx) = Wp(pcoord[0] - cx);
            auto cy = coord[Y] + i;
            WpY(idx) = Wp(pcoord[1] - cy);
            auto cz = coord[Z] + i;
            WpZ(idx) = Wp(pcoord[2] - cz);
            idx++;
        }

        int idk = 0;
        for (int k = W1_li; k <= W1_hi; k++) {
            auto cz = coord[Z] + k;
            cz = AMREX_D_PICK(0, 0, cz);
            int idj = 0;
            for (int j = W1_li; j <= W1_hi; j++) {
                int idi = 0;
                auto cy = coord[Y] + j;
                cy = AMREX_D_PICK(0, cy, cy);
                for (int i = W1_li; i <= W1_hi; i++) {
                    auto cx = coord[X] + i;
                    dvx += E(cx, cy, cz, X) * WpX(idi) * W1Y(idj) * W1Z(idk);
                    dvy += E(cx, cy, cz, Y) * W1X(idi) * WpY(idj) * W1Z(idk);
                    dvz += E(cx, cy, cz, Z) * W1X(idi) * W1Y(idj) * WpZ(idk);
                    idi++;
                }
                idj++;
            }
            idk++;
        }
        particles[i].rdata(+ChargedParticle::ParticleIdxR::VX) += dvx * coef;
        particles[i].rdata(+ChargedParticle::ParticleIdxR::VY) += dvy * coef;
        particles[i].rdata(+ChargedParticle::ParticleIdxR::VZ) += dvz * coef;
    });
}

void Symplectic::Theta_B(MFP* mfp, Real time, Real dt) const
{
    BL_PROFILE("Symplectic::Theta_B");

    const amrex::BoxArray& grids = mfp->boxArray();
    const DistributionMapping& dmap = mfp->DistributionMap();
    const Geometry& geom = mfp->Geom();

    MultiFab fields_local(grids,
                          dmap,
                          +FieldDef::ConsIdx::NUM,
                          interpolation_range + 1,
                          MFInfo(),
                          mfp->Factory());
    mfp->FillPatch(*mfp,
                   fields_local,
                   interpolation_range + 1,
                   time,
                   field->data_idx,
                   0,
                   +FieldDef::ConsIdx::NUM);
    MultiFab B = MultiFab(fields_local, MakeType::make_alias, +FieldDef::ConsIdx::Bx, 3);
    MultiFab E = MultiFab(fields_local, MakeType::make_alias, +FieldDef::ConsIdx::Dx, 3);

    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);
    MultiFab& fields_global = mfp->get_new_data(field->data_idx);

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {
        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.validbox();
        auto const& B4 = B.array(mfi);
        auto const& E4 = E.array(mfi);
        push_E_B(geom, box, E4, B4, dt);

        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }

    // only copy what we have updated
    MultiFab::Copy(fields_global,
                   fields_local,
                   +FieldDef::ConsIdx::Dx,
                   +FieldDef::ConsIdx::Dx,
                   3,
                   0);
}

void Symplectic::Theta_E(MFP* mfp, Real time, Real dt) const
{
    BL_PROFILE("Symplectic::Theta_E");

    const amrex::BoxArray& grids = mfp->boxArray();
    const DistributionMapping& dmap = mfp->DistributionMap();
    const Geometry& geom = mfp->Geom();
    const int level = mfp->get_level();

    MultiFab fields_local(grids,
                          dmap,
                          +FieldDef::ConsIdx::NUM,
                          interpolation_range + 1,
                          MFInfo(),
                          mfp->Factory());
    mfp->FillPatch(*mfp,
                   fields_local,
                   interpolation_range + 1,
                   time,
                   field->data_idx,
                   0,
                   +FieldDef::ConsIdx::NUM);
    MultiFab B = MultiFab(fields_local, MakeType::make_alias, +FieldDef::ConsIdx::Bx, 3);
    MultiFab E = MultiFab(fields_local, MakeType::make_alias, +FieldDef::ConsIdx::Dx, 3);

    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);
    MultiFab& fields_global = mfp->get_new_data(field->data_idx);

    for (auto& P : species) {
        for (CParIterType pti(*(P->particles.get()), level); pti.isValid(); ++pti) {
            Real wt = ParallelDescriptor::second();

            CParticleType* AMREX_RESTRICT pc = &(pti.GetArrayOfStructs()[0]);
            const int np = pti.numParticles();
            Array4<Real const> const& E4 = E.const_array(pti);
            push_V_E(pc, np, geom, E4, dt);

            Box box = pti.fabbox();
            wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
            cost[pti].plus(wt, box);
        }
    }

    for (MFIter mfi(E); mfi.isValid(); ++mfi) {
        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.validbox();
        Array4<Real const> const& E4 = E.const_array(mfi);
        Array4<Real> const& B4 = B.array(mfi);
        push_B_E(geom, box, B4, E4, dt);

        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }

    // only copy what we have updated
    MultiFab::Copy(fields_global,
                   fields_local,
                   +FieldDef::ConsIdx::Bx,
                   +FieldDef::ConsIdx::Bx,
                   3,
                   0);
}

AMREX_GPU_HOST_DEVICE inline bool
  segment_reflect(int isPeriodic, int smallEnd, int bigEnd, Real* seg_points, int* seg_idx)
{
    if (!isPeriodic && ((seg_idx[1] == smallEnd + WRANGE) || (seg_idx[1] == bigEnd - WRANGE))) {
        seg_idx[1] = seg_idx[0];
        seg_points[2] = 2 * seg_points[1] - seg_points[2];
        return true;
    }
    return false;
}

AMREX_GPU_HOST_DEVICE inline void particle_reflect(CParticleType* p, Real* seg_points, int comp)
{
    p->pos(comp) = seg_points[2];
    p->rdata(comp + 2) *= -1;
}

void Theta(CParticleType* particles,
           long np,
           const Geometry geom,
           Array4<Real> const& E,
           Array4<Real> const& B,
           Box bx,
           Real dt,
           int comp,
           bool generate_field)
{
    const auto _Ics = geom.InvCellSize();
    constexpr int W1_r = WRANGE * 2;
    constexpr int W1_li = -WRANGE + 1;
    constexpr int W1_hi = WRANGE;
    constexpr int Wp_li = W1_li;
    constexpr int Wp_hi = W1_hi - 1;
    constexpr int Wp_r = W1_r - 1;

    auto comp_u = (comp + 1) % 3;
    auto comp_l = (comp + 2) % 3;
    auto comp_Cs = geom.CellSize(comp);
    GpuArray<Real, 3> lower_bound = {0, 0, 0};
    GpuArray<Real, 3> Ics = {1, 1, 1};
    AMREX_D_TERM(lower_bound[X] = geom.ProbLo(X); Ics[X] = _Ics[X];
                 , lower_bound[Y] = geom.ProbLo(Y);
                 Ics[Y] = _Ics[Y];
                 , lower_bound[Z] = geom.ProbLo(Z);
                 Ics[Z] = _Ics[Z];)

    int isPeriodic, smallEnd, bigEnd;
    if (comp < AMREX_SPACEDIM) {
        isPeriodic = geom.isPeriodic(comp);
        smallEnd = geom.Domain().smallEnd(comp);
        bigEnd = geom.Domain().bigEnd(comp);
    } else {
        isPeriodic = 1;
        smallEnd = 0;
        bigEnd = 1;
    }

    // Remote particle grid lower corner
    //    ParallelFor(np,
    //                [=] AMREX_GPU_DEVICE (long i)
    for (int i = 0; i < np; ++i) {
        const Real m = particles[i].rdata(+ChargedParticle::ParticleIdxR::Mass);
        const Real q = particles[i].rdata(+ChargedParticle::ParticleIdxR::Charge);
        const Real B_coef = q / m * comp_Cs;
        const Real E_coef = q * Ics[X] * Ics[Y] * Ics[Z] * comp_Cs;

        Real old_pos;
        if (comp < AMREX_SPACEDIM) {
            old_pos = particles[i].pos(comp);
        } else {
            old_pos = 0.0;
        }

        Real new_pos = old_pos + dt * particles[i].rdata(comp + 2);

        Array<int, 3> coord = {0, 0, 0};
        AMREX_D_TERM(coord[X] = floor((particles[i].pos(X) - lower_bound[X]) * Ics[X]);
                     , coord[Y] = floor((particles[i].pos(Y) - lower_bound[Y]) * Ics[Y]);
                     , coord[Z] = floor((particles[i].pos(Z) - lower_bound[Z]) * Ics[Z]);)

        Real res_c1 = 0;
        Real res_c2 = 0;

        Array1D<Real, 0, W1_r> comp_uW1 = {0};
        Array1D<Real, 0, W1_r> comp_lW1 = {0};
        Array1D<Real, 0, W1_r> comp_uWp = {0};
        Array1D<Real, 0, W1_r> comp_lWp = {0};

        Real nl;
        if (comp_l < AMREX_SPACEDIM) {
            nl = (particles[i].pos(comp_l) - lower_bound[comp_l]) * Ics[comp_l];
        } else {
            nl = 0.0;
        }
        int idx = 0;
        for (int l = W1_li; l <= W1_hi; l++) {
            auto cl = coord[comp_l] + (l);
            comp_lW1(idx) = W1(nl - cl);
            idx++;
        }
        idx = 0;
        for (int l = Wp_li; l <= Wp_hi; l++) {
            auto cl = coord[comp_l] + (l);
            comp_lWp(idx) = Wp(nl - cl);
            idx++;
        }
        idx = 0;
        Real nu;
        if (comp_u < AMREX_SPACEDIM) {
            nu = (particles[i].pos(comp_u) - lower_bound[comp_u]) * Ics[comp_u];
        } else {
            nu = 0.0;
        }
        for (int u = W1_li; u <= W1_hi; u++) {
            auto cu = coord[comp_u] + (u);
            comp_uW1(idx) = W1(nu - cu);
            idx++;
        }
        idx = 0;

        for (int u = Wp_li; u <= Wp_hi; u++) {
            auto cu = coord[comp_u] + (u);
            comp_uWp(idx) = Wp(nu - cu);
            idx++;
        }

        Array<Real, 3> seg_points;
        Array<int, 2> seg_idx;

        // calculate the coordinates and cell indexes of the particle path
        seg_idx[0] = coord[comp];                                       // starting cell index
        seg_idx[1] = floor((new_pos - lower_bound[comp]) * Ics[comp]);  // ending cell index

        int diff = seg_idx[1] - seg_idx[0];

        if (std::abs(diff) > 1) Abort("Particles has moved too many cells!");

        int num_segments = amrex::Math::abs(diff) + 1;

        seg_points[0] = old_pos;

        if (num_segments == 2) {
            seg_points[1] = lower_bound[comp] + (seg_idx[0] + (diff + 1) / 2) * comp_Cs;
            seg_points[2] = new_pos;
        } else {
            seg_points[1] = new_pos;
        }

        // handle periodicity/reflection of the particle
        bool out_not_periodic = false;
        if ((isPeriodic == 0) && ((seg_idx[1] < smallEnd) || (seg_idx[1] > bigEnd))) {
            seg_idx[1] = seg_idx[0];
            seg_points[2] = 2 * seg_points[1] - seg_points[2];
            out_not_periodic = true;
        }

        for (int seg = 0; seg < num_segments; seg++) {
            coord[comp] = seg_idx[seg];
            auto i_s = seg_points[seg];
            auto i_e = seg_points[seg + 1];
            Array1D<Real, 0, Wp_r> compI_W12 = {0};

            auto ncs = (i_s - lower_bound[comp]) * Ics[comp];
            auto nce = (i_e - lower_bound[comp]) * Ics[comp];
            idx = 0;
            for (int c = Wp_li; c <= Wp_hi; c++) {
                auto cc = coord[comp] + (c);
                compI_W12(idx) = I_Wp(ncs - cc, nce - cc);
                idx++;
            }

            int idl = 0;
            for (int l = W1_li; l <= W1_hi; l++) {
                int idu = 0;
                for (int u = W1_li; u <= W1_hi; u++) {
                    int idc = 0;
                    Real mul = comp_lW1(idl) * comp_uW1(idu);
                    Real mulu12l1 = comp_uWp(idu) * comp_lW1(idl);
                    Real mull12u1 = comp_lWp(idl) * comp_uW1(idu);
                    for (int c = Wp_li; c <= Wp_hi; c++) {
                        int cx, cy, cz;

                        if (comp == 0) {
                            cx = coord[X] + (c);
                            cy = coord[Y] + (u);
                            cz = coord[Z] + (l);
                        } else if (comp == 1) {
                            cx = coord[X] + (l);
                            cy = coord[Y] + (c);
                            cz = coord[Z] + (u);
                        } else {
                            cx = coord[X] + (u);
                            cy = coord[Y] + (l);
                            cz = coord[Z] + (c);
                        }

                        cy = AMREX_D_PICK(0, cy, cy);
                        cz = AMREX_D_PICK(0, 0, cz);

                        if (generate_field) {
                            Gpu::Atomic::Add(&E(cx, cy, cz, comp), -E_coef * mul * compI_W12(idc));
                        }
                        res_c1 += B(cx, cy, cz, comp_u) * mull12u1 * compI_W12(idc);
                        res_c2 -= B(cx, cy, cz, comp_l) * compI_W12(idc) * mulu12l1;
                        idc++;
                    }
                    idu++;
                }
                idl++;
            }
        }

        // Update position
        if (out_not_periodic) {
            // reflect particle position and velocity
            particles[i].pos(comp) = seg_points[2];
            particles[i].rdata(comp + 2) *= -1;
        }

        // .Redistribute should wrap the particles if periodic
        else if (comp < AMREX_SPACEDIM) {
            particles[i].pos(comp) += dt * particles[i].rdata(2 + comp);
        }
        // B Vel update
        particles[i].rdata((comp + 2) % 3 + 2) += B_coef * res_c1;
        particles[i].rdata((comp + 1) % 3 + 2) += B_coef * res_c2;

    }  //);
}

void Symplectic::Theta_R(MFP* mfp, Real time, Real dt, int comp) const
{
    BL_PROFILE("Symplectic::Theta_R");

    const amrex::BoxArray& grids = mfp->boxArray();
    const DistributionMapping& dmap = mfp->DistributionMap();
    const Geometry& geom = mfp->Geom();
    const int level = mfp->get_level();

    MultiFab fields_local(grids,
                          dmap,
                          +FieldDef::ConsIdx::NUM,
                          interpolation_range + 1,
                          MFInfo(),
                          mfp->Factory());
    mfp->FillPatch(*mfp,
                   fields_local,
                   interpolation_range + 1,
                   time,
                   field->data_idx,
                   0,
                   +FieldDef::ConsIdx::NUM);
    MultiFab B = MultiFab(fields_local, MakeType::make_alias, +FieldDef::ConsIdx::Bx, 3);
    MultiFab E = MultiFab(fields_local, MakeType::make_alias, +FieldDef::ConsIdx::Dx, 3);

    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);
    MultiFab& fields_global = mfp->get_new_data(field->data_idx);

    for (auto& P : species) {
        for (CParIterType pti(*(P->particles.get()), level); pti.isValid(); ++pti) {
            Real wt = ParallelDescriptor::second();

            auto box = E.box(pti.index());

            // Each grid,tile has a their own local particle container
            CParticleType* AMREX_RESTRICT p = &(pti.GetArrayOfStructs()[0]);
            const int np = pti.numParticles();

            FArrayBox& bfab = B[pti];

            Array4<Real> const& B_loc = bfab.array();
            Array4<Real> const& E_loc = E[pti].array();
            Theta(p, np, geom, E_loc, B_loc, box, dt, comp, allow_field_generation);

            wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
            cost[pti].plus(wt, box);
        }

        P->redistribute();
    }

    // only copy what we have updated
    MultiFab::Copy(fields_global,
                   fields_local,
                   +FieldDef::ConsIdx::Dx,
                   +FieldDef::ConsIdx::Dx,
                   3,
                   0);
}
#endif
