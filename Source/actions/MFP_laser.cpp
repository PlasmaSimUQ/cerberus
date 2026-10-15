#include "MFP_laser.H"

#include "MFP.H"
#include "MFP_state.H"

#include <AMReX_ParmParse.H>

std::string Laser::tag = "laser";
bool Laser::registered = GetActionFactory().Register(Laser::tag, ActionBuilder<Laser>);

Laser::Laser() {}
Laser::~Laser() {}

Laser::Laser(const int idx, const sol::table& def)
{
    BL_PROFILE("Laser::Laser");

    action_idx = idx;
    name = def["name"];

    const std::string state_name = def.get_or<std::string>("state", "");
    if (state_name.empty() || (MFP::get_state(state_name).get_type() != State::StateType::Hydro))
        Abort("Action '" + name + "' requires 'state' to name a hydro state");
    state = &HydroState::get_state(state_name);
    state_indexes.push_back(state->global_idx);

    {
        ParmParse pp("amr");
        int max_level = 0;
        pp.query("max_level", max_level);
        if (max_level > 0) Abort("Action '" + name + "': laser supports amr.max_level = 0 only");
    }

    // geometry (lab units in, code units out)
    const sol::table o = def["origin"];
    const sol::table d = def["direction"];
    Real norm = 0.0;
    for (int i = 0; i < 3; ++i) {
        origin[i] = o[i + 1].get_or(0.0) / MFP::x_ref;
        dir[i] = d[i + 1].get_or(0.0);
        norm += dir[i] * dir[i];
        if ((i >= AMREX_SPACEDIM) && (dir[i] != 0.0))
            Abort("Action '" + name + "': 'direction' has a component outside the simulation dimension");
    }
    if (norm <= 0.0) Abort("Action '" + name + "': 'direction' must be a non-zero vector");
    norm = std::sqrt(norm);
    for (int i = 0; i < 3; ++i) dir[i] /= norm;

    radius = def.get_or("radius", 0.0) / MFP::x_ref;
    if (radius < 0.0) Abort("Action '" + name + "': 'radius' must be >= 0");
    const std::string profile = def.get_or<std::string>("profile", "tophat");
    gaussian = (profile == "gaussian");
    if (!gaussian && (profile != "tophat"))
        Abort("Action '" + name + "': 'profile' must be 'tophat' or 'gaussian'");
    extent = gaussian ? 2.0 * radius : radius;  // gaussian truncated at e^-4
    rays_per_cell = def.get_or("rays_per_cell", 2);
    if (rays_per_cell < 1) Abort("Action '" + name + "': 'rays_per_cell' must be >= 1");

    get_udf(def["intensity"], intensity, 0.0);
    // I [W m^-2] * dA [x_ref^2] -> power in code units (prs_ref x_ref^3 / t_ref)
    I_scale = MFP::t_ref / (MFP::prs_ref * MFP::x_ref);

    // plasma side
    const HydroGas& gas = *state->gas;
    electron_state = true;
    for (const Real q : gas.charge) electron_state = electron_state && (q < 0.0);
    Z_eff = def.get_or("Z_eff", 1.0);

    const Real eps0 = 8.8541878128e-12, m_e = 9.1093837015e-31, e = 1.602176634e-19,
               c = 299792458.0, kB = 1.38064852e-23;
    const Real wavelength = def["wavelength"];
    const Real omega = 2.0 * M_PI * c / wavelength;
    n_c = eps0 * m_e * omega * omega / (e * e) / MFP::n_ref;
    crit_frac = def.get_or("critical_fraction", 1.0);
    if (crit_frac <= 0.0) Abort("Action '" + name + "': 'critical_fraction' must be > 0");

    kT_eV = kB * MFP::T_ref / e;
    n_ref_cm3 = MFP::n_ref * 1.0e-6;
    // kappa_code = C_kappa * lnLambda * Z2n_code * T_code^-3/2 * x / sqrt(1 - x)
    C_kappa = 2.91e-12 * MFP::n_ref * std::pow(kT_eV, -1.5) * MFP::x_ref / c;

    const sol::object lnL = def["lnLambda"];
    lnLambda_nrl = false;
    lnLambda = 0.0;
    if (lnL.is<Real>()) {
        lnLambda = lnL.as<Real>();
    } else if (lnL.is<std::string>() && (lnL.as<std::string>() == "NRL")) {
        lnLambda_nrl = true;
    } else {
        Abort("Action '" + name + "': 'lnLambda' must be a number or 'NRL'");
    }
}

void Laser::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("Laser::get_data");
    Vector<Array<int, 2>> options = {{state->global_idx, 0}};
    Action::get_data(mfp, options, update, time);
}

// transverse lattice of ray offsets in the plane normal to `dir`
void Laser::build_lattice(const Real* dx)
{
    lattice.clear();
    auto weight = [&](const Real r2) {
        if (gaussian) return std::exp(-r2 / (radius * radius));
        return (r2 <= radius * radius) ? 1.0 : 0.0;
    };

#if AMREX_SPACEDIM == 1
    lattice.push_back({{0.0, 0.0, 0.0}, 1.0});  // plane wave
#else
    Real dx_min = dx[0];
    for (int i = 1; i < AMREX_SPACEDIM; ++i) dx_min = std::min(dx_min, dx[i]);
    const Real s = dx_min / rays_per_cell;
    if (extent <= 0.0) {
        lattice.push_back({{0.0, 0.0, 0.0}, 1.0});  // single ray, unit area
        return;
    }
    // midpoint lattice: offsets (a + 1/2) s, so a beam axis on a cell face or
    // centre never puts rays exactly on faces (which would alias with the grid)
    const int n = std::ceil(extent / s);
    #if AMREX_SPACEDIM == 2
    const Array<Real, 3> p = {-dir[1], dir[0], 0.0};
    for (int a = -n; a < n; ++a) {
        const Real u = (a + 0.5) * s, w = weight(u * u);
        if (w > 0.0) lattice.push_back({{u * p[0], u * p[1], 0.0}, w * s});
    }
    #else
    // p1 = dir x (axis least aligned with dir), p2 = dir x p1
    int ax = 0;
    for (int i = 1; i < 3; ++i) if (std::abs(dir[i]) < std::abs(dir[ax])) ax = i;
    Array<Real, 3> a = {0.0, 0.0, 0.0};
    a[ax] = 1.0;
    Array<Real, 3> p1 = {dir[1] * a[2] - dir[2] * a[1], dir[2] * a[0] - dir[0] * a[2], dir[0] * a[1] - dir[1] * a[0]};
    const Real p1n = std::sqrt(p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]);
    for (Real& v : p1) v /= p1n;
    const Array<Real, 3> p2 = {dir[1] * p1[2] - dir[2] * p1[1], dir[2] * p1[0] - dir[0] * p1[2], dir[0] * p1[1] - dir[1] * p1[0]};
    for (int b = -n; b < n; ++b) {
        for (int c = -n; c < n; ++c) {
            const Real u = (b + 0.5) * s, v = (c + 0.5) * s, w = weight(u * u + v * v);
            if (w > 0.0)
                lattice.push_back({{u * p1[0] + v * p2[0], u * p1[1] + v * p2[1], u * p1[2] + v * p2[2]}, w * s * s});
        }
    }
    #endif
#endif
}

void Laser::trace(Array<Real, 3> o,
                  Real P,
                  const Array4<const Real>& f4,
                  const Array4<Real>& d4,
                  const Real* dx,
                  const Real* plo,
                  const Real* phi,
                  const Box& dom,
                  Real& absorbed,
                  Real& dumped,
                  Real& transmitted) const
{
    constexpr Real inf = std::numeric_limits<Real>::max();

    // ray parameter range inside the domain (slab test)
    Real s0 = 0.0, s1 = inf;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        if (dir[i] == 0.0) {
            if ((o[i] < plo[i]) || (o[i] > phi[i])) { transmitted += P; return; }
        } else {
            const Real t1 = (plo[i] - o[i]) / dir[i], t2 = (phi[i] - o[i]) / dir[i];
            s0 = std::max(s0, std::min(t1, t2));
            s1 = std::min(s1, std::max(t1, t2));
        }
    }
    if (s1 <= s0) { transmitted += P; return; }

    // Amanatides-Woo setup at the entry point
    int ijk[3] = {0, 0, 0}, step[3] = {0, 0, 0};
    Real tMax[3] = {inf, inf, inf}, tDelta[3] = {inf, inf, inf};
    Real V = 1.0;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const Real pos = o[i] + s0 * dir[i];
        ijk[i] = std::clamp((int)std::floor((pos - plo[i]) / dx[i]), dom.smallEnd(i), dom.bigEnd(i));
        V *= dx[i];
        if (dir[i] > 0.0) {
            step[i] = 1;
            tMax[i] = s0 + (plo[i] + (ijk[i] + 1) * dx[i] - pos) / dir[i];
            tDelta[i] = dx[i] / dir[i];
        } else if (dir[i] < 0.0) {
            step[i] = -1;
            tMax[i] = s0 + (plo[i] + ijk[i] * dx[i] - pos) / dir[i];
            tDelta[i] = -dx[i] / dir[i];
        }
    }

    const Real P0 = P;
    Real s = s0;
    while (true) {
        int ax = 0;
        for (int i = 1; i < AMREX_SPACEDIM; ++i) if (tMax[i] < tMax[ax]) ax = i;
        const Real s_next = std::min(tMax[ax], s1);
        const Real ds = s_next - s;

        const IntVect iv(AMREX_D_DECL(ijk[0], ijk[1], ijk[2]));
        const Real n_e = f4(iv, 0), Z2n = f4(iv, 1), T = f4(iv, 2);
        const Real x = n_e / n_c;

        if (x >= crit_frac) {  // critical surface: deposit the remainder here
            d4(iv) += P / V;
            dumped += P;
            return;
        }
        if ((ds > 0.0) && (x > 0.0) && (x < 1.0) && (T > 0.0) && (Z2n > 0.0)) {
            Real lnL = lnLambda;
            if (lnLambda_nrl) {  // NRL electron-ion Coulomb logarithm, n_e cm^-3, T_e eV
                const Real ne = n_e * n_ref_cm3, Te = T * kT_eV;
                lnL = (Te < 10.0 * Z_eff * Z_eff) ? 23.0 - std::log(std::sqrt(ne) * Z_eff * std::pow(Te, -1.5))
                                                  : 24.0 - std::log(std::sqrt(ne) / Te);
            }
            const Real kappa = C_kappa * lnL * Z2n * std::pow(T, -1.5) * x / std::sqrt(1.0 - x);
            const Real dP = P * (1.0 - std::exp(-kappa * ds));
            d4(iv) += dP / V;
            P -= dP;
            absorbed += dP;
            if (P < 1.0e-12 * P0) {  // numerically extinguished
                d4(iv) += P / V;
                absorbed += P;
                return;
            }
        }

        s = s_next;
        if (s >= s1) { transmitted += P; return; }
        ijk[ax] += step[ax];
        tMax[ax] += tDelta[ax];
        if ((ijk[ax] < dom.smallEnd(ax)) || (ijk[ax] > dom.bigEnd(ax))) { transmitted += P; return; }
    }
}

void Laser::calc_time_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt)
{
    BL_PROFILE("Laser::calc_time_derivative");

    if (mfp->get_level() > 0) return;

    const Geometry& geom = mfp->Geom();
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    const Real* phi = geom.ProbHi();
    const Box& dom = geom.Domain();

    if (lattice.empty()) build_lattice(dx);

    MultiFab& U = update[state->data_idx].U;
    MultiFab& dU = update[state->data_idx].dU;
    update[state->data_idx].dU_status = UpdateData::Status::Changed;

    // 1. n_e, sum Z^2 n_i, T_e on this rank's boxes
    MultiFab fields(U.boxArray(), U.DistributionMap(), 3, 0);
    const HydroGas& gas = *state->gas;
    Vector<Real> Uc(state->n_cons()), Qc(state->n_prim()), alpha;
    for (MFIter mfi(fields); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const Array4<const Real>& U4 = U.const_array(mfi);
        const Array4<Real>& f4 = fields.array(mfi);
        amrex::LoopOnCpu(box, [&](int i, int j, int k) {
            for (size_t l = 0; l < Uc.size(); ++l) Uc[l] = U4(i, j, k, l);
            Real n_e = 0.0, Z2n = 0.0, T = 0.0;
            const Real rho = Uc[+HydroDef::ConsIdx::Density];
            if (rho > state->effective_zero) {
                gas.cons2prim(Uc, Qc);
                gas.get_alpha_fractions_from_prim(Qc, alpha);
                T = Qc[+HydroDef::PrimIdx::Temp];
                for (size_t n = 0; n < gas.n_species(); ++n) {
                    const Real ni = alpha[n] * rho / gas.mass[n], q = gas.charge[n];
                    if (electron_state) {
                        n_e += ni;
                    } else if (q > 0.0) {
                        n_e += q * ni;  // quasi-neutral
                        Z2n += q * q * ni;
                    }
                }
                if (electron_state) Z2n = Z_eff * n_e;
            }
            f4(i, j, k, 0) = n_e;
            f4(i, j, k, 1) = Z2n;
            f4(i, j, k, 2) = T;
        });
    }

    // 2. gather the level onto rank 0
    const BoxArray ba_all(dom);
    const DistributionMapping dm_all(Vector<int>{0});
    MultiFab all(ba_all, dm_all, 3, 0);
    all.ParallelCopy(fields);
    MultiFab dep_all(ba_all, dm_all, 1, 0);
    dep_all.setVal(0.0);

    // 3. march every ray (only rank 0 owns a box)
    Real P_in = 0.0, absorbed = 0.0, dumped = 0.0, transmitted = 0.0;
    const Real I_now = I_scale * intensity({{"t", time}});
    for (MFIter mfi(all); mfi.isValid(); ++mfi) {
        const Array4<const Real>& f4 = all.const_array(mfi);
        const Array4<Real>& d4 = dep_all.array(mfi);
        for (const auto& ray : lattice) {
            const Real P = I_now * ray.second;
            P_in += P;
            trace({origin[0] + ray.first[0], origin[1] + ray.first[1], origin[2] + ray.first[2]},
                  P, f4, d4, dx, plo, phi, dom, absorbed, dumped, transmitted);
        }
    }

    // 4. scatter the deposited power density back and apply
    MultiFab dep(U.boxArray(), U.DistributionMap(), 1, 0);
    dep.ParallelCopy(dep_all);
    MultiFab::Saxpy(dU, dt, dep, 0, +HydroDef::ConsIdx::Eden, 1, 0);

    if (MFP::verbosity >= 2) {
        amrex::Print().SetPrecision(16)
          << "[laser] t = " << time << " P_in = " << P_in << " absorbed = " << absorbed
          << " dumped = " << dumped << " transmitted = " << transmitted << "\n";
    }
}
