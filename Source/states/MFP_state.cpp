#include "MFP_state.H"

#include <math.h>
#include "MFP_global.H"
#include "MFP_Riemann_solvers.H"
#include "MFP_diagnostics.H"
#include "MFP_viscous.H"

#ifdef AMREX_USE_EB
#include "MFP_fillcc_eb.H"
#endif

#include <AMReX_YAFluxRegister.H>
using CellType = YAFluxRegister::CellType;


using GD = GlobalData;

State::State()
{
    name = "null";
    global_idx = -1;

    shock_idx = -1;

    refine_grad_max_lev = -1;

    reflux = 1;

    project_divergence = 0;

    mass = {0,0};
    charge = {0,0};
    gamma = {0,0};

    num_grow = 0;

    pressure_relaxation_rate = 0.0;

    particle_index = -1;

}
State::~State(){}

void State::init_from_lua()
{
    BL_PROFILE("State::init_from_lua");
    //
    // get the initial flow state
    //
    set_udf();

    //
    // reflux
    //
    set_reflux();

    //
    // get reconstruction
    //

    set_reconstruction();

    //
    // field flux solver
    //

    set_flux();

    //
    // refinement
    //
    set_grad_refinement();

    //
    // number of ghost cells
    //
    set_num_grow();

#ifdef AMREX_USE_EB
    //
    // how to handle boundary divergence
    //
    set_eb_divergence();
#endif

}

void State::post_init_from_lua()
{
    // fix periodicity
    for (BCRec &bcr : boundary_conditions.fill_bc) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            if (GD::periodic[d]) {
                bcr.setHi(d, BCType::int_dir);
                bcr.setLo(d, BCType::int_dir);
            }
        }
    }

    // update viscosity, this is required as some of them need to have fully defined source terms
    // so that they know what states are linked to it via sources
    if (viscous)
        viscous->update_linked_states();

}

void State::set_udf()
{
    sol::state& lua = GD::lua;

    bool success;

    sol::table state_def = lua["states"][name];

    // check if we have 'value' defined
    sol::table value = state_def["value"].get_or(sol::table());

    if ((!value.valid() || value.empty()) && ParallelDescriptor::IOProcessor())
        Warning("WARNING: State '"+name+"' does not have 'value' defined for initial conditions, using defaults");

    // get a list of any initialisation functions that need to be called during the run

    sol::table dynamic = state_def["dynamic"].get_or(sol::table());

    const auto ignore = init_with_value();

    const Vector<std::string> &prim_names = get_prim_names();
    for (int i = 0; i<prim_names.size(); ++i) {

        const std::string &comp = prim_names[i];

        Optional3D1VFunction v;

        if (!value.valid() || value.empty()) {
            v.set_value(0.0);
            success = false;
        } else {
            success = get_udf(value[comp], v, 0.0);
        }

        // set default values
        if (!success) {
            for (const auto &j : ignore) {
                if (i == j.first) {
                    v.set_value(j.second);
                    break;
                }
            }
        }

        functions[comp] = v;

        if (dynamic.valid()) {
            for (const auto &d : dynamic) {
                if (d.second.as<std::string>().compare(comp) == 0) {
                    dynamic_functions[i] = &functions[comp];
                }
            }
        }
    }

    return;
}

void State::set_reconstruction()
{
    PhysicsFactory<Reconstruction> rfact = GetReconstructionFactory();

    sol::table state_def = GD::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string rec = state_def["reconstruction"].get_or<std::string>("null");

    state_def["reconstruction"] = rec; // consistency when using default option "null"

    if (is_transported() && rec == "null")
        Abort("Reconstruction option required for state '"+name+"'. Options are "+vec2str(rfact.getKeys()));

    reconstruction = rfact.Build(rec, state_def);

    if (!reconstruction)
        Abort("Invalid reconstruction option '"+rec+"'. Options are "+vec2str(rfact.getKeys()));

    return;

}

void State::set_flux()
{

    if (!is_transported()) return;

    sol::state& lua = GD::lua;

    PhysicsFactory<RiemannSolver> rfact = GetRiemannSolverFactory();

    sol::table state_def = GD::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string flux = state_def["flux"].get_or<std::string>("null");

    if (flux == "null")
        Abort("Flux option required for state '"+name+"'. Options are "+vec2str(rfact.getKeys()));

    flux_solver = rfact.Build(flux, state_def);

    if (!flux_solver)
        Abort("Invalid flux solver option '"+flux+"'. Options are "+vec2str(rfact.getKeys()));

    return;

}

void State::set_reflux()
{
    if (AMREX_SPACEDIM < 2) {
        reflux = 0;
    } else {
        reflux = GD::lua["states"][name]["reflux"].get_or(1);
    }
}

void State::set_viscosity()
{

    //
    // viscous terms coefficients
    //

    PhysicsFactory<Viscous> vfact = GetViscousFactory();

    sol::table state_def = GD::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string visc = state_def["viscosity"]["type"].get_or<std::string>("");

    viscous = vfact.Build(visc, state_def);

    if (!visc.empty() && !viscous)
        Abort("Invalid viscosity option '"+visc+"'. Options are "+vec2str(vfact.getKeys()));

}

void State::set_shock_detector()
{

    PhysicsFactory<ShockDetector> sdfact = GetShockDetectorFactory();

    sol::table sd_def = GD::lua["states"][name]["shock_detector"].get_or(sol::table());

    if (!sd_def.valid()) return;

    sd_def["global_idx"] = global_idx;

    std::string sd_name = sd_def["name"].get_or<std::string>("");

    shock_detector = sdfact.Build(sd_name, sd_def);

    if (!sd_name.empty() && !shock_detector)
        Abort("Invalid shock_detector option '"+sd_name+"'. Options are "+vec2str(sdfact.getKeys()));

    shock_idx = 0;
}


void State::set_grad_refinement()
{

    // grab the level at which to stop refining
    refine_grad_max_lev = GD::lua["states"][name]["refine_max_level"].get_or(-1);


    // get the threshold gradients for any of the specified primitives
    sol::table thresh = GD::lua["states"][name]["refine_grad_threshold"].get_or(sol::table());

    if (!thresh.valid()) {
        return;
    }

    const Vector<std::string> &cons = get_cons_names();
    const Vector<std::string> &prim = get_prim_names();

    Real min_value = thresh["min_value"].get_or(0.0);

    std::pair<bool, int> found;
    for (const auto& t : thresh) {

        std::string tag_name = t.first.as<std::string>();
        Real value = t.second.as<Real>();

        if (tag_name.compare("min_value")==0)
            continue;

        // look through conserved
        found = findInVector(cons, tag_name);
        if (found.first) {
            refine_grad_threshold.push_back({false, found.second, value, min_value});
        } else {
            // look through primitives
            found = findInVector(prim, tag_name);
            if (found.first) {
                refine_grad_threshold.push_back({true, found.second, value, min_value});
            } else {
                Abort("refine threshold selection "+tag_name+" is not a valid selection. Choose one of " + vec2str(cons)+" (preferred) or "+vec2str(prim));
            }
        }
    }
}

//void State::set_num_grow()
void State::set_num_grow(int n)
{
    num_grow = std::max(num_grow, n);
    num_grow = std::max(num_grow, reconstruction->get_num_grow());
}

void State::get_state_vector(const FArrayBox& U, const int i, const int j, const int k, Vector<Real>& S)
{
    BL_PROFILE("State::get_state_vector");

    Array4<const Real> const& U4 = U.array();

    for (int n=0; n<U.nComp(); ++n) {
        S[n] = U4(i,j,k,n);
    }
}

void State::get_state_values(const Box& box,
                             const FArrayBox& src,
                             std::map<std::string,FArrayBox>& out,
                             Vector<std::string>& updated
                             EB_OPTIONAL(,const FArrayBox& vfrac)
                             ) const
{
    return;
}

void State::calc_velocity(const Box& box,
                          FArrayBox& src,
                          FArrayBox& vel
                          EB_OPTIONAL(,const EBCellFlagFab& flag)
                          ) const
{
    return;
}

void State::calc_primitives(const Box& box,
                            FArrayBox& cons,
                            FArrayBox& prim,
                            const Real* dx,
                            const Real t,
                            const Real* prob_lo
                            EB_OPTIONAL(,const FArrayBox& vfrac)
                            ) const
{
    BL_PROFILE("State::calc_primitives");

    Vector<Real> U(n_cons()), Q(n_prim());

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);
    Array4<Real> const& s4 = cons.array();
    Array4<Real> const& p4 = prim.array();

#ifdef AMREX_USE_EB
    Array4<const Real> const& vfrac4 = vfrac.array();

    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);
#endif

    Real x, y, z;
    const int nc = n_cons();
    const int np = n_prim();

    for     (int k = lo.z; k <= hi.z; ++k) {
        z = prob_lo[2] + (k + 0.5)*dx[2];
        for   (int j = lo.y; j <= hi.y; ++j) {
            y = prob_lo[1] + (j + 0.5)*dx[1];
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                x = prob_lo[0] + (i + 0.5)*dx[0];

#ifdef AMREX_USE_EB
                if (vfrac4(i,j,k) == 0.0) {

                    // iterate over all neighbouring cells checking if it has valid data
                    // from these calculate the volume weighted average to populate the
                    // covered cell
                    Real vtot = 0.0;

                    std::fill(U.begin(), U.end(), 0.0);
                    for (const auto& index : grab) {

                        const int ii = i+index[0];
                        const int jj = j+index[1];
                        const int kk = k+index[2];

                        // make sure our stencil is within bounds
                        if ((lo.x > ii) || (ii > hi.x) ||
                                (lo.y > jj) || (jj > hi.y) ||
                                (lo.z > kk) || (kk > hi.z)) continue;

                        const Real vf = vfrac4(ii,jj,kk);
                        if (vf > 0.0) {
                            for (int n=0; n<nc; ++n) {
                                U[n] += vf*s4(ii,jj,kk,n);
                            }

                            vtot += vf;
                        }
                    }

                    // if we were close enough to a boundary to have valid data we
                    // average out the volume fraction weighted contributions, otherwise,
                    // fill in the primitives with zeros
                    if (vtot > 0.0) {
                        for (int n=0; n<nc; ++n) {
                            U[n] /= vtot;
                        }
                    } else {
                        for (int n=0; n<np; ++n) {
                            p4(i,j,k,n) = 0.0;
                        }
                        continue;
                    }
                } else {
#endif
                    // grab the conserved variables
                    for (int n=0; n<nc; ++n) {
                        U[n] = s4(i,j,k,n);
                    }
#ifdef AMREX_USE_EB
                }
#endif


                // convert to primitive
                cons2prim(U, Q);

                // modify the primitives vector if needed and upload back to
                // the conserved vector
                if (!dynamic_functions.empty()) {
                    for (const auto &f : dynamic_functions) {
                        Q[f.first] = (*f.second)(x, y, z, t);//TODO uodate this for new commit, which passes in a
                    }

                    // copy into primitive
                    for (int n=0; n<np; ++n) {
                        p4(i,j,k,n) = Q[n];
                    }

                    // convert primitive to conserved
                    prim2cons(Q, U);

                    // copy back into conserved array
                    for (int n=0; n<nc; ++n) {
                        s4(i,j,k,n) = U[n];
                    }
                }

                // copy into primitive
                for (int n=0; n<np; ++n) {
                    p4(i,j,k,n) = Q[n];
                }
            }
        }
    }

    return;
}



void State::calc_reconstruction(const Box& box,
                                FArrayBox &prim,
                                Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                EB_OPTIONAL(,const EBCellFlagFab& flag)
                                EB_OPTIONAL(,const FArrayBox &vfrac)
                                ) const
{
    BL_PROFILE("State::calc_reconstruction");

    // don't need to do this iota, could just use a loop, but useful for
    // when we don't need to do all components
    Vector<int> index(prim.nComp());
    std::iota (index.begin(), index.end(), 0);

    int N = index.size();

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& src4 = prim.array();

#ifdef AMREX_USE_EB

    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);


    Array4<const EBCellFlag> const& f4 = flag.array();
    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    Vector<Real> stencil(reconstruction->stencil_length);
    int offset = reconstruction->stencil_length/2;
    Array<int,3> stencil_index;
    Real lo_face, hi_face;


    // cycle over dimensions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        // make sure our arrays for putting lo and hi reconstructed values into
        // are the corect size
        rlo[d].resize(box, N);
        rhi[d].resize(box, N);

#ifdef AMREX_USE_EB
        if (check_eb) {
            rlo[d].copy(prim,box);
            rhi[d].copy(prim,box);
        }
#endif

        Array4<Real> const& lo4 = rlo[d].array();
        Array4<Real> const& hi4 = rhi[d].array();

        // cycle over all components
        for (const int& n : index) {

            for     (int k = lo.z; k <= hi.z; ++k) {
                for   (int j = lo.y; j <= hi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                            for (int i = lo.x; i <= hi.x; ++i) {



#ifdef AMREX_USE_EB
                        if (check_eb) {
                            // covered cell doesn't need calculating
                            if (f4(i,j,k).isCovered() || check_covered_stencil(f4, i, j, k, d, reconstruction->stencil_length)) {
                                continue;
                            }
                        }
#endif
                        stencil_index.fill(0);
                        for (int s=0; s<reconstruction->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            stencil[s] = src4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], n);
                        }

                        // perform reconstruction
                        reconstruction->get_face_values(stencil, lo_face, hi_face);

                        lo4(i,j,k,n) = lo_face;
                        hi4(i,j,k,n) = hi_face;
                    }
                }
            }
        }
    }


    return;
}


void State::calc_slope(const Box& box,
                       const FArrayBox& src,
                       FArrayBox &slope,
                       EB_OPTIONAL(const EBCellFlagFab &flag,)
                       const Real *dx,
                       int index,
                       int dim,
                       Reconstruction &reco)
{
    BL_PROFILE("State::calc_slope");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& src4 = src.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
    bool check_eb = flag.getType() != FabType::regular;

#endif

    Vector<Real> stencil(reco.stencil_length);
    int offset = reco.stencil_length/2;
    Array<int,3> stencil_index;



    // make sure our array is the corect size
    slope.resize(box);

    Array4<Real> const& s4 = slope.array();



    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {


                stencil_index.fill(0);

#ifdef AMREX_USE_EB
                if (check_eb) {
                    // covered cell doesn't need calculating
                    if (f4(i,j,k).isCovered() || check_covered_stencil(f4, i, j, k, dim, reco.stencil_length)) {
                        s4(i,j,k) = 0.0;
                        continue;
                    }
                }
#endif

                for (int s=0; s<reco.stencil_length; ++s) {
                    stencil_index[dim] = s - offset;
                    stencil[s] = src4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], index);
                }

                // perform reconstruction
                s4(i,j,k) = reco.get_slope(stencil)/dx[dim]; // account for local cell size
            }
        }
    }


    return;
}

void State::calc_time_averaged_faces(const Box& box,
                                     const FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                                     EB_OPTIONAL(const EBCellFlagFab& flag,)
                                     const Real* dx,
                                     Real dt) const
{
    BL_PROFILE("State::calc_time_averaged_faces");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& p4 = prim.array();

    const int np = prim.nComp();

    Vector<Real> Q(np);

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    Real lo_face, centre, hi_face, a6;
    Real dt_2dx;

    const Real ft = 4.0/3.0;

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    continue;
                }
#endif

                for (int n=0; n<np; ++n) {
                    Q[n] = p4(i,j,k,n);
                }

                const RealArray local_c = get_speed_from_prim(Q);


                for (int d=0; d<AMREX_SPACEDIM; ++d) {

                    Array4<Real> const& lo4 = rlo[d].array();
                    Array4<Real> const& hi4 = rhi[d].array();

                    dt_2dx = 0.5*dt/dx[d];

                    // cycle over all components
                    for (int n=0; n<np; ++n) {

                        lo_face = lo4(i,j,k,n);
                        hi_face = hi4(i,j,k,n);
                        centre = p4(i,j,k,n);

                        a6 = 6.0*centre - 3*(lo_face + hi_face);

                        lo4(i,j,k,n) += local_c[d]*dt_2dx*(hi_face - lo_face + (1 - local_c[d]*ft*dt_2dx)*a6);
                        hi4(i,j,k,n) -= local_c[d]*dt_2dx*(hi_face - lo_face - (1 - local_c[d]*ft*dt_2dx)*a6);

                    }
                }
            }
        }
    }

    return;
}

Vector<BoundaryInfo> State::get_bc_limits(const Box& box,
                                          const Geometry &geom) const
{
    BL_PROFILE("State::get_bc_limits");

    Vector<BoundaryInfo> limits;

    if (AMREX_D_TERM(geom.isPeriodic(0), && geom.isPeriodic(1), && geom.isPeriodic(2)))
        return limits;


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
            BoundaryInfo &b = limits.back();
            b.imin = lo.x;
            b.imax = ilo - 1;
            b.dir = 0;
            b.lo_hi = 0;
        }

        if (hi.x > ihi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.imin = ihi+1;
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
            BoundaryInfo &b = limits.back();
            b.jmin = lo.y;
            b.jmax = jlo - 1;
            b.dir = 1;
            b.lo_hi = 0;
        }

        if (hi.y > jhi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.jmin = jhi+1;
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
            BoundaryInfo &b = limits.back();
            b.kmin = lo.z;
            b.kmax = klo - 1;
            b.dir = 2;
            b.lo_hi = 0;
        }

        if (hi.z > khi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.kmin = khi+1;
            b.kmax = hi.z;
            b.dir = 2;
            b.lo_hi = 1;
        }
    }
#endif

    return limits;

}

#ifdef AMREX_USE_EB

void State::update_eb_vfrac(const Geometry &geom,
                      FArrayBox& vfrac) const
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
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = problo[i] + dx[i]*(lo[i]-dom_lo[i]);
    }

//    plot_FAB_2d(vfrac, 0, "vfrac before", false, false);

    Vector<BoundaryInfo> limits = get_bc_limits(vfrac.box(), geom);
    if (!limits.empty()) {
        Array4<Real> const& vfrac4 = vfrac.array();

        // if the cell we are copying from has no solid in it, don't copy
        const Real check = 1.0;

        for (size_t i=0; i<limits.size(); ++i) {
            const auto L = limits[i];

            const IntVect box_lo(AMREX_D_DECL(L.imin,L.jmin,L.kmin));
            const IntVect box_hi(AMREX_D_DECL(L.imax,L.jmax,L.kmax));
            const Box box(box_lo, box_hi);
            fab_filcc_eb(box, vfrac4, 1, geom.Domain(), dx, xlo, &bc, check);
        }
    }
//        plot_FAB_2d(vfrac, 0, "vfrac after", false, true);

}

void State::update_eb_flags(const Geometry &geom,
                      EBCellFlagFab& flags) const
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
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = problo[i] + dx[i]*(lo[i]-dom_lo[i]);
    }

    Vector<BoundaryInfo> limits = get_bc_limits(flags.box(), geom);
    if (!limits.empty()) {

        Array4<EBCellFlag> const& flags4 = flags.array();

        // if the cell we are copying is regular then don't copy
        EBCellFlag check; check.setRegular();

        for (size_t i=0; i<limits.size(); ++i) {
            const auto L = limits[i];

            const IntVect box_lo(AMREX_D_DECL(L.imin,L.jmin,L.kmin));
            const IntVect box_hi(AMREX_D_DECL(L.imax,L.jmax,L.kmax));
            const Box box(box_lo, box_hi);
            fab_filcc_eb(box, flags4, 1, geom.Domain(), dx, xlo, &bc, check);
        }

    }

}


#endif

void State::update_boundary_cells(const Box& box,
                                  const Geometry &geom,
                                  FArrayBox &prim,
                                  EB_OPTIONAL(const FArrayBox& vfrac,)
                                  const Real time) const
{
    BL_PROFILE("State::update_boundary_cells");
    Vector<BoundaryInfo> limits = get_bc_limits(box, geom);

    if (limits.empty())
        return;

    const int ncomp = prim.nComp();
    const Box& domain = geom.Domain();
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    const Vector<std::string> prims = get_prim_names();

    Real x, y, z;
    std::map<std::string, Real> Q{{"t", time}};

    Array<int,3> grab;

    Array4<Real> const& p4 = prim.array();

#ifdef AMREX_USE_EB
    Array4<const Real> const& vfrac4 = vfrac.array();
#endif

    const BCRec &bc = boundary_conditions.where_is_inflow;
    for (const auto &L : limits) {
        if (((L.lo_hi == 0) && (bc.lo(L.dir) == BCType::ext_dir))||
                ((L.lo_hi == 1) && (bc.hi(L.dir) == BCType::ext_dir))) {

            if (L.lo_hi == 0) {
                grab[L.dir] = domain.smallEnd(L.dir);
            } else {
                grab[L.dir] = domain.bigEnd(L.dir);
            }

            for (int k=L.kmin; k<=L.kmax; ++k) {
                z = prob_lo[2] + (k + 0.5)*dx[2];
                Q["z"] = z;
                if (L.dir != 2) {
                    grab[2] = k;
                }
                for (int j=L.jmin; j<=L.jmax; ++j) {
                    y = prob_lo[1] + (j + 0.5)*dx[1];
                    Q["y"] = y;
                    if (L.dir != 1) {
                        grab[1] = j;
                    }
                    for (int i=L.imin; i<=L.imax; ++i) {
                        x = prob_lo[0] + (i + 0.5)*dx[0];
                        Q["x"] = x;
                        if (L.dir != 0) {
                            grab[0] = i;
                        }
#ifdef AMREX_USE_EB
                        if (vfrac4(grab[0],grab[1],grab[2]) == 0.0) {
                            continue;
                        }
#endif

                        // get data from closest internal cell
                        for (int n=0; n<ncomp; ++n) {
                            Q[prims[n]] = p4(grab[0],grab[1],grab[2],n);
                        }

                        // update the primitives from our UDFs, but only those that are valid functions
                        for (int n=0; n<ncomp; ++n) {
                            const Optional3D1VFunction &f = boundary_conditions.get(L.lo_hi, L.dir, prims[n]);
                            if (f.is_valid()) {
                                p4(i,j,k,n) = f(Q);
                            }
                        }
                    }
                }
            }
        }
    }
}

/*
* the data passed to this function is indexed by cell
* but resides at the interface
*/

void State::face_bc(const int dir,
                    Box const& box,
                    const FArrayBox& src,
                    FArrayBox& dest,
                    const int numcomp,
                    const Geometry &geom,
                    EB_OPTIONAL(const EBCellFlagFab& flag,)
                    const Real time,
                    const bool do_all) const
{
    BL_PROFILE("State::face_bc");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    const Box& domain = geom.Domain();
    const Dim3 domlo = amrex::lbound(domain);
    const Dim3 domhi = amrex::ubound(domain);

    const Real* prob_lo = geom.ProbLo();
    const Real* dx = geom.CellSize();

    Real x, y, z;
    std::map<std::string, Real> Q{{"t", time}};

    const Vector<std::string> prims = get_prim_names();

    // define the limits of our operations
    Vector<BoundaryInfo> limits;

    BoundaryInfo bi;
    bi.imin = lo.x;
    bi.imax = hi.x;
    bi.jmin = lo.y;
    bi.jmax = hi.y;
    bi.kmin = lo.z;
    bi.kmax = hi.z;



    if (dir == 0) {

        int ilo = domlo.x;
        int ihi = domhi.x + 1; // account for face indexing

        if (lo.x == ilo) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.imin = lo.x;
            b.imax = lo.x;
            b.lo_hi = 0;
        }

        if (hi.x == ihi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.imin = hi.x;
            b.imax = hi.x;
            b.lo_hi = 1;
        }
    }

#if AMREX_SPACEDIM >= 2
    if (dir == 1) {
        int jlo = domlo.y;
        int jhi = domhi.y + 1; // account for face indexing

        if (lo.y == jlo) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.jmin = lo.y;
            b.jmax = lo.y;
            b.lo_hi = 0;
        }

        if (hi.y == jhi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.jmin = hi.y;
            b.jmax = hi.y;
            b.lo_hi = 1;
        }
    }
#endif

#if AMREX_SPACEDIM == 3
    if (dir == 2) {
        int klo = domlo.z;
        int khi = domhi.z + 1; // account for face indexing

        if (lo.z == klo) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.kmin = lo.z;
            b.kmax = lo.z;
            b.lo_hi = 0;
        }

        if (hi.z == khi) {
            limits.push_back(bi);
            BoundaryInfo &b = limits.back();
            b.kmin = hi.z;
            b.kmax = hi.z;
            b.lo_hi = 1;
        }
    }
#endif

    Array4<Real> const& d4 = dest.array();
    Array4<const Real> const& s4 = src.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    if (do_all) {
        // first do the usual boundary conditions
        for (int n=0; n<numcomp; ++n) {
            const BCRec &bc = boundary_conditions.fill_bc[n];

            for (const auto &L : limits) {

                if (
                        ((L.lo_hi == 0) && (bc.lo(dir) == BCType::foextrap))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::foextrap))||

                        ((L.lo_hi == 0) && (bc.lo(dir) == BCType::hoextrap))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::hoextrap))||

                        ((L.lo_hi == 0) && (bc.lo(dir) == BCType::reflect_even))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::reflect_even))
                        ) {
                    for (int k=L.kmin; k<=L.kmax; ++k) {
                        for (int j=L.jmin; j<=L.jmax; ++j) {
                            for (int i=L.imin; i<=L.imax; ++i) {

#ifdef AMREX_USE_EB
                                if (f4(i,j,k).isCovered()) {
                                    continue;
                                }
#endif


                                d4(i,j,k,n) = s4(i,j,k,n);
                            }
                        }
                    }
                }
                if (((L.lo_hi == 0) && (bc.lo(dir) == BCType::reflect_odd))||
                        ((L.lo_hi == 1) && (bc.hi(dir) == BCType::reflect_odd))) {
                    for (int k=L.kmin; k<=L.kmax; ++k) {
                        for (int j=L.jmin; j<=L.jmax; ++j) {
                            for (int i=L.imin; i<=L.imax; ++i) {

#ifdef AMREX_USE_EB
                                if (f4(i,j,k).isCovered()) {
                                    continue;
                                }
#endif


                                d4(i,j,k,n) = -s4(i,j,k,n);
                            }
                        }
                    }
                }
            }
        }
    }

    // now do any dirichlet boundaries
    Array<Real,3> offset = {0.5, 0.5, 0.5};
    offset[dir] = 0.0;
    const BCRec &bc = boundary_conditions.where_is_inflow;
    for (const auto &L : limits) {
        if (((L.lo_hi == 0) && (bc.lo(dir) == BCType::ext_dir))||
                ((L.lo_hi == 1) && (bc.hi(dir) == BCType::ext_dir))) {

            for (int k=L.kmin; k<=L.kmax; ++k) {
                z = prob_lo[2] + (k+offset[2])*dx[2];
                Q["z"] = z;
                for (int j=L.jmin; j<=L.jmax; ++j) {
                    y = prob_lo[1] + (j+offset[1])*dx[1];
                    Q["y"] = y;
                    for (int i=L.imin; i<=L.imax; ++i) {
                        x = prob_lo[0] + (i+offset[0])*dx[0];
                        Q["x"] = x;

#ifdef AMREX_USE_EB
                        if (f4(i,j,k).isCovered()) {
                            continue;
                        }
#endif

                        // load data from the other side of the face
                        for (int n=0; n<numcomp; ++n) {
                            Q[prims[n]] = s4(i,j,k,n);
                        }

                        // update the primitives from our UDFs, but only those that are valid functions
                        for (int n=0; n<numcomp; ++n) {
                            const Optional3D1VFunction &f = boundary_conditions.get(L.lo_hi, dir, prims[n]);
                            if (f.is_valid()) {
                                d4(i,j,k,n) = f(Q);
                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

void State::update_face_prim(const Box& box, const Geometry& geom,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                             Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                             EB_OPTIONAL(const EBCellFlagFab& flag,)
                             const Real time, const bool do_all) const
{
    BL_PROFILE("State::update_face_prim");
    const Box domain = geom.Domain();
    const int*     domainlo    = domain.loVect();
    const int*     domainhi    = domain.hiVect();



    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (geom.isPeriodic(d)) {
            continue;
        }

        if (box.smallEnd(d) <= domainlo[d]) {
            Box lo = box;
            lo.setBig(d,  domainlo[d]);
            lo.setSmall(d,domainlo[d]);

            // get the left and right reconstructed primitives for all faces
            FArrayBox &L = r_hi[d];
            FArrayBox &R = r_lo[d];

            // adjust the boxes so that we can use the same index
            L.shift(d,1);

            face_bc(d,
                    lo,
                    R,
                    L,
                    L.nComp(),
                    geom,
                    EB_OPTIONAL(flag,)
                    time,
                    do_all);

            // change it back
            L.shift(d,-1);
        }

        if (box.bigEnd(d) >= domainhi[d]) {
            Box hi = box;
            hi.setBig(d,  domainhi[d]+1);
            hi.setSmall(d,domainhi[d]+1);

            // get the left and right reconstructed primitives for all faces
            FArrayBox &L = r_hi[d];
            FArrayBox &R = r_lo[d];

            // adjust the boxes so that we can use the same index
            L.shift(d,1);

            face_bc(d,
                    hi,
                    L,
                    R,
                    R.nComp(),
                    geom,
                    EB_OPTIONAL(flag,)
                    time,
                    do_all);

            // change it back
            L.shift(d,-1);
        }
    }

}

// given all of the available face values load the ones expected by the flux calc into a vector
void State::load_state_for_flux(Vector<Array4<const Real>> &face,
                                        int i, int j, int k, Vector<Real> &S) const
{
    BL_PROFILE("State::load_state_for_flux");
    const int np = n_prim();

    Array4<const Real> const &f4 = face[global_idx];

    for (int n=0; n<np; ++n ) {
        S[n] = f4(i,j,k,n);
    }
}

void State::calc_fluxes(const Box& box,
                        Vector<Array<FArrayBox, AMREX_SPACEDIM>> &r_lo,
                        Vector<Array<FArrayBox, AMREX_SPACEDIM>> &r_hi,
                        Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                        EB_OPTIONAL(const EBCellFlagFab& flag,)
                        const Real *dx,
                        const Real dt,
                        FArrayBox *shock) const
{
    BL_PROFILE("State::calc_fluxes");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    int ns = r_lo.size();
    int nc = n_cons();
    int np = n_prim();
    int nf = n_flux();
    Vector<int> rotate_idx_flux = get_flux_vector_idx();
    Vector<int> rotate_idx_cons = get_cons_vector_idx();

    Array<int, 3> index;
    Vector<Real> L(nf), R(nf), F(nc);

    // stuff for shock detection
    Real shk;

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
#endif

    // cycle over dimensions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        FArrayBox& flux = fluxes[d];
        flux.setVal(0.0);
        Array4<Real> const& flux4 = flux.array();

        const Box &fbox = flux.box();
        const Dim3 flo = amrex::lbound(fbox);
        const Dim3 fhi = amrex::ubound(fbox);

        index.fill(0);
        index[d] = 1;

        Vector<Array4<const Real>> lo4(ns), hi4(ns);

        for (int nsi=0; nsi<ns; ++nsi) {
            lo4[nsi] = r_lo[nsi][d].array();
            hi4[nsi] = r_hi[nsi][d].array();
        }

        for     (int k = flo.z; k <= fhi.z; ++k) {
            for   (int j = flo.y; j <= fhi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = flo.x; i <= fhi.x; ++i) {

#ifdef AMREX_USE_EB

                    // both cells have to be covered before we skip
                    if (f4(i,j,k).isCovered() && f4(i-index[0],j-index[1],k-index[2]).isCovered())
                        continue;
#endif

                    // get left and right states
                    load_state_for_flux(hi4, i-index[0],j-index[1],k-index[2], L);
                    load_state_for_flux(lo4, i,j,k, R);


                    // rotate the vectors
                    transform_global2local(L, d, rotate_idx_flux);
                    transform_global2local(R, d, rotate_idx_flux);

                    // only compute shock detector when necessary
                    if (shock_detector) {
                        shk = shock_detector->solve(L, R);
                        Array4<Real> const& shock4 = shock->array();
                        shock4(i-index[0],j-index[1],k-index[2]) = std::max(shock4(i-index[0],j-index[1],k-index[2]), shk);
                        shock4(i,j,k) = std::max(shock4(i,j,k), shk);
                    }

                    flux_solver->solve(L, R, F, &shk);

                    // rotate the flux back to local frame
                    transform_local2global(F, d, rotate_idx_cons);

                    // load the flux into the array
                    for (int n=0; n<nc; ++n) {
                        flux4(i,j,k,n) += F[n];
                    }

                }
            }
        }
    }

    return;
}

void State::correct_face_prim(const Box& box,
                              Array<FArrayBox, AMREX_SPACEDIM> &r_lo,
                              Array<FArrayBox, AMREX_SPACEDIM> &r_hi,
                              const Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                              EB_OPTIONAL(const EBCellFlagFab &flag,)
                              const Real *dx,
                              const Real dt) const
{
    BL_PROFILE("State::correct_face_prim");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    int nc = n_cons();
    int np = n_prim();
    //    Vector<int> rotate_idx_prim = get_prim_vector_idx();
    //    Vector<int> rotate_idx_cons = get_cons_vector_idx();

    Vector<Real> L_prim(np), R_prim(np), F(nc);
    Vector<Real> L_cons(np), R_cons(np);

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();

    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    // the below loops selects the alternative dimensions to calculate the corrections from
    // e.g. 1D: d1 = 0, d2 = -
    // e.g. 2D: d1 = 0, d2 = 1
    //          d1 = 1, d2 = 0
    // e.g. 3D: d1 = 0, d2 = 1, 2
    //          d1 = 1, d2 = 2, 0
    //          d1 = 2, d2 = 0, 1

    Array<int, 3> idx1, idx2;

    Real hiF, loF;
    int d1, dd, d2;

    for (d1=0; d1<AMREX_SPACEDIM; ++d1) { // face direction
        idx1.fill(0);
        idx1[d1] = 1;

        Array4<Real> const &lo4 = r_lo[d1].array();
        Array4<Real> const &hi4 = r_hi[d1].array();

        for (dd=1; dd<AMREX_SPACEDIM; ++dd) {
            d2 = (d1+dd)%AMREX_SPACEDIM; // flux direction

            idx2.fill(0);
            idx2[d2] = 1;

            Array4<const Real> const &flux4 = fluxes[d2].array();

            for     (int k = lo.z; k <= hi.z + idx1[2]; ++k) {
                for   (int j = lo.y; j <= hi.y + idx1[1]; ++j) {
                    AMREX_PRAGMA_SIMD
                            for (int i = lo.x; i <= hi.x + idx1[0]; ++i) {

#ifdef AMREX_USE_EB
                        if (check_eb) {
                            // only calculate corrections for faces where the stencil of fluxes is full

                            // a) this interface has a valid flux
                            if (f4(i,j,k).isDisconnected(-idx1[0],-idx1[1],-idx1[2]))
                                continue;

                            // b) right cell has valid fluxes hi and lo
                            if (f4(i,j,k).isDisconnected( idx2[0], idx2[1], idx2[2]))
                                continue;
                            if (f4(i,j,k).isDisconnected(-idx2[0],-idx2[1],-idx2[2]))
                                continue;

                            // c) left cell has valid fluxes hi and lo
                            if (f4(i-idx1[0],j-idx1[1],k-idx1[2]).isDisconnected( idx2[0], idx2[1], idx2[2]))
                                continue;
                            if (f4(i-idx1[0],j-idx1[1],k-idx1[2]).isDisconnected(-idx2[0],-idx2[1],-idx2[2]))
                                continue;
                        }
#endif

                        // get left and right states
                        for (int n=0; n<np; ++n ) {
                            L_prim[n] = hi4(i-idx1[0],j-idx1[1],k-idx1[2],n);
                            R_prim[n] = lo4(i,j,k,n);
                        }

                        // convert to conserved
                        prim2cons(L_prim, L_cons);
                        prim2cons(R_prim, R_cons);

                        /* example for x-direction in 2D
                                   * L & R = reconstructed x- face values
                                   * hiF & loF = fluxes in y- direction on hi and lo faces
                                   *   relative to the reconstructed values
                                   * diagonal lines represent the contributing factors to L & R
                                   *
                                   * L* = L + 0.5*dt*dy*(loF_L - hiF_L)
                                   * R* = R + 0.5*dt*dy*(loF_R - hiF_R)
                                   *
                                   *
                                   * | - - hiF_L - - | - - hiF_R - - |
                                   * |        \      |     /         |
                                   * |          \    |   /           |
                                   * |             L | R             |
                                   * |          /    |   \           |
                                   * |        /      |     \         |
                                   * | - - loF_L - - | - - loF_R - - |
                                   */

                        // left side
                        for (int n=0; n<nc; ++n ) {
                            loF = flux4(i-idx1[0],j-idx1[1],k-idx1[2],n);
                            hiF = flux4(i-idx1[0]+idx2[0],j-idx1[1]+idx2[1],k-idx1[2]+idx2[2],n);

                            L_cons[n] += 0.5*dt/dx[d2]*(loF - hiF);
                        }

                        // right side
                        for (int n=0; n<nc; ++n ) {
                            loF = flux4(i,j,k,n);
                            hiF = flux4(i+idx2[0],j+idx2[1],k+idx2[2],n);

                            R_cons[n] += 0.5*dt/dx[d2]*(loF - hiF);
                        }

                        // convert to primitive
                        cons2prim(L_cons, L_prim);
                        cons2prim(R_cons, R_prim);

                        // update reconstruction
                        for (int n=0; n<np; ++n ) {
                            hi4(i-idx1[0],j-idx1[1],k-idx1[2],n) = L_prim[n];
                            lo4(i,j,k,n) = R_prim[n];
                        }
                    }
                }
            }
        }
    }
}

void State::calc_divergence(const Box& box,
                            Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                            FArrayBox& du,
                            const Real *dx,
                            const Real dt) const
{
    BL_PROFILE("State::calc_divergence");
    // make sure du is empty
    du.setVal(0.0);

    int N = du.nComp();

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<int,3> index = {0,0,0};
    Array4<Real> const& du4 = du.array();
    Real dxinv;

    Array<Array4<const Real>,AMREX_SPACEDIM> flux4;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        flux4[d] = fluxes[d].array();
    }

    for (int n=0; n<N; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

                    // cycle over dimensions
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        index[d] = 1;
                        dxinv = dt/dx[d];

                        const Real& flux_lo  = flux4[d](i,j,k,n);

                        const Real& flux_hi  = flux4[d](i+index[0], j+index[1], k+index[2], n);

                        du4(i,j,k,n) += dxinv*(flux_lo - flux_hi);


                        index[d] = 0;
                    }
                }
            }
        }
    }

    return;
}


// note that this function is adding to du!
void State::calc_divergence(const Box& box,
                            FArrayBox& U,
                            FArrayBox& du,
                            EB_OPTIONAL(const EBCellFlagFab &flag,)
                            const Real *dx) const
{
    BL_PROFILE("State::calc_divergence");
#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();

    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<int,3> index = {0,0,0};
    Array4<Real> const& U4 = U.array();
    Array4<Real> const& du4 = du.array();
    Real dxinv;

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    continue;
                }
#endif
                // cycle over dimensions
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    index[d] = 1;
                    dxinv = 1.0/dx[d];

#ifdef AMREX_USE_EB
                    // left cell missing
                    if (f4(i,j,k).isDisconnected(-index[0],-index[1],-index[2])) {
                        du4(i,j,k) += dxinv*(U4(i+index[0],j+index[1],k+index[2],d) - U4(i,j,k,d));
                        index[d] = 0;
                        continue;
                    }

                    // right cell missing
                    if (f4(i,j,k).isDisconnected(index[0],index[1],index[2])) {
                        du4(i,j,k) += dxinv*(U4(i,j,k,d) - U4(i-index[0],j-index[1],k-index[2],d));
                        index[d] = 0;
                        continue;
                    }
#endif
                    du4(i,j,k) += 0.5*dxinv*(U4(i+index[0],j+index[1],k+index[2],d) - U4(i-index[0],j-index[1],k-index[2],d));
                    //                    du4(i,j,k) += dxinv*(U4(i,j,k,d) - U4(i-index[0],j-index[1],k-index[2],d));
                    index[d] = 0;
                }
            }
        }
    }

    return;
}


void State::write_info(nlohmann::json &js) const
{

    // write out stuff that is common to all states

    js["name"] = name;

    int tp = get_type();
    js["type_idx"] = tp;

    js["type"] = get_tag();

    js["global_idx"] = global_idx;

    if (shock_idx > -1) {
        js["shock_idx"] = shock_idx;
        shock_detector->write_info(js);
    }

    js["refine_grad_max_lev"] = refine_grad_max_lev;

    js["reflux"] = reflux;


    if (!refine_grad_threshold.empty()) {
        Vector<int> prim;
        Vector<int> idx;
        Vector<Real> thr;
        for (const auto &tup : refine_grad_threshold) {

            prim.push_back((int)std::get<0>(tup));
            idx.push_back(std::get<1>(tup));
            thr.push_back(std::get<2>(tup));
        }

        nlohmann::json& grp = js["refine_grad_threshold"];


        grp["is_primitive"] = prim;;
        grp["index"] = idx;
        grp["threshold"] = thr;
    }
}

std::string State::str() const
{
    std::stringstream msg;

    msg <<  "* "<< get_tag() <<" state '" << name << "' active:\n";

    msg << "    reconstruction : "<<reconstruction->get_tag()<<"\n";

    if (flux_solver) {
        msg << "    flux solver    : "<<flux_solver->get_tag()<<"\n";
    }

#ifdef AMREX_USE_EB
    msg << "    eb divergence  : "<<eb_div->str()<<"\n";
#endif

    msg << "    boundary conditions: \n";
    msg << "      names=" << vec2str(get_prim_names()) << "\n";
    msg << boundary_conditions.str("      ");

    msg << "    refinement: ";

    if (refine_grad_threshold.empty())
        msg << "nil";

    msg << "\n";

    for (const auto &v : refine_grad_threshold) {
        msg << "      ";
        if (std::get<bool>(v)) {
            msg << get_prim_names()[std::get<int>(v)];
        } else {
            msg << get_cons_names()[std::get<int>(v)];
        }
        msg << "(>" << num2str(std::get<3>(v)) <<") " << num2str(std::get<2>(v)) <<"\n";
    }



    msg << "    initial conditions: \n";
    std::pair<bool, int > check;
    for (const auto &f : functions) {
        msg << "      " << f.first << "=" << f.second;

        // check if this is also a dynamic function
        check = findInVector(get_prim_names(), f.first); // get the index
        if (check.first) {
            if ( dynamic_functions.find(check.second) != dynamic_functions.end() ) {
                msg << " (dynamic)";
            }
        }
        msg << "\n";
    }

    if (!associated_state.empty()) {
        msg << "    associated states:\n";
        for (const auto &a : associated_state) {
            msg << "      ";

            switch(a.first) {
                case AssociatedType::Ion:
                    msg << "ion";
                    break;
                case AssociatedType::Electron:
                    msg << "electron";
                    break;
                case AssociatedType::Neutral:
                    msg << "neutral";
                    break;
                case AssociatedType::Field:
                    msg << "field";
                    break;
                case AssociatedType::MHD:
                    msg << "MHD";
                    break;
                case AssociatedType::isNull:
                    msg << "null";
                    break;
            }

            Vector<std::string> names;

            for (const int& idx : a.second) {
                names.push_back(GD::get_state(idx).name);
            }

            msg << "=" << vec2str(names);
            msg << "\n";
        }
    }

    return msg.str();
}


PhysicsFactory<State>& GetStateFactory() {
    static PhysicsFactory<State> F;
    return F;
}
