#include "MFP_field.H"

#include "MFP_global.H"
#include "MFP_Riemann_solvers.H"
#include "MFP_field_bc.H"

using GD = GlobalData;

Vector<std::string> FieldState::cons_names = {
    "x_D",
    "y_D",
    "z_D",
    "x_B",
    "y_B",
    "z_B",
    "phi",
    "psi",
    "mu",
    "ep",
};

Vector<std::string> FieldState::input_names = cons_names;

Vector<set_bc> FieldState::bc_set = {
    &set_x_D_bc,
    &set_y_D_bc,
    &set_z_D_bc,
    &set_x_B_bc,
    &set_y_B_bc,
    &set_z_B_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc,
    &set_scalar_bc,
};

Vector<int> FieldState::flux_vector_idx = {+FieldState::FluxIdx::Bx, +FieldState::FluxIdx::Dx};
Vector<int> FieldState::cons_vector_idx = {+FieldState::ConsIdx::Bx, +FieldState::ConsIdx::Dx};
Vector<int> FieldState::prim_vector_idx = {+FieldState::PrimIdx::Bx, +FieldState::PrimIdx::Dx};

std::string FieldState::tag = "field";
bool FieldState::registered = GetStateFactory().Register(FieldState::tag, StateBuilder<FieldState>);


FieldState::FieldState() : State(){}
FieldState::FieldState(const sol::table &def)
{
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
    sum_cons.resize(+ConsIdx::NUM);
}
FieldState::~FieldState(){}

void FieldState::init_from_lua()
{
    BL_PROFILE("FieldState::init_from_lua");
    sol::state& lua = GD::lua;

    const sol::table state_def = lua["states"][name];

    is_static = (bool) state_def["static"].get_or(0);

    State::init_from_lua();

    //
    // boundary conditions
    //

    const Vector<std::string> dir_name = {"x", "y", "z"};
    const Vector<std::string> side_name = {"lo", "hi"};

    // mapping between name and index for various groupings
    std::map<std::string,std::map<std::string, int>> bc2index;
    for (int i=+FieldState::PrimIdx::Dx; i<=+FieldState::PrimIdx::Dz; ++i) {
        bc2index["fill_D_bc"][input_names[i]] = i;
    }

    for (int i=+FieldState::PrimIdx::Bx; i<=+FieldState::PrimIdx::Bz; ++i) {
        bc2index["fill_B_bc"][input_names[i]] = i;
    }

    bc2index["fill_ep_bc"][input_names[+FieldState::PrimIdx::ep]] = +FieldState::PrimIdx::ep;
    bc2index["fill_mu_bc"][input_names[+FieldState::PrimIdx::mu]] = +FieldState::PrimIdx::mu;

    bc2index["fill_psi_bc"] = {{"psi",+FieldState::PrimIdx::psi}};
    bc2index["fill_phi_bc"] = {{"phi",+FieldState::PrimIdx::phi}};

    BoundaryState &bs = boundary_conditions;


    bs.fill_bc.resize(+FieldState::PrimIdx::NUM);
    bs.phys_fill_bc.resize(+FieldState::PrimIdx::NUM);

    for (int ax = 0; ax < AMREX_SPACEDIM; ++ax) {


        for (int lh=0; lh<2; ++lh) {

#ifdef AMREX_USE_EB
        bool is_symmetry = false;
#endif
            for (const auto &bc : bc2index) {

                // get the base boundary condition for cell centered values
                std::string side_bc = state_def["bc"][dir_name[ax]][side_name[lh]][bc.first].get_or<std::string>("outflow");
                int i_side_bc = bc_names.at(side_bc);
#ifdef AMREX_USE_EB
                if (i_side_bc == PhysBCType::symmetry || i_side_bc == 4) is_symmetry = true;
#endif
                // fill in the bc list for AMReX as well as gather any custom values/functions
                for (const auto &var : bc.second) {

                    if (lh==0) {
                        bs.phys_fill_bc[var.second].setLo(ax,i_side_bc);
                    } else {
                        bs.phys_fill_bc[var.second].setHi(ax,i_side_bc);
                    }

                    const sol::object bcv = state_def["bc"][dir_name[ax]][side_name[lh]][var.first].get_or(sol::object());

                    Optional3D1VFunction v;

                    // special case for phi and psi (set to zero in boundary)
                    if (var.second == +PrimIdx::phi || var.second == +PrimIdx::psi) {
                        get_udf(bcv,v,0.0);
                    } else {
                        v = get_udf(bcv);
                    }
                    bs.set(ax,input_names[var.second],lh,v);

                    // special case for inflow condition
                    if (i_side_bc == PhysBCType::inflow && !v.is_valid()) {
                        Abort("Setting '"+bc.first+" = inflow' requires all primitive variables to be defined, '" + var.first + "' is not defined");
                    }
                }
            }

#ifdef AMREX_USE_EB
            if (lh==0) {
                bs.eb_bc.setLo(ax,is_symmetry ? BCType::reflect_even : BCType::foextrap);
            } else {
                bs.eb_bc.setHi(ax,is_symmetry ? BCType::reflect_even : BCType::foextrap);
            }
#endif

        }
    }

    // check validity of inflow bc
    boundary_conditions.post_init();

    //
    // divergence handling
    //

    relative_div_speed = state_def["div_transport"].get_or(0.0);

    div_speed = relative_div_speed*GD::lightspeed;

    fastest_speed = std::max(GD::lightspeed, div_speed);

    // or we can use the projection method for divergence error control
    project_divergence = state_def["project_divergence"].get_or(0);


    return;
}

const Vector<std::string>& FieldState::get_cons_names() const
{
    return cons_names;
}

const Vector<std::string>& FieldState::get_prim_names() const
{
    return input_names;
}

const Vector<set_bc>& FieldState::get_bc_set() const
{
    return bc_set;
}

bool FieldState::cons2prim(Vector<Real>& U) const {return true;}
void FieldState::prim2cons(Vector<Real>& U) const {}

Real FieldState::get_energy_from_cons(const Vector<Real> &U) const
{
    return 0.5*(U[+ConsIdx::Bx]*U[+ConsIdx::Bx] + U[+ConsIdx::By]*U[+ConsIdx::By] + U[+ConsIdx::Bz]*U[+ConsIdx::Bz]
            + U[+ConsIdx::Dx]*U[+ConsIdx::Dx] + U[+ConsIdx::Dy]*U[+ConsIdx::Dy] + U[+ConsIdx::Dz]*U[+ConsIdx::Dz]);
}

RealArray FieldState::get_speed_from_cons(const Vector<Real>& U) const
{
    return get_speed_from_prim(U);
}

RealArray FieldState::get_speed_from_prim(const Vector<Real> &Q) const
{
    BL_PROFILE("FieldState::get_speed_from_prim");
    RealArray s;
    Real cf;
    if (is_static) {
        cf = 0.0;
    } else {
        Real mu = Q[+FieldState::PrimIdx::mu];
        Real ep = Q[+FieldState::PrimIdx::ep];
        cf = GD::lightspeed/std::sqrt(mu*ep); // local 'lightspeed'
    }

    s.fill(cf);

    return s;
}

void FieldState::get_state_values(const Box& box,
                                  const FArrayBox& src,
                                  std::map<std::string,FArrayBox>& out,
                                  Vector<std::string>& updated
                                  EB_OPTIONAL(,const FArrayBox& vfrac)
                                ) const
{
    BL_PROFILE("FieldState::get_state_values");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

#ifdef AMREX_USE_EB
    Array4<const Real> const& vf4 = vfrac.array();
#endif

    updated.resize(0);

    // check conserved variables
    std::map<std::string,int> cons_tags;
    for (int i=0; i<n_cons(); ++i) {
        const std::string s = cons_names[i];
        const std::string var_name = s+"-"+name;
        if ( out.find(var_name) == out.end() ) continue;
        cons_tags[var_name] = i;
        updated.push_back(var_name);
    }

    // additional variables
    Vector<std::string> other;

#ifdef AMREX_USE_EB
    const std::string vfrac_name = "vfrac-"+name;
    bool load_vfrac = out.find(vfrac_name) != out.end();
    if (load_vfrac) other.push_back(vfrac_name);
#endif

    updated.insert(updated.end(), other.begin(), other.end());

    std::map<std::string,Array4<Real>> out4;
    for (const std::string& s : updated) {
        out[s].resize(box, 1);
        out4[s] = out[s].array();
    }

    // temporary storage for retrieving the state data
    Vector<Real> S;

    Array4<const Real> const& src4 = src.array();

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (vf4(i,j,k) == 0.0) {
                    for (const std::string& s : updated) {
                        out4[s](i,j,k) = 0.0;
                    }
                    continue;
                }
#endif

                S = get_state_vector(src, i, j, k);

                if (!cons_tags.empty()) {
                    for (const auto& var : cons_tags) {
                        out4[var.first](i,j,k) = S[var.second];
                    }
                }

            #ifdef AMREX_USE_EB
                if (load_vfrac)  out4[vfrac_name](i,j,k)  = vf4(i,j,k);
            #endif
            }
        }
    }


    return;
}

void FieldState::calc_primitives(const Box& box,
                            FArrayBox& cons,
                            FArrayBox& prim,
                            const Real* dx,
                            const Real t,
                            const Real* prob_lo
                            EB_OPTIONAL(,const FArrayBox& vfrac)
                            ) const
{
    BL_PROFILE("FieldState::calc_primitives");
    Vector<Real> U;

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
    int nc = n_cons();
    int np = n_prim();

    for     (int k = lo.z; k <= hi.z; ++k) {
        z = prob_lo[2] + (k + 0.5)*dx[2];
        for   (int j = lo.y; j <= hi.y; ++j) {
            y = prob_lo[1] + (j + 0.5)*dx[1];
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                x = prob_lo[0] + (i + 0.5)*dx[0];

                U.resize(nc);

#ifdef AMREX_USE_EB
                if (vfrac4(i,j,k) <= 0.0) {

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
                        U[+FieldState::ConsIdx::psi] = 0.0;
                        U[+FieldState::ConsIdx::phi] = 0.0;
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
                cons2prim(U);

                // modify the primitives vector if needed and upload back to
                // the conserved vector
                if (!dynamic_functions.empty()) {
                    for (const auto &f : dynamic_functions) {
                        U[f.first] = (*f.second)(x, y, z, t);
                    }

                    // copy into primitive
                    for (int n=0; n<np; ++n) {
                        p4(i,j,k,n) = U[n];
                    }

                    // convert primitive to conserved
                    prim2cons(U);

                    // copy back into conserved array
                    for (int n=0; n<nc; ++n) {
                        s4(i,j,k,n) = U[n];
                    }
                }

                // copy into primitive
                for (int n=0; n<np; ++n) {
                    p4(i,j,k,n) = U[n];
                }
            }
        }
    }

    return;
}



int rollover(const int a, const int b)
{
    int c = a;
    for (int d = 0; d<b; ++d) {
        c += 1;
        if (c > 2) {
            c = 0;
        }
    }
    return c;
}

void FieldState::calc_reconstruction(const Box& box,
                                     FArrayBox &prim,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                                     Array<FArrayBox, AMREX_SPACEDIM> &rhi
                                     EB_OPTIONAL(,const EBCellFlagFab& flag)
                                     EB_OPTIONAL(,const FArrayBox &vfrac)
                                     ) const
{
    BL_PROFILE("FieldState::calc_reconstruction");
    // only operate over the components that need transporting
    const Array<int,8> transport_index = {+FieldState::PrimIdx::Dx,
                                          +FieldState::PrimIdx::Dy,
                                          +FieldState::PrimIdx::Dz,
                                          +FieldState::PrimIdx::Bx,
                                          +FieldState::PrimIdx::By,
                                          +FieldState::PrimIdx::Bz,
                                          +FieldState::PrimIdx::phi,
                                          +FieldState::PrimIdx::psi};

    const Array<int,2> stationary_index = {+FieldState::PrimIdx::ep,
                                           +FieldState::PrimIdx::mu,
                                          };

    int N = prim.nComp();

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    const Array4<Real> &p4 = prim.array();

    Vector<Real> stencil(reconstruction->stencil_length);
    int offset = reconstruction->stencil_length/2;
    Array<int,3> stencil_index;
    Real lo_face, hi_face;

#ifdef AMREX_USE_EB
    std::vector<std::array<int,3>> grab;
    multi_dim_index({-1,AMREX_D_PICK(0,-1,-1),AMREX_D_PICK(0,0,-1)},
    {1,AMREX_D_PICK(0, 1, 1),AMREX_D_PICK(0,0, 1)},
                    grab, false);

    Array4<const EBCellFlag> const& f4 = flag.array();
    // do we need to check our stencil for covered cells?
    bool check_eb = flag.getType() != FabType::regular;
#endif

    // get a copy of the D/B form
    FArrayBox db(prim, MakeType::make_deep_copy, 0, prim.nComp());
    const Array4<const Real> c4 = db.array();

    // get E/H from D/B
    for     (int k = lo.z-AMREX_D_PICK(0,0,offset); k <= hi.z+AMREX_D_PICK(0,0,offset); ++k) {
        for   (int j = lo.y-AMREX_D_PICK(0,offset,offset); j <= hi.y+AMREX_D_PICK(0,offset,offset); ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x-offset; i <= hi.x+offset; ++i) {

#ifdef AMREX_USE_EB
                if (f4(i,j,k).isCovered()) {
                    continue;
                }
#endif

                for (int ii = 0; ii<3; ++ii) {
                    p4(i,j,k,+FieldState::PrimIdx::Dx+ii) /= p4(i,j,k,+FieldState::PrimIdx::ep);
                    p4(i,j,k,+FieldState::PrimIdx::Bx+ii) /= p4(i,j,k,+FieldState::PrimIdx::mu);
                }
            }
        }
    }

    int cnt;
    //    Real ave;
    Array<int,3> of;
    Array4<const Real> s4;

    // cycle over dimensions
    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        of.fill(0);
        of[d] = 1;

        // make sure our arrays for putting lo and hi reconstructed values into
        // are the corect size
        rlo[d].resize(box, N);
        rhi[d].resize(box, N);



        Array4<Real> const& lo4 = rlo[d].array();
        Array4<Real> const& hi4 = rhi[d].array();

        // straight copy for those quantities not being transported
        for (const int &n : stationary_index) {
            rlo[d].copy(prim, box, n, box, n, 1);
            rhi[d].copy(prim, box, n, box, n, 1);
        }

#ifdef AMREX_USE_EB
        if (check_eb) {
            for (const int &n : transport_index) {
                rlo[d].copy(prim, box, n, box, n, 1);
                rhi[d].copy(prim, box, n, box, n, 1);
            }
        }
#endif

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

                    cnt = 0;
                    for (const int &n : transport_index) {

                        // if the component direction aligns with the dimension then use D/B
                        // else use E/H
                        if (cnt == d) {
                            s4 = c4;
                        } else {
                            s4 = p4;
                        }

                        stencil_index.fill(0);
                        for (int s=0; s<reconstruction->stencil_length; ++s) {
                            stencil_index[d] = s - offset;
                            stencil[s] = s4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], n);
                        }

                        // perform reconstruction
                        reconstruction->get_face_values(stencil, lo_face, hi_face);

                        lo4(i,j,k,n) = lo_face;
                        hi4(i,j,k,n) = hi_face;

                    }
                }
            }

            // increment and/or reset
            cnt++;
            if (cnt == 3)
                cnt = 0;

        }

        // for each dimension apply the correction to retrieve D/B from E/H
        // rollover essentially increments the provided index through the
        // looping sequence [0, 1, 2]

        rlo[d].mult(rlo[d], +FieldState::PrimIdx::ep, rollover(+FieldState::PrimIdx::Dy,d), 1);
        rlo[d].mult(rlo[d], +FieldState::PrimIdx::ep, rollover(+FieldState::PrimIdx::Dz,d), 1);
        rlo[d].mult(rlo[d], +FieldState::PrimIdx::mu,  rollover(+FieldState::PrimIdx::By,d), 1);
        rlo[d].mult(rlo[d], +FieldState::PrimIdx::mu,  rollover(+FieldState::PrimIdx::Bz,d), 1);

        rhi[d].mult(rhi[d], +FieldState::PrimIdx::ep, rollover(+FieldState::PrimIdx::Dy,d), 1);
        rhi[d].mult(rhi[d], +FieldState::PrimIdx::ep, rollover(+FieldState::PrimIdx::Dz,d), 1);
        rhi[d].mult(rhi[d], +FieldState::PrimIdx::mu,  rollover(+FieldState::PrimIdx::By,d), 1);
        rhi[d].mult(rhi[d], +FieldState::PrimIdx::mu,  rollover(+FieldState::PrimIdx::Bz,d), 1);

    }

    // get D/B from  E/H (return to original form)
    for (int i = 0; i<3; ++i) {
        prim.mult(prim, +FieldState::PrimIdx::ep, +FieldState::PrimIdx::Dx+i, 1);
        prim.mult(prim, +FieldState::PrimIdx::mu, +FieldState::PrimIdx::Bx+i, 1);
    }


    return;
}

// given all of the available face values load the ones expected by the flux calc into a vector
void FieldState::load_state_for_flux(Vector<Array4<const Real>> &face,
                                               int i, int j, int k, Vector<Real> &S) const
{
    BL_PROFILE("FieldState::load_state_for_flux");

    // first get the primitives of this state
    Array4<const Real> const &f4 = face[global_idx];

    S[+FluxIdx::Dx] = f4(i,j,k,+PrimIdx::Dx);
    S[+FluxIdx::Dy] = f4(i,j,k,+PrimIdx::Dy);
    S[+FluxIdx::Dz] = f4(i,j,k,+PrimIdx::Dz);
    S[+FluxIdx::Bx] = f4(i,j,k,+PrimIdx::Bx);
    S[+FluxIdx::By] = f4(i,j,k,+PrimIdx::By);
    S[+FluxIdx::Bz] = f4(i,j,k,+PrimIdx::Bz);
    S[+FluxIdx::phi] = f4(i,j,k,+PrimIdx::phi);
    S[+FluxIdx::psi] = f4(i,j,k,+PrimIdx::psi);
    S[+FluxIdx::mu] = f4(i,j,k,+PrimIdx::mu);
    S[+FluxIdx::ep] = f4(i,j,k,+PrimIdx::ep);

    return;
}

void FieldState::write_info(nlohmann::json &js) const
{

    State::write_info(js);

    js["div_transport"] = relative_div_speed;

}
