#include "MFP_mhd.H"

#include "Eigen"
#include "MFP.H"
#include "MFP_diagnostics.H"
#include "MFP_lua.H"
#include "MFP_mhd_refine.H"
#include "MFP_transforms.H"

Vector<std::string> MHDState::cons_names = {
  "rho", "x_mom", "y_mom", "z_mom", "nrg", "x_B", "y_B", "z_B", "psi"};

Vector<std::string> MHDState::prim_names = {
  "rho", "x_vel", "y_vel", "z_vel", "p", "x_B", "y_B", "z_B", "psi"};

Array<int, 2> MHDState::cons_vector_idx = {+MHDDef::ConsIdx::Xmom, +MHDDef::ConsIdx::Bx};
Array<int, 2> MHDState::prim_vector_idx = {+MHDDef::PrimIdx::Xvel, +MHDDef::PrimIdx::Bx};

std::map<std::string, int> MHDState::bc_names = {{"interior", PhysBCType::interior},
                                                 {"inflow", PhysBCType::inflow},
                                                 {"outflow", PhysBCType::outflow},
                                                 {"symmetry", PhysBCType::symmetry},
                                                 {"slipwall", PhysBCType::slipwall},
                                                 {"noslipwall", PhysBCType::noslipwall}};

Vector<set_bc> MHDState::bc_set = {
  &set_scalar_bc,
  &set_x_vel_bc,
  &set_y_vel_bc,
  &set_z_vel_bc,
  &set_scalar_bc,
  &set_x_B_bc,
  &set_y_B_bc,
  &set_z_B_bc,
  &set_scalar_bc,
};

std::string MHDState::tag = "mhd";
bool MHDState::registered = GetStateFactory().Register(MHDState::tag, StateBuilder<MHDState>);

MHDState::MHDState() {}

MHDState::MHDState(const sol::table& def)
{
    BL_PROFILE("MHDState::MHDState");

    num_grow = 0;
    name = def.get<std::string>("name");
    global_idx = def.get<int>("global_idx");
}

MHDState::~MHDState() {}

void MHDState::set_udf()
{
    BL_PROFILE("MHDState::set_udf");

    using namespace std::placeholders;

    sol::state& lua = MFP::lua;

    sol::table state_def = lua["states"][name];

    // check if we have 'value' defined
    const sol::table value = state_def["value"].get_or(sol::table());

    if (!value.valid())
        Abort("State " + name + " does not have 'value' defined for initial conditions");

    functions.resize(+MHDDef::PrimIdx::NUM);

    // now for the primitives
    for (int i = 0; i < +MHDDef::PrimIdx::NUM; ++i) {
        std::string comp = prim_names[i];

        Optional3D1VFunction v;

        get_udf(value[comp], v, 0.0);

        functions[i] = v;
    }

    return;
}

void MHDState::set_flux()
{
    BL_PROFILE("MHDState::set_flux");

    if (!is_transported()) return;

    ClassFactory<MHDRiemannSolver> rfact = GetMHDRiemannSolverFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;
    state_def["gamma"] = gamma;
    state_def["div_transport"] = div_transport;

    std::string flux = state_def["flux"].get_or<std::string>("null");

    if (flux == "null")
        Abort("Flux option required for state '" + name + "'. Options are " +
              vec2str(rfact.getKeys()));

    flux_solver = rfact.Build(flux, state_def);

    if (!flux_solver)
        Abort("Invalid flux solver option '" + flux + "'. Options are " + vec2str(rfact.getKeys()));

    return;
}

void MHDState::set_shock_detector()
{
    BL_PROFILE("MHDState::set_shock_detector");

    ClassFactory<MHDShockDetector> sdfact = GetMHDShockDetectorFactory();

    sol::table sd_def = MFP::lua["states"][name]["shock_detector"].get_or(sol::table());

    if (!sd_def.valid()) return;

    sd_def["global_idx"] = global_idx;

    std::string sd_name = sd_def["name"].get_or<std::string>("");

    shock_detector = sdfact.Build(sd_name, sd_def);

    if (!sd_name.empty() && !shock_detector)
        Abort("Invalid shock_detector option '" + sd_name + "'. Options are " +
              vec2str(sdfact.getKeys()));
}

void MHDState::set_refinement()
{
    BL_PROFILE("MHDState::set_refinement");

    ClassFactory<Refinement> rfact = GetMHDRefinementFactory();

    sol::table r_def = MFP::lua["states"][name]["refinement"].get_or(sol::table());

    if (!r_def.valid()) return;

    r_def["global_idx"] = global_idx;

    std::string r_name = r_def["name"].get_or<std::string>("");

    refine = rfact.Build(r_name, r_def);

    if (!r_name.empty() && !refine)
        Abort("Invalid refinement option '" + r_name + "'. Options are " +
              vec2str(rfact.getKeys()));
}

void MHDState::init_from_lua()
{
    BL_PROFILE("MHDState::init_from_lua");

    EulerianState::init_from_lua();

    sol::state& lua = MFP::lua;

    const sol::table state_def = lua["states"][name];

    //
    // get gamma
    //

    gamma = state_def["gamma"];

    //
    // user defined functions
    //
    set_udf();

    //
    // domain boundary conditions
    //

    const Vector<std::string> dir_name = {"x", "y", "z"};
    const Vector<std::string> side_name = {"lo", "hi"};
    const Vector<std::string>& mhd_var = prim_names;
    const int N = mhd_var.size();

    BoundaryState& bs = boundary_conditions;
    bs.phys_fill_bc.resize(+MHDDef::PrimIdx::NUM);

    for (int ax = 0; ax < AMREX_SPACEDIM; ++ax) {
        for (int lh = 0; lh < 2; ++lh) {
            std::string side_bc =
              state_def["bc"][dir_name[ax]][side_name[lh]]["fill_mhd_bc"].get_or<std::string>(
                "outflow");
            int i_side_bc = bc_names.at(side_bc);

            // get any custom values/functions
            for (int j = 0; j < N; ++j) {
                if (lh == 0) {
                    bs.phys_fill_bc[j].setLo(ax, i_side_bc);
                } else {
                    bs.phys_fill_bc[j].setHi(ax, i_side_bc);
                }

                const sol::object v =
                  state_def["bc"][dir_name[ax]][side_name[lh]][mhd_var[j]].get_or(sol::object());
                Optional3D1VFunction f = get_udf(v);
                bs.set(ax, mhd_var[j], lh, f);

                // special case for inflow condition
                if (i_side_bc == PhysBCType::inflow && !f.is_valid()) {
                    Abort("Setting 'fill_mhd_bc = inflow' requires all primitive variables to be "
                          "defined, '" +
                          mhd_var[j] + "' is not defined");
                }
            }

#ifdef AMREX_USE_EB
            bool is_symmetry =
              (i_side_bc == PhysBCType::symmetry) || (i_side_bc == PhysBCType::slipwall);
            if (lh == 0) {
                bs.eb_bc.setLo(ax, is_symmetry ? BCType::reflect_even : BCType::foextrap);
            } else {
                bs.eb_bc.setHi(ax, is_symmetry ? BCType::reflect_even : BCType::foextrap);
            }
#endif
        }
    }

    // check validity of inflow bc
    boundary_conditions.post_init();

    //
    // divergence handling
    //

    div_transport = state_def["div_transport"].get_or(0.0);

    //
    // riemann solver
    //
    set_flux();

    //
    // shock detector
    //
    set_shock_detector();

    //
    // refinement
    //

    set_refinement();
}

void MHDState::variable_setup(Vector<int> periodic)
{
    BL_PROFILE("MHDState::variable_setup");

    boundary_conditions.fill_bc.resize(+MHDDef::PrimIdx::NUM);

    for (int icomp = 0; icomp < +MHDDef::ConsIdx::NUM; ++icomp) {
        set_bc s = bc_set[icomp];  // the function that sets the bc

        // make sure our periodicity isn't being overwritten
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (periodic[d]) {
                boundary_conditions.phys_fill_bc[icomp].setLo(d, PhysBCType::interior);
                boundary_conditions.phys_fill_bc[icomp].setHi(d, PhysBCType::interior);
            }
        }

        // grab the per component BCRec and apply the set bc function to it
        (*s)(boundary_conditions.fill_bc[icomp], boundary_conditions.phys_fill_bc[icomp]);
    }

    Vector<std::string> comp_names(+MHDDef::ConsIdx::NUM);
    for (int icomp = 0; icomp < +MHDDef::ConsIdx::NUM; ++icomp) {
        comp_names[icomp] = cons_names[icomp] + "-" + name;
    }

    int ng = num_grow;

#ifdef AMREX_USE_EB
    Interpolater* interp = &eb_cell_cons_interp;
#else
    Interpolater* interp = &cell_cons_interp;
#endif

    bool state_data_extrap = false;
    bool store_in_checkpoint = true;

    DescriptorList& desc_lst = MFP::get_desc_lst();

    data_idx = desc_lst.size();

    desc_lst.addDescriptor(data_idx,
                           IndexType::TheCellType(),
                           StateDescriptor::Point,
                           ng,
                           +MHDDef::ConsIdx::NUM,
                           interp,
                           state_data_extrap,
                           store_in_checkpoint);

    desc_lst.setComponent(data_idx, 0, comp_names, boundary_conditions.fill_bc, FillBC());

    if (MFP::verbosity >= 1) { Print() << str(); }
}

// in place conversion from conserved to primitive
bool MHDState::cons2prim(Vector<Real>& U, Vector<Real>& Q) const
{
    BL_PROFILE("MHDState::cons2prim");

    Real rho = U[+MHDDef::ConsIdx::Density];
    Real rhoinv = 1 / rho;
    Real u = U[+MHDDef::ConsIdx::Xmom] * rhoinv;
    Real v = U[+MHDDef::ConsIdx::Ymom] * rhoinv;
    Real w = U[+MHDDef::ConsIdx::Zmom] * rhoinv;

    Real Bx = U[+MHDDef::ConsIdx::Bx];
    Real By = U[+MHDDef::ConsIdx::By];
    Real Bz = U[+MHDDef::ConsIdx::Bz];

    Real nrg = U[+MHDDef::ConsIdx::Eden] -
               0.5 * (rho * (u * u + v * v + w * w) + Bx * Bx + By * By + Bz * Bz);

    Real p = nrg * (gamma - 1.0);

    Q[+MHDDef::PrimIdx::Density] = rho;
    Q[+MHDDef::PrimIdx::Xvel] = u;
    Q[+MHDDef::PrimIdx::Yvel] = v;
    Q[+MHDDef::PrimIdx::Zvel] = w;
    Q[+MHDDef::PrimIdx::Prs] = p;
    Q[+MHDDef::PrimIdx::Bx] = Bx;
    Q[+MHDDef::PrimIdx::By] = By;
    Q[+MHDDef::PrimIdx::Bz] = Bz;
    Q[+MHDDef::PrimIdx::psi] = U[+MHDDef::ConsIdx::psi];

    return prim_valid(Q);
}

void MHDState::prim2cons(Vector<Real>& Q, Vector<Real>& U) const
{
    BL_PROFILE("MHDState::prim2cons");

    Real rho = Q[+MHDDef::PrimIdx::Density];
    Real mx = rho * Q[+MHDDef::PrimIdx::Xvel];
    Real my = rho * Q[+MHDDef::PrimIdx::Yvel];
    Real mz = rho * Q[+MHDDef::PrimIdx::Zvel];
    Real p = Q[+MHDDef::PrimIdx::Prs];

    Real Bx = Q[+MHDDef::PrimIdx::Bx];
    Real By = Q[+MHDDef::PrimIdx::By];
    Real Bz = Q[+MHDDef::PrimIdx::Bz];

    U[+MHDDef::ConsIdx::Density] = rho;
    U[+MHDDef::ConsIdx::Xmom] = mx;
    U[+MHDDef::ConsIdx::Ymom] = my;
    U[+MHDDef::ConsIdx::Zmom] = mz;
    U[+MHDDef::ConsIdx::Eden] = p / (gamma - 1) + 0.5 * (mx * mx + my * my + mz * mz) / rho +
                                0.5 * (Bx * Bx + By * By + Bz * Bz);
    U[+MHDDef::ConsIdx::Bx] = Bx;
    U[+MHDDef::ConsIdx::By] = By;
    U[+MHDDef::ConsIdx::Bz] = Bz;
    U[+MHDDef::ConsIdx::psi] = Q[+MHDDef::PrimIdx::psi];
}

bool MHDState::prim_valid(const Vector<Real>& Q) const
{
    BL_PROFILE("MHDState::prim_valid");

    if ((Q[+MHDDef::PrimIdx::Density] <= 0.0) || (Q[+MHDDef::PrimIdx::Prs] <= 0.0)) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

bool MHDState::cons_valid(const Vector<Real>& U) const
{
    BL_PROFILE("MHDState::cons_valid");

    if ((U[+MHDDef::ConsIdx::Density] <= 0.0) || (U[+MHDDef::ConsIdx::Eden] <= 0.0)) {
        //        amrex::Abort("Primitive values outside of physical bounds!!");
        return false;
    }
    return true;
}

RealArray MHDState::get_speed_from_cons(const Vector<Real>& U) const
{
    BL_PROFILE("MHDState::get_speed_from_cons");

    Real rho = U[+MHDDef::ConsIdx::Density];
    Real ed = U[+MHDDef::ConsIdx::Eden];
    Real rhoinv = 1 / rho;

    Real ux = U[+MHDDef::ConsIdx::Xmom] * rhoinv;
    Real uy = U[+MHDDef::ConsIdx::Ymom] * rhoinv;
    Real uz = U[+MHDDef::ConsIdx::Zmom] * rhoinv;

    Real kineng = 0.5 * rho * (ux * ux + uy * uy + uz * uz);
    Real p = (ed - kineng) * (gamma - 1);

    Real Bx = U[+MHDDef::ConsIdx::Bx];
    Real By = U[+MHDDef::ConsIdx::By];
    Real Bz = U[+MHDDef::ConsIdx::Bz];

    Real B = std::sqrt(Bx * Bx + By * By + Bz * Bz);

    Real cf = std::sqrt((gamma * p + B) / rho);

    RealArray s = {AMREX_D_DECL(cf + std::abs(ux), cf + std::abs(uy), cf + std::abs(uz))};

    return s;
}

RealArray MHDState::get_speed_from_prim(const Vector<Real>& Q) const
{
    BL_PROFILE("MHDState::get_speed_from_prim");

    Real rho = Q[+MHDDef::PrimIdx::Density];
    Real ux = Q[+MHDDef::PrimIdx::Xvel];
    Real uy = Q[+MHDDef::PrimIdx::Yvel];
    Real uz = Q[+MHDDef::PrimIdx::Zvel];
    Real p = Q[+MHDDef::PrimIdx::Prs];

    Real Bx = Q[+MHDDef::PrimIdx::Bx];
    Real By = Q[+MHDDef::PrimIdx::By];
    Real Bz = Q[+MHDDef::PrimIdx::Bz];

    Real B = std::sqrt(Bx * Bx + By * By + Bz * Bz);

    Real cf = std::sqrt((gamma * p + B) / rho);

    RealArray s = {AMREX_D_DECL(cf + std::abs(ux), cf + std::abs(uy), cf + std::abs(uz))};

    return s;
}

void MHDState::calc_velocity(const Box& box,
                             FArrayBox& cons,
                             FArrayBox& prim
#ifdef AMREX_USE_EB
                             ,
                             const FArrayBox& vfrac
#endif
) const
{
    BL_PROFILE("MHDState::calc_velocity");

    prim.resize(box, AMREX_SPACEDIM);

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);
    Array4<Real> const& s4 = cons.array();
    Array4<Real> const& p4 = prim.array();

#ifdef AMREX_USE_EB
    Array4<const Real> const& vfrac4 = vfrac.array();
#endif

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
#ifdef AMREX_USE_EB
                if (vfrac4(i, j, k) == 0.0) continue;
#endif

                const Real irho = 1 / s4(i, j, k, +MHDDef::ConsIdx::Density);

                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    p4(i, j, k, n) = irho * s4(i, j, k, +MHDDef::ConsIdx::Xmom + n);
                }
            }
        }
    }

    return;
}

void MHDState::write_info(nlohmann::json& js) const
{
    BL_PROFILE("MHDState::write_info");

    EulerianState::write_info(js);

    // write out stuff that is common to all states

    js["type"] = tag;
    js["type_idx"] = +get_type();

    js["gamma"] = gamma;
}
