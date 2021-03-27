#include "MFP_reactions.H"
#include "MFP_global.H"

using GD = GlobalData;

std::string Ionisation::tag = "ionisation";
bool Ionisation::registered = GetSourceTermFactory().Register(Ionisation::tag, SourceTermBuilder<Ionisation>);


Ionisation::Ionisation(){}

Ionisation::Ionisation(const sol::table& def) {

  name = def.get<std::string>("name");

  if (!Ionisation::valid_solver(def["solver"])) {
      Abort("Error: Source '"+name+"' needs a different solver");
  }

  Vector<int> index;
  Vector<std::string> includes;

  get_includes(def, &Ionisation::valid_state, includes, index);

  offsets.resize(index.size());
  for (int idx=0; idx<index.size(); ++idx) {
      offsets[idx].local = idx;
      offsets[idx].global = index[idx];
  }

  if (index.size() < 3) {
      Abort("Error: For source `"+name+"`, you must define the states in order: electron, ion, neutral");
  }


  Real ionz_binding_nrg = def["ionz_binding_nrg"];
  if (!ionz_binding_nrg) Abort("Error: ionz_binding_nrg must be defined for the following ionisation source term: `"+name+"`.");

  std::string ionz_dens_fit_file =
      def["ionz_dens_fit_file"];

  if (ionz_dens_fit_file.length() == 0) {
    Abort(
        "Error: ionz_dens_fit_file must be defined for the "
        "following ionisation source term: `"+name+"`.");
  }

  dens_interp = Interp2D(ionz_dens_fit_file);

  sol::optional<std::string> ionz_mom_fit_file =  def["ionz_mom_fit_file"];

  if (ionz_mom_fit_file) {
      mom_interp = Interp2D(ionz_mom_fit_file.value());

      sol::optional<std::string> ionz_nrg_fit_file =  def["ionz_nrg_fit_file"];
      if (ionz_nrg_fit_file) {
          nrg_interp = Interp2D(ionz_nrg_fit_file.value());
      } else {
          Warning("Warning: ionz_nrg_fit_file for source `"+name+"` not defined. Ignoring the energy source terms");
      }

  } else {
      Warning("Warning: ionz_mom_fit_file for source `"+name+"` not defined. Ignoring the momentum and energy source terms");
  }

}

Ionisation::~Ionisation()
{
    // do nothing
}

Vector<Real> Ionisation::ionisation(const Vector<Real> &y0,
                                    const Vector<OffsetIndex> &offsets,
                                    const Interp2D &dens_interp,
                                    const Interp2D &mom_interp,
                                    const Interp2D &nrg_interp,
                                    const Real &bind_nrg) {
  Real Debye = GD::Debye;
  Real n0 = GD::n0;

  int num_hydro = offsets.size();

  // vector for hydro primitive values
  Vector<Vector<Real>> hydro_prim(offsets.size());

  int nc;
  for (const auto &idx : offsets) {
    State &istate = GD::get_state(idx.global);

    // get a copy of the conserved variables
    Vector<Real> &U = hydro_prim[idx.local];
    U.resize(istate.n_cons());
    for (int i = 0; i < U.size(); ++i) {
      U[i] = y0[idx.solver + i];
    }

    // convert to primitive
    istate.cons2prim(U);
  }

  Vector<Real> ydot(y0.size());

  // grab properties for electron, ion, neutral
  const OffsetIndex &offset_e = offsets[0];
  State &state_e = GD::get_state(offset_e.global);
  Vector<Real> &Q_e = hydro_prim[offset_e.local];

  Real rho_e = Q_e[+HydroState::PrimIdx::Density];
  Real alpha_e = Q_e[+HydroState::PrimIdx::Alpha];
  Real m_e = state_e.get_mass(alpha_e);
  Real p_e = Q_e[+HydroState::PrimIdx::Prs];
  Real n_e = rho_e / m_e;
  Real T_e = p_e / n_e;

  const OffsetIndex &offset_i = offsets[1];
  State &state_i = GD::get_state(offset_i.global);
  Vector<Real> &Q_i = hydro_prim[offset_i.local];

  Real rho_i = Q_i[+HydroState::PrimIdx::Density];
  Real alpha_i = Q_i[+HydroState::PrimIdx::Alpha];
  Real m_i = state_i.get_mass(alpha_i);
  Real p_i = Q_i[+HydroState::PrimIdx::Prs];
  Real n_i = rho_i / m_i;
  Real T_i = p_i / n_i;

  const OffsetIndex &offset_n = offsets[2];
  State &state_n = GD::get_state(offset_n.global);
  Vector<Real> &Q_n = hydro_prim[offset_n.local];

  Real rho_n = Q_n[+HydroState::PrimIdx::Density];
  Real alpha_n = Q_n[+HydroState::PrimIdx::Alpha];
  Real m_n = state_n.get_mass(alpha_n);
  Real p_n = Q_n[+HydroState::PrimIdx::Prs];
  Real n_n = rho_n / m_n;
  Real T_n = p_n / n_n;

  Real M = m_n;
  Real mu = m_e;
  Real muin = m_i * m_n / (m_i + m_n);

  // calculate x,y,z velocities and dot products
  Real du2 = 0.0;
  Real du_dot_U0 = 0.0;
  Real U02 = 0.0;
  Array<Real, 3> du;
  Array<Real, 3> U0;
  for (int i = 0; i < 3; ++i) {  // Xvel is 0, Yvel is 1, Zvel is 2
    du[i] = Q_e[+HydroState::PrimIdx::Xvel + i] -
            Q_n[+HydroState::PrimIdx::Xvel + i];


    du2 += du[i] * du[i];

    U0[i] = (m_e / M) * Q_e[+HydroState::PrimIdx::Xvel + i] +
            (m_n / M) * Q_n[+HydroState::PrimIdx::Xvel + i];

    du_dot_U0 += du[i] * U0[i];

    U02 += U0[i] * U0[i];
  }

  Real lam = 0.5 * mu * du2 / T_e;

  // TODO: move this stuff so MFP_config converts the dimensional tables to non-dimensional
  Real T_ref = GD::T_ref;
  Real x_ref = GD::x_ref;
  Real n_ref = GD::n_ref;
  Real u_ref = GD::u_ref;

  // Calculate the reaction rate coefficient based on the fit pars,
  Real T_e_real = T_e * T_ref;

  Real reaction_rate = n_e * n_n * dens_interp.interpolate(T_e_real, lam) * x_ref / u_ref * n_ref;

  if (mom_interp.is_valid) {
    // Calculate the momentum rate coefficient based on the fit pars,
    Real mom_rate = n_e * n_n * mom_interp.interpolate(T_e_real, lam) * x_ref / u_ref * n_ref;

    Real kappa = reaction_rate - mom_rate;
    Real TnTeTemu = (T_n - T_e) / T_e * mu;

    for (int i = 0; i < 3; ++i) {
      ydot[offset_e.solver + +HydroState::ConsIdx::Xmom + i] +=
          -mu * mom_rate * du[i];
      ydot[offset_i.solver + +HydroState::ConsIdx::Xmom + i] +=
          M * reaction_rate * U0[i] + TnTeTemu * kappa * du[i];
      ydot[offset_n.solver + +HydroState::ConsIdx::Xmom + i] +=
          -M * reaction_rate * U0[i] - TnTeTemu * kappa * du[i] +
          mu * mom_rate * du[i];
    }

    if (nrg_interp.is_valid) {
      // Calculate the energy rate coefficient based on the fit pars,
      Real en_rate = n_e * n_n * nrg_interp.interpolate(T_e_real, lam) * x_ref / u_ref * n_ref;


      Real J = en_rate - lam * mom_rate;
      Real W = en_rate - 2 * lam * mom_rate + lam * reaction_rate;

      // COM Temperature and Total Energy
      Real Tstar = M * T_e * T_n / (m_e * T_n + m_n * T_e);
      Real Estar = 0.5 * M * U02 + 1.5 * Tstar;

      ydot[offset_n.solver + +HydroState::ConsIdx::Eden] +=
          -reaction_rate * Estar - TnTeTemu * (T_n - T_e) * W / M +
          TnTeTemu * kappa * du_dot_U0 + mu * mom_rate * du_dot_U0 -
          2 * mu * (T_n - T_e) * J / M;
      ydot[offset_e.solver + +HydroState::ConsIdx::Eden] +=
          -reaction_rate * bind_nrg - mu * mom_rate * du_dot_U0 +
          2 * mu * (T_n - T_e) * J / M;
      ydot[offset_i.solver + +HydroState::ConsIdx::Eden] +=
          reaction_rate * Estar + TnTeTemu * (T_n - T_e) * W / M -
          TnTeTemu * kappa * du_dot_U0;
    }
  }

  ydot[offset_n.solver + +HydroState::ConsIdx::Density] -= m_n * reaction_rate;
  ydot[offset_i.solver + +HydroState::ConsIdx::Density] += m_i * reaction_rate;
  ydot[offset_e.solver + +HydroState::ConsIdx::Density] += m_e * reaction_rate;

  return ydot;
}

int Ionisation::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{

    const int n_terms = y0.size();

    // copy to eigen vector
    Vector<Real> yy(n_terms);
    for (int i=0; i<n_terms; ++i) {
        yy[i] = y0[i];
    }

    // call source function
    ydot = ionisation(yy, offsets, dens_interp, mom_interp, nrg_interp, bind_nrg);

    return 0;
}

int Ionisation::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{

    const int n_terms = y0.size();

    Vector<Real> ydot(n_terms);
    num_jac(x, y, z, t, y0, ydot, J);

    return 0;
}

bool Ionisation::valid_state(const int global_idx)
{

    State &istate = GD::get_state(global_idx);

    switch (istate.get_type()) {
    case +StateType::isHydro:
        return true;
    case +StateType::isHydro2P:
        return true;
    default:
        return false;
    }
}

bool Ionisation::valid_solver(const int solve_idx)
{
    return true;
}

std::string Recombination::tag = "recombination";
bool Recombination::registered = GetSourceTermFactory().Register(Recombination::tag, SourceTermBuilder<Recombination>);


Recombination::Recombination(){}

Recombination::Recombination(const sol::table &def)
{

    name = def.get<std::string>("name");

    if (!Recombination::valid_solver(def["solver"])) {
        Abort("Error: Source '"+name+"' needs a different solver");
    }

    Vector<int> index;
    Vector<std::string> includes;

    get_includes(def, &Ionisation::valid_state, includes, index);

    offsets.resize(index.size());
    for (int idx=0; idx<index.size(); ++idx) {
        offsets[idx].local = idx;
        offsets[idx].global = index[idx];
    }

    if (includes.size() < 3) {
        Abort("Error: For source `"+name+"`, you must define the states in order: electron, ion, neutral");
    }

    bind_nrg = def["recomb_binding_nrg"];
    if (!bind_nrg) Abort("Error: recomb_binding_nrg must be defined for the following recombination source term: `"+name+"`");

    std::string recomb_dens_fit_file = def["recomb_dens_fit_file"];

    if (recomb_dens_fit_file.length() == 0) {
        Abort("Error: recomb_dens_fit_file must be defined for the recombination source term: `"+name+"`");
    }

    std::string recomb_R_fit_file = def["recomb_R_fit_file"];

    if (recomb_R_fit_file.length() == 0) {
        Abort("Error: recomb_R_fit_file must be defined for the recombination source term: `"+name+"`");
    }
    std::string recomb_K_fit_file = def["recomb_K_fit_file"];

    if (recomb_K_fit_file.length() == 0) {
        Abort("Error: recomb_K_fit_file must be defined for the recombination source term: `"+name+"`");
    }

    std::string recomb_W_fit_file = def["recomb_W_fit_file"];

    if (recomb_W_fit_file.length() == 0) {
        Abort("Error: recomb_W_fit_file must be defined for the recombination source term: `"+name+"`");
    }

    std::string recomb_J_fit_file = def["recomb_J_fit_file"];

    if (recomb_J_fit_file.length() == 0) {
        Abort("Error: recomb_J_fit_file must be defined for the recombination source term: `"+name+"`");
    }

    dens_interp = Interp2D(recomb_dens_fit_file);
    R_interp = Interp2D(recomb_R_fit_file);
    K_interp = Interp2D(recomb_K_fit_file);
    W_interp = Interp2D(recomb_W_fit_file);
    J_interp = Interp2D(recomb_J_fit_file);

}

Recombination::~Recombination()
{
    // do nothing
}

Vector<Real> Recombination::recombination(const Vector<Real> &y0, const Vector<OffsetIndex> &offsets, const Interp2D &dens_interp, const Interp2D &R_interp, const Interp2D &K_interp, const Interp2D &W_interp, const Interp2D &J_interp, const Real &bind_nrg)
{

    Real Debye = GD::Debye;
    Real n0 = GD::n0;

    int num_hydro = offsets.size();

    // vector for hydro primitive values
    Vector<Vector<Real>> hydro_prim(offsets.size());


    int nc;
    for (const auto &idx : offsets) {

        State &istate = GD::get_state(idx.global);

        // get a copy of the conserved variables
        Vector<Real> &U = hydro_prim[idx.local];
        U.resize(istate.n_cons());
        for (int i=0; i<U.size(); ++i) {
            U[i] = y0[idx.solver+i];
        }

        // convert to primitive
        istate.cons2prim(U);
    }

    Vector<Real> ydot(y0.size());

    // grab properties for electron, ion, neutral
    const OffsetIndex &offset_e = offsets[0];
    State &state_e = GD::get_state(offset_e.global);
    Vector<Real> &Q_e = hydro_prim[offset_e.local];

    Real rho_e = Q_e[+HydroState::PrimIdx::Density];
    Real alpha_e = Q_e[+HydroState::PrimIdx::Alpha];
    Real m_e = state_e.get_mass(alpha_e);
    Real p_e = Q_e[+HydroState::PrimIdx::Prs];
    Real n_e = rho_e/m_e;
    Real T_e = p_e/n_e;

    const OffsetIndex &offset_i = offsets[1];
    State &state_i = GD::get_state(offset_i.global);
    Vector<Real> &Q_i = hydro_prim[offset_i.local];

    Real rho_i = Q_i[+HydroState::PrimIdx::Density];
    Real alpha_i = Q_i[+HydroState::PrimIdx::Alpha];
    Real m_i = state_i.get_mass(alpha_i);
    Real p_i = Q_i[+HydroState::PrimIdx::Prs];
    Real n_i = rho_i/m_i;
    Real T_i = p_i/n_i;

    const OffsetIndex &offset_n = offsets[2];
    State &state_n = GD::get_state(offset_n.global);
    Vector<Real> &Q_n = hydro_prim[offset_n.local];

    Real rho_n = Q_n[+HydroState::PrimIdx::Density];
    Real alpha_n = Q_n[+HydroState::PrimIdx::Alpha];
    Real m_n = state_n.get_mass(alpha_n);
    Real p_n = Q_n[+HydroState::PrimIdx::Prs];
    Real n_n = rho_n/m_n;
    Real T_n = p_n/n_n;

    Real M = m_n;
    Real mu = m_e;

    // calculate x,y,z velocities and dot products
    Real du2 = 0.0;
    Real du_dot_U1 = 0.0;
    Real U12 = 0.0;
    Array<Real,3> du;
    Array<Real,3> u_e;
    Array<Real,3> u_i;
    Array<Real,3> U1;
    for (int i=0; i<3; ++i) { // Xvel is 0, Yvel is 1, Zvel is 2
        u_e[i] = Q_e[+HydroState::PrimIdx::Xvel + i];
        u_e[i] = Q_i[+HydroState::PrimIdx::Xvel + i];
        du[i] = (-m_e / M)*Q_e[+HydroState::PrimIdx::Xvel + i] - (m_i/M)*Q_i[+HydroState::PrimIdx::Xvel + i];
        du2 += du[i]*du[i];

        U1[i] = 2 * (m_e/M) * Q_e[+HydroState::PrimIdx::Xvel + i] + (m_i/M) * Q_i[+HydroState::PrimIdx::Xvel + i];
        U12 += U1[i]*U1[i];

        du_dot_U1 += du[i] * U1[i];

    }

    Real lam = 0.5 * mu * du2 / T_e;
    Real T_ref = GD::T_ref;
    Real x_ref = GD::x_ref;
    Real n_ref = GD::n_ref;
    Real u_ref = GD::u_ref;

    // Calculate the reaction rate coefficient based on the fit pars,
    Real T_e_real = T_e * T_ref;

    Real reaction_rate = n_e * n_e * n_i * dens_interp.interpolate(T_e_real, lam) * n_ref * n_ref * x_ref / u_ref;
    Real K = K_interp.interpolate(T_e_real, lam) * n_ref * n_ref * x_ref / u_ref * n_e * n_e * n_i;
    Real R = R_interp.interpolate(T_e_real, lam) * n_ref * n_ref * x_ref / u_ref * n_e * n_e * n_i;
    Real W = W_interp.interpolate(T_e_real, lam) * n_ref * n_ref * x_ref / u_ref * n_e * n_e * n_i;
    Real J = J_interp.interpolate(T_e_real, lam) * n_ref * n_ref * x_ref / u_ref * n_e * n_e * n_i;

    ydot[offset_n.solver + +HydroState::ConsIdx::Density] +=  m_n * reaction_rate;
    ydot[offset_i.solver + +HydroState::ConsIdx::Density] -=  m_i * reaction_rate;
    ydot[offset_e.solver + +HydroState::ConsIdx::Density] -=  m_e * reaction_rate;

    for (int i = 0; i < 3; ++i) {
      ydot[offset_i.solver + +HydroState::ConsIdx::Xmom + i] += -M * reaction_rate * U1[i] - ((T_i - T_e) / T_e) * mu * K * du[i] + mu * R * du[i];
      ydot[offset_n.solver + +HydroState::ConsIdx::Xmom + i] += M * reaction_rate * U1[i] + ((T_i - T_e) / T_e) * mu* K *du[i];
      ydot[offset_e.solver + +HydroState::ConsIdx::Xmom + i] += -mu*R*du[i];
    }

    Real Tstar = T_i;
    Real Estar = 0.5 * M * U12 + 1.5 * Tstar;

    ydot[offset_n.solver + +HydroState::ConsIdx::Eden] += reaction_rate * Estar + mu/M * (T_i - T_e)*(T_i-T_e)/T_e*W+(T_i-T_e)/T_e*mu*K*du_dot_U1;
    ydot[offset_i.solver + +HydroState::ConsIdx::Eden] += -reaction_rate * Estar - mu/M * (T_i - T_e)*(T_i-T_e)/T_e*W-(T_i-T_e)/T_e*mu*K*du_dot_U1 + mu*R*du_dot_U1 - 2*mu/M * (T_i-T_e)*J;
    ydot[offset_e.solver + +HydroState::ConsIdx::Eden] += reaction_rate * bind_nrg - mu*R*du_dot_U1 + 2*mu/M * (T_i-T_e) * J;

    return ydot;
}

int Recombination::fun_rhs(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &ydot, Real dt) const
{

    const int n_terms = y0.size();

    // copy to eigen vector
    Vector<Real> yy(n_terms);
    for (int i=0; i<n_terms; ++i) {
        yy[i] = y0[i];
    }

    // call source function
    ydot = recombination(yy, offsets, dens_interp, R_interp, K_interp, W_interp, J_interp, bind_nrg);

    return 0;
}

int Recombination::fun_jac(Real x, Real y, Real z, Real t, Vector<Real> &y0, Vector<Real> &J) const
{

    const int n_terms = y0.size();

    Vector<Real> ydot(n_terms);
    num_jac(x, y, z, t, y0, ydot, J);

    return 0;
}

bool Recombination::valid_state(const int global_idx)
{

    State &istate = GD::get_state(global_idx);

    switch (istate.get_type()) {
    case +StateType::isHydro:
        return true;
    case +StateType::isHydro2P:
        return true;
    default:
        return false;
    }
}

bool Recombination::valid_solver(const int solve_idx)
{
    return true;
}
