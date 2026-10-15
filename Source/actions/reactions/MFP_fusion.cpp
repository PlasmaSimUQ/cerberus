#include "MFP_fusion.H"

#include "MFP.H"
#include "MFP_utility.H"

std::string Fusion::tag = "fusion";
bool Fusion::registered = GetActionFactory().Register(Fusion::tag, ActionBuilder<Fusion>);

Fusion::Fusion() {}
Fusion::~Fusion() {}

Fusion::Fusion(const int idx, const sol::table& def) : Collisions(idx, def)
{
    BL_PROFILE("Fusion::Fusion");

    Vector<std::string> names;
    for (const SpeciesInfo& s : species_info) names.push_back(s.name);

    auto lookup = [&](const std::string& sname, const char* key) {
        const std::pair<bool, int> found = findInVector(names, sname);
        if (!found.first)
            Abort("Action '" + name + "': '" + key + "' entry '" + sname + "' is not in 'species'");
        return found.second;
    };

    const sol::table fuel = def["fuel"];
    if (fuel.size() != 2)
        Abort("Action '" + name + "': 'fuel' must name exactly two species (D-T only)");

    amrex::Print() << "\n#=============================================================#\n"
                   << "#  WARNING [" << name << "]: fusion action is D + T -> He4 + n ONLY.  #\n"
                   << "#  No other reactions; fuel/product accounting assumes this.  #\n"
                   << "#=============================================================#\n\n";

    d_idx = lookup(fuel[1], "fuel");
    t_idx = lookup(fuel[2], "fuel");
    p_idx = lookup(def["product"], "product");
    e_idx = electron_species();

    if ((species_info[d_idx].state_idx != species_info[t_idx].state_idx) ||
        (species_info[d_idx].state_idx != species_info[p_idx].state_idx))
        Abort("Action '" + name + "': fuel and product species must belong to the same state");

    alpha_fraction = def.get_or("alpha_fraction", 1.0);
    ion_fraction = def.get_or("ion_fraction", 0.0);
    max_burn_fraction = def.get_or("max_burn_fraction", 0.5);
    if ((alpha_fraction < 0.0) || (alpha_fraction > 1.0) || (ion_fraction < 0.0) || (ion_fraction > 1.0))
        Abort("Action '" + name + "': alpha_fraction and ion_fraction must lie in [0, 1]");

    const std::string partition = def.get_or<std::string>("alpha_partition", "fixed");
    fraley = (partition == "fraley");
    if (!fraley && (partition != "fixed"))
        Abort("Action '" + name + "': alpha_partition must be 'fixed' or 'fraley'");
    fraley_constant = def.get_or("fraley_constant", 32.0);

    const Real kB = 1.38064852e-23, e = 1.602176634e-19;
    kT_keV = kB * MFP::T_ref / e / 1.0e3;
    C_sv = 3.68e-18 * MFP::n_ref * MFP::t_ref;              // R_code = n_D n_T C_sv f(T_keV)
    E_alpha = 3.5e6 * e / (MFP::m_ref * MFP::u_ref * MFP::u_ref);

    if (species_info[d_idx].m + species_info[t_idx].m <= species_info[p_idx].m)
        Abort("Action '" + name + "': product heavier than the fuel, check the species masses");
}

void Fusion::calc_update()
{
    SpeciesInfo& D = species_info[d_idx];
    SpeciesInfo& T = species_info[t_idx];
    SpeciesInfo& P = species_info[p_idx];

    const Real T_keV = D.T * kT_keV;
    if ((T_keV <= 0.0) || (D.n <= 0.0) || (T.n <= 0.0)) return;

    const Real x = std::pow(T_keV, -1.0 / 3.0);
    Real R = D.n * T.n * C_sv * x * x * std::exp(-19.94 * x);  // reactions per volume per time

    // cold cells: the Gamow factor underflows to a denormal or zero; a burn
    // time above 1e30 code units is no burn, and n_min / R would overflow
    const Real n_min = std::min(D.n, T.n);
    if (!(R > 1.0e-30 * n_min)) return;

    if (dt_stage > 0.0) R = std::min(R, max_burn_fraction * n_min / dt_stage);

    // mass transfer; products inherit the state velocity v (shared by all sub-species)
    const Real v2 = D.vel[0] * D.vel[0] + D.vel[1] * D.vel[1] + D.vel[2] * D.vel[2];
    auto move = [&](SpeciesInfo& s, const Real dm) {
        s.delta[+HydroDef::ConsIdx::Density] += dm;
        s.delta[+HydroDef::ConsIdx::Xmom] += dm * D.vel[0];
        s.delta[+HydroDef::ConsIdx::Ymom] += dm * D.vel[1];
        s.delta[+HydroDef::ConsIdx::Zmom] += dm * D.vel[2];
        s.delta[+HydroDef::ConsIdx::Eden] += 0.5 * v2 * dm;
    };
    move(D, -D.m * R);
    move(T, -T.m * R);
    move(P, P.m * R);  // net state loss (m_D + m_T - m_He4) R leaves with the neutron

    const Real E_dep = alpha_fraction * E_alpha * R;
    Real f_i = ion_fraction;
    if (fraley) {
        const Real Te_keV = ((e_idx >= 0) ? species_info[e_idx].T : D.T) * kT_keV;
        f_i = 1.0 / (1.0 + fraley_constant / Te_keV);
    }
    P.delta[+HydroDef::ConsIdx::Eden] += f_i * E_dep;
    species_info[(e_idx >= 0) ? e_idx : p_idx].delta[+HydroDef::ConsIdx::Eden] += (1.0 - f_i) * E_dep;

    tau_min = std::min(tau_min, n_min / R);
    if (E_dep > 0.0) tau_min = std::min(tau_min, 1.5 * (D.n + T.n + P.n) * D.T / E_dep);
}
