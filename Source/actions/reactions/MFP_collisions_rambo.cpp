#include "MFP_collisions_rambo.H"
#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"

std::string CollisionsRambo::tag = "collisions_rambo";
bool CollisionsRambo::registered = GetActionFactory().Register(CollisionsRambo::tag, ActionBuilder<CollisionsRambo>);

CollisionsRambo::CollisionsRambo(){}
CollisionsRambo::~CollisionsRambo(){}

CollisionsRambo::CollisionsRambo(const int idx, const sol::table &def) : Collisions(idx, def)
{

    // grab any collision cross sections
    sol::table cross_sections = def["cross_sections"].get_or(sol::table());

    if (cross_sections.empty()) {

        // only need ccs with a neutral species
        for (const SpeciesInfo& info : species_info) {
            if (states[info.state_idx]->gas->charge[info.alpha_idx] != 0.0) {
                Abort("Error: collision cross section for neutral species '"+info.name+"' is required by action '"+name+"'");
            }
        }
    }

    const size_t n_species = species_info.size();
    Vector<std::string> species_names(n_species);
    ccs.resize(n_species);
    for (int a=0; a<n_species; ++a) {
        species_names[a] = species_info[a].name;
        ccs[a].resize(n_species,0);
    }

    // load the hydro cross sections
    std::pair<bool, int> index_a, index_b;
    Real sigma;
    for (const auto& a : cross_sections) {
        index_a = findInVector(species_names, a.first.as<std::string>());
        if (!index_a.first) {continue;}
        for (const auto& b : a.second.as<sol::table>()) {
            index_b = findInVector(species_names, b.first.as<std::string>());
            if (!index_b.first) {continue;}
            sigma = b.second.as<Real>()/(MFP::x_ref*MFP::x_ref); // non-dimensionalise
            ccs[index_a.second][index_b.second] = sigma;
            ccs[index_b.second][index_a.second] = sigma; // cross load as sigma_ab = sigma_ba
        }
    }

    Debye = MFP::Debye;
    n0 = MFP::n0;

    // coefficients
    c1 = 0.02116454531141366/(n0*pow(Debye,4)); // sqrt(2)/(12*pi^(3/2))*...
    c2 = 0.4135669939329334; // (2/(9*pi))**(1/3.0)
    c3 = 2.127692162140974*n0; // (4/3)*n0*sqrt(8/pi)
    lnC = 10.0; // Coulomb logarithm

    return;
}

void CollisionsRambo::calc_update()
{

    size_t n_species = species_info.size();

    for (int a = 0; a < n_species; ++a) {

        SpeciesInfo& sda = species_info[a];

        for (int b = a+1; b < n_species; ++b) {

            SpeciesInfo& sdb = species_info[b];

            const Real m_ab = (sda.m*sdb.m)/(sda.m + sdb.m);

            Real du2 = 0.0;
            Array<Real,3> du;
            for (int i=0; i<3; ++i) {
                du[i] = sdb.vel[i] - sda.vel[i];
                du2 += du[i]*du[i];
            }

            // collision frequency
            Real nu;
            if ((sda.q != 0) && (sdb.q != 0)) {
                const Real coeff_1 = c1*((sda.q2*sdb.q2*lnC)/(m_ab*sda.m));
                nu = sdb.n*coeff_1*pow(c2*du2 + sda.T/sda.m + sdb.T/sdb.m, -1.5);
            } else {
                const Real coeff_2 = (sdb.m*ccs[sda.idx][sdb.idx])/(sda.m + sdb.m);
                nu = sdb.n*coeff_2*c3*sqrt(sda.T/sda.m + sdb.T/sdb.m);
            }


            const Real S = 0.0;

            // energy exchange
            Real Q = m_ab*sda.n*nu*du2 + 3*sda.rho*nu/(sda.m + sdb.m)*(sdb.T - sda.T);

            // momentum exchange
            Array<Real,3> R;
            for (int n=0; n<3; ++n) {
                R[n] = sda.rho*nu*du[n];
                Q += R[n]*sda.vel[n];
            }

            // contribution to change of density
            sda.delta[+HydroDef::ConsIdx::Density] += S;
            sdb.delta[+HydroDef::ConsIdx::Density] -= S;

            // contribution to change of momentum
            for (int n=0; n<3; ++n) {
                sda.delta[+HydroDef::ConsIdx::Xmom+n] += R[n];
                sdb.delta[+HydroDef::ConsIdx::Xmom+n] -= R[n];
            }

            // contribution to change of energy
            sda.delta[+HydroDef::ConsIdx::Eden] += Q;
            sdb.delta[+HydroDef::ConsIdx::Eden] -= Q;
        }
    }
}
