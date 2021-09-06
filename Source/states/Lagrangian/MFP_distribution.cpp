#include "MFP_distribution.h"
#include "MFP.H"

// ============================================================================

ParticleDistribution::ParticleDistribution(){}

ParticleDistribution::~ParticleDistribution(){}

ClassFactory<ParticleDistribution>& GetParticleDistributionFactory()
{
    static ClassFactory<ParticleDistribution> F;
    return F;
}

void ParticleDistribution::write_info(nlohmann::json& js) const
{
}

// ============================================================================

std::string DistributionSpecies::tag = "distribution";
bool DistributionSpecies::registered = GetParticleDistributionFactory().Register(DistributionSpecies::tag, ParticleDistributionBuilder<DistributionSpecies>);

DistributionSpecies::DistributionSpecies(){}

DistributionSpecies::DistributionSpecies(const sol::table& def)
{

    charge = def["charge"];
    mass = def["mass"];

    number_density = get_udf(def["nd"]);
    particles_per_cell = get_udf(def["particles_per_cell"]);

    ClassFactory<PDF> pdf_fact = GetPDFFactory();

    velocity = pdf_fact.Build(def["distribution"],def);

    if (!velocity)
        Abort("Failed to read velocity distribution, must be one of "+vec2str(pdf_fact.getKeys()));

}

void DistributionSpecies::write_info(nlohmann::json& js) const
{
    ParticleDistribution::write_info(js);

    std::stringstream ss;

    ss << number_density;
    js["number_density"] = ss.str();
    ss.str(""); ss.clear();

    ss << particles_per_cell;
    js["particles_per_cell"] = ss.str();
    ss.str(""); ss.clear();

    js["velocity"] = velocity->get_tag();

    js["mass"] = mass;
    js["charge"] = charge;

}

// ============================================================================

std::string DefinedSpecies::tag = "defined";
bool DefinedSpecies::registered = GetParticleDistributionFactory().Register(DefinedSpecies::tag, ParticleDistributionBuilder<DefinedSpecies>);

DefinedSpecies::DefinedSpecies(){}

DefinedSpecies::DefinedSpecies(const sol::table& def)
{

    charge = def["charge"];
    mass = def["mass"];

    // lua table of particle position and velocity
    // {{x,y,z,u,v,w}, ...}
    sol::table particles = def["particles"];

    for (auto& p : particles) {
        sol::table dat = p.second;
        pos.push_back(RealArray({AMREX_D_DECL(dat[1],dat[2],dat[3])}));
        vel.push_back({AMREX_D_DECL(dat[4],dat[5],dat[6])});
    }

}

void DefinedSpecies::write_info(nlohmann::json& js) const
{
    ParticleDistribution::write_info(js);

    js["position"] = pos;
    js["velocity"] = vel;

    js["mass"] = mass;
    js["charge"] = charge;
}
