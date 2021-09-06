#ifndef SPECIES_H
#define SPECIES_H

#include "MFP_pdf.h"

#include "MFP_fillbc.H"

#include "json.hpp"

// ============================================================================

enum class ParticleDistributionInit : int {
    Distribution,
    Defined,
    NUM
};

class ParticleDistribution
{
public:

    ParticleDistribution();
    virtual ~ParticleDistribution();

    virtual const std::string& get_tag() const = 0;
    virtual const ParticleDistributionInit get_type() const = 0;
    virtual void write_info(nlohmann::json& js) const;
};

template <typename D>
std::unique_ptr<ParticleDistribution> ParticleDistributionBuilder(const sol::table& def)
{
    if (def["type"] == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

ClassFactory<ParticleDistribution>& GetParticleDistributionFactory();

// ============================================================================



class DistributionSpecies : public ParticleDistribution{
public:
    DistributionSpecies();
    DistributionSpecies(const sol::table& def);

    virtual const std::string& get_tag() const override {return tag;}
    virtual const ParticleDistributionInit get_type() const override {return ParticleDistributionInit::Distribution;}
    virtual void write_info(nlohmann::json& js) const override;

    static std::string tag;
    static bool registered;

    Optional3D1VFunction number_density;
    Optional3D1VFunction particles_per_cell;
    std::unique_ptr<PDF> velocity;
    Real charge, mass;
};

// ============================================================================

class DefinedSpecies : public ParticleDistribution{
public:
    DefinedSpecies();
    DefinedSpecies(const sol::table& def);

    virtual const std::string& get_tag() const override {return tag;}
    virtual const ParticleDistributionInit get_type() const override {return ParticleDistributionInit::Defined;}
    int num_particles() const {return pos.size();}
    virtual void write_info(nlohmann::json& js) const override;

    static std::string tag;
    static bool registered;

    Vector<RealArray> pos;
    Vector<Array<Real,3>> vel;
    Real charge, mass;
};


#endif // SPECIES_H
