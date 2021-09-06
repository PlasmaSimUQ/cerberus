#ifndef CHARGEDPARTICLE_H
#define CHARGEDPARTICLE_H
#ifdef AMREX_PARTICLES

#include "MFP_lagrangian.H"
#include "MFP_utility.H"
#include "MFP_distribution.h"

using namespace amrex;

class ChargedParticle : public LagrangianState
{
public:
    ChargedParticle();
    ChargedParticle(const sol::table& def);
    virtual void init_data(MFP* mfp, const Real time) override;
    virtual void init(MFP* mfp, bool make_particles=true) override;
    virtual Real get_allowed_time_step(MFP* mfp) const override;
    virtual void checkpoint(const std::string& dir) override;
    virtual void restart(const std::string& dir) override;
    virtual void redistribute(int level=0, int finest_level=-1, int ngrow=0, int local=0) override;
    virtual int num_source() override {return 3;}
    virtual void calculate_source(MFIter& mfi, FArrayBox& S, Geometry &geom, int level) const override;
    virtual void clear() override;

    void push_particles(MFIter& mfi,
                        const FArrayBox& prim,
                        Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                        Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                        const int E_idx,
                        const int B_idx,
                        const Real dt,
                        const Geometry geom,
                        const int level
                    #ifdef AMREX_USE_EB
                        ,const EBCellFlagFab& flag
                    #endif
                        );

    virtual Vector<std::string> get_plot_output_names() const override;
    virtual void get_plot_output(MFP *mfp, MultiFab &plot_data, std::vector<std::string> &plot_names) const override;

    virtual StateType get_type() const override {return StateType::ChargedParticle;}

#ifdef AMREX_USE_EB
    virtual void set_eb_bc(const sol::table &bc_def) override;
#endif

    virtual void write_info(nlohmann::json& js) const override;

    enum class ParticleIdxR : int {
        Mass=0,
        Charge,
        VX,
        VY,
        VZ,
        NUM
    };

    enum class ParticleIdxI : int {
        NUM=0
    };

    static std::string tag;
    static bool registered;

    int state_idx;
    int verbosity;

    std::unique_ptr<ParticleDistribution> initial;

    static Vector<std::string> particle_real_names;
    static Vector<std::string> particle_int_names;

    std::unique_ptr<AmrParticleContainer<+ParticleIdxR::NUM, +ParticleIdxI::NUM>> particles;
};

typedef AmrParticleContainer<+ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM> AmrCParContType;
typedef ParIter< +ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM, 0, 0 > CParIterType;
typedef ParticleContainer< +ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM, 0, 0 > CParContainerType;
typedef Particle< +ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM > CParticleType;
typedef ArrayOfStructs< +ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM > CParAoSType;
typedef ParticleTile< +ChargedParticle::ParticleIdxR::NUM, +ChargedParticle::ParticleIdxI::NUM, 0, 0 > CParTileType;

#endif
#endif // CHARGEDPARTICLE_H
