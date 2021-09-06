#ifndef TRACER_H
#define TRACER_H
#ifdef AMREX_PARTICLES

#include "MFP_lagrangian.H"

using namespace amrex;

class TracerParticle : public LagrangianState
{
public:
    TracerParticle();
    TracerParticle(const sol::table& def);
    virtual void init_data(MFP* mfp, const Real time) override;
    virtual void init(MFP* mfp, bool make_particles=true) override;

    virtual Real get_allowed_time_step(MFP* mfp) const override {return std::numeric_limits<Real>::max();}
    void push_particles(const int level,
                        MFIter& mfi,
                        const FArrayBox& prim,
                        const Geometry geom,
                        const Real dt
                    #ifdef AMREX_USE_EB
                        ,const EBCellFlagFab& flag
                    #endif
                        ) const;

    virtual Vector<std::string> get_plot_output_names() const override;
    virtual void get_plot_output(MFP *mfp, MultiFab &plot_data, std::vector<std::string> &plot_names) const override;

    virtual void checkpoint(const std::string& dir) override;
    virtual void restart(const std::string& dir) override;
    virtual void redistribute(int level=0, int finest_level=-1, int ngrow=0, int local=0) override;
    virtual void clear() override;

    virtual StateType get_type() const override {return StateType::TracerParticle;}
    virtual void write_info(nlohmann::json& js) const override;

    enum class ParticleIdxR : int {
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

    static Vector<std::string> particle_real_names;
    static Vector<std::string> particle_int_names;

    Vector<RealArray> initial_positions;
    int verbosity;



    std::unique_ptr<AmrParticleContainer<+ParticleIdxR::NUM, +ParticleIdxI::NUM>> particles;
};

typedef AmrParticleContainer<+TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM> AmrTParContType;
typedef ParIter< +TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM, 0, 0 > TParIterType;
typedef ParticleContainer< +TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM, 0, 0 > TParContainerType;
typedef Particle< +TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM > TParticleType;
typedef ArrayOfStructs< +TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM > TParAoSType;
typedef ParticleTile< +TracerParticle::ParticleIdxR::NUM, +TracerParticle::ParticleIdxI::NUM, 0, 0 > TParTileType;

#endif
#endif // TRACER_H
