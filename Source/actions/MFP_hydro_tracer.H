#ifndef HYDROTRACER_H
#define HYDROTRACER_H
#ifdef AMREX_PARTICLES

    #include "MFP_action.H"
    #include "MFP_hydro.H"
    #include "MFP_tracer.H"

class HydroTracer : public Action
{
  public:
    HydroTracer();
    HydroTracer(const int idx, const sol::table& def);
    ~HydroTracer();

    void push_particles(TParTileType& ptile,
                        const FArrayBox& prim,
                        const Geometry geom,
                        const Real dt
    #ifdef AMREX_USE_EB
                        ,
                        const EBCellFlagFab& flag
    #endif
    );

    virtual void apply_spatial_derivative(MFP* mfp, const Real time, const Real dt) override;

    virtual std::string get_tag() const override { return tag; }
    virtual ActionType get_type() const override { return ActionType::HydroTracer; }

    HydroState* hydro_state;
    TracerParticle* tracer_state;

    static bool registered;
    static std::string tag;
};

#endif
#endif  // TRACER_H
