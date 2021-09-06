#ifndef LAGRANGIANSTATE_H
#define LAGRANGIANSTATE_H
#ifdef AMREX_PARTICLES

#include <AMReX_AmrParticles.H>

#include "MFP.H"
#include "MFP_state.H"


class LagrangianState : public State
{
public:
    LagrangianState();
    ~LagrangianState();

    std::map<std::string,Optional3D1VFunction> functions;

    virtual StateClassification get_classification() const override {return StateClassification::Lagrangian;}

    virtual void init_from_lua() override;

    virtual bool is_transported() const {return true;}

    virtual StateType get_type() const = 0;

    virtual void get_plot_output(MFP* mfp, MultiFab &plot_data, std::vector<std::string> &plot_names) const {return;}

    virtual void init(MFP* mfp, bool make_particles=true) = 0;
    virtual void checkpoint(const std::string& dir) = 0;
    virtual void restart(const std::string& dir) = 0;
    virtual void redistribute(int level=0, int finest_level=-1, int ngrow=0, int local=0) = 0;
    virtual int num_source() {return 0;}
    virtual void calculate_source(MFIter& mfi, FArrayBox& S, Geometry &geom, int level) const {return;}
    virtual void clear() = 0;




    virtual void write_info(nlohmann::json& js) const override;

    static LagrangianState& get_state_global(const int idx) {return static_cast<LagrangianState&>(MFP::get_state(idx));}
    static LagrangianState& get_state(const int idx) {return static_cast<LagrangianState&>(MFP::get_state(MFP::eulerian_states[idx]));}

};
#endif
#endif // LAGRANGIANSTATE_H
