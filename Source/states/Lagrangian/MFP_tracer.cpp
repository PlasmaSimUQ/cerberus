#ifdef AMREX_PARTICLES
#include "MFP_tracer.H"
#include "sol.hpp"

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

#include "MFP_diagnostics.H"

Vector<std::string> TracerParticle::particle_real_names = {"x_vel", "y_vel", "z_vel"};
Vector<std::string> TracerParticle::particle_int_names = {};

std::string TracerParticle::tag = "tracer";
bool TracerParticle::registered = GetStateFactory().Register(TracerParticle::tag, StateBuilder<TracerParticle>);

TracerParticle::TracerParticle()
{

}

TracerParticle::TracerParticle(const sol::table& def)
{
    name = def["name"];
    global_idx = def["global_idx"];
    verbosity = def.get_or("verbosity",0);

    // grab the initial positions
    sol::table particle_def = def["particles"];

    for (auto& pd : particle_def) {
        sol::table dat = pd.second;
        initial_positions.push_back({AMREX_D_DECL(dat[1],dat[2],dat[3])});
    }
}

void TracerParticle::init_data(MFP* mfp, const Real time)
{
    if (mfp->get_level() == 0)
        init(mfp);
}

// this should only be called at level 0
void TracerParticle::init(MFP* mfp, bool make_particles)
{
    AmrCore* amr_core = mfp->get_parent();

    // generate the particle container
    particles = std::unique_ptr<AmrTParContType>(new AmrTParContType(amr_core));
    particles->SetVerbose(verbosity);

    if (make_particles) {
        // now make the particles

        constexpr int level = 0;

        const Geometry& geom = particles->Geom(level);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        IntVect pos_idx;

        Vector<int> done(initial_positions.size(), 0);

        // iterate over all of the boxes on this level and make particles if they fit into one
        for(MFIter mfi = particles->MakeMFIter(level); mfi.isValid(); ++mfi) {
            Box box = mfi.validbox();

            // get the tile of particles for the local box
            TParTileType& pc = particles->DefineAndReturnParticleTile(level, mfi.index(), mfi.LocalTileIndex());

            for (int pi=0; pi < initial_positions.size(); ++pi) {
                if (done[pi]) continue;

                RealArray& pos = initial_positions[pi];

                // convert position to index
                for (int dim=0; dim<AMREX_SPACEDIM; ++dim) {
                    pos_idx[dim] = std::floor((pos[dim] - plo[dim])*dxi[dim]);
                }

                if (box.contains(pos_idx)) {
                    TParticleType p;
                    p.id()   = TParticleType::NextID();
                    p.cpu()  = ParallelDescriptor::MyProc();
                    AMREX_D_TERM(
                                p.pos(0) = pos[0];,
                            p.pos(1) = pos[1];,
                    p.pos(2) = pos[2];
                    )

                    AMREX_D_TERM(
                                p.rdata(+ParticleIdxR::VX) = 0.0;,
                            p.rdata(+ParticleIdxR::VY) = 0.0;,
                    p.rdata(+ParticleIdxR::VZ) = 0.0;
                    )

                    pc.push_back(p);

                    done[pi] = 1;
                }
            }
        }

        particles->Redistribute();
    }
}

void TracerParticle::checkpoint(const std::string& dir)
{
    particles->Checkpoint(dir, "Particles_"+name, true, particle_real_names, particle_int_names);
}

void TracerParticle::restart(const std::string& dir)
{
    particles->Restart(dir, "Particles_"+name);
}

void TracerParticle::redistribute(int level, int finest_level, int ngrow, int local)
{
    particles->Redistribute(level, finest_level, ngrow, local);
}

void TracerParticle::clear()
{
    particles.reset();
}

// TODO: more advanced velocity interpolation!!
void TracerParticle::push_particles(const int level,
                                    MFIter& mfi,
                                    const FArrayBox& prim,
                                    const Geometry geom,
                                    const Real dt
                                    #ifdef AMREX_USE_EB
                                    ,const EBCellFlagFab& flag
                                    #endif
                                    ) const
{
    BL_PROFILE("TracerParticle::push_particles");

    TParTileType& ptile = particles->ParticlesAt(level, mfi);

    // advance particle locations

    //    const Real          strttime = amrex::second();
    const auto          plo      = geom.ProbLoArray();
    const auto          dxi      = geom.InvCellSizeArray();

    for (int ipass = 0; ipass < 2; ipass++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif


        auto& aos  = ptile.GetArrayOfStructs();
        const int n          = aos.numParticles();


        const auto p4 = prim.array();
        auto  p_pbox = aos().data();

#ifdef AMREX_USE_EB

        const auto f4 = flag.array();
#endif
        amrex::ParallelFor(n,
                           [=] AMREX_GPU_DEVICE (int i)
        {
            TParticleType& p  = p_pbox[i];
            if (p.id() <= 0) return;
            Real v[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0,0.0,0.0)};

            // implement particle boundary conditions here??

            // calculate where we are in index space and cell local space
            Array<int,3> iloc = {0,0,0};

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                const Real loc = (p.pos(d) - plo[d]) * dxi[d];
                iloc[d] = static_cast<int>(amrex::Math::floor(loc));
            }
#ifdef AMREX_USE_EB
            if (f4(iloc[0], iloc[1], iloc[2]).isCovered()) {
                p.id() = -p.id();
                return;
            }
#endif

            for (int d1=0; d1<AMREX_SPACEDIM; ++d1) {
                v[d1] = p4(iloc[0], iloc[1], iloc[2],d1); // cell centre value
            }

            // apply updates to particle position and velocity

            if (ipass == 0) {
                for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                    p.rdata(dim) = p.pos(dim);
                    p.pos(dim) += 0.5*dt*v[dim];
                }
            } else  {
                for (int dim=0; dim < AMREX_SPACEDIM; dim++) {
                    p.rdata(dim) = p.rdata(dim) + dt*v[dim];
                    p.rdata(dim) = v[dim];
                }
            }
        });
    }

    return;
}

Vector<std::string> TracerParticle::get_plot_output_names() const
{
    Vector<std::string> plot_names = {"count", "x_vel", "y_vel", "z_vel"};

    return plot_names;
}

void TracerParticle::get_plot_output(MFP* mfp, MultiFab &plot_data, std::vector<std::string> &plot_names) const
{

    std::pair<bool,int> get_count = findInVector(plot_names, "count-"+name);
    std::pair<bool,int> get_u = findInVector(plot_names, "x_vel-"+name);
    std::pair<bool,int> get_v = findInVector(plot_names, "y_vel-"+name);
    std::pair<bool,int> get_w = findInVector(plot_names, "z_vel-"+name);

    if (!(get_count.first || get_u.first || get_v.first ||get_w.first)) return;

    // zero out the data
    if (get_count.first) plot_data.setVal(0.0, get_count.second, 1, 0);
    if (get_u.first) plot_data.setVal(0.0, get_u.second, 1, 0);
    if (get_v.first) plot_data.setVal(0.0, get_v.second, 1, 0);
    if (get_w.first) plot_data.setVal(0.0, get_w.second, 1, 0);

    const auto& P = particles;
    const int level = mfp->get_level();
    const Geometry& geom = mfp->Geom();

    const auto plo = geom.ProbLoArray();
    const Real* dxi = geom.InvCellSize();

    // just do a simple count of how many particles in each cell
    for(amrex::MFIter mfi= P->MakeMFIter(level) ;mfi.isValid();++mfi) {

        const Box& box = mfi.fabbox();

        FArrayBox& plot = plot_data[mfi];

        Array4<Real> const& plot4 = plot.array();

        // Each grid,tile has a their own local particle container
        TParTileType& pc = P->DefineAndReturnParticleTile(level, mfi.index(), mfi.LocalTileIndex());

        auto& aos  = pc.GetArrayOfStructs();
        const int n  = aos.numParticles();

        auto  p_pbox = aos().data();

        IArrayBox counter(box);
        counter.setVal(0);

        Array4<int> const& counter4 = counter.array();

        amrex::ParallelFor(n,[=] AMREX_GPU_DEVICE (int i)
        {
            TParticleType& p  = p_pbox[i];
            if (p.id() <= 0) return;

            // calculate where we are in index space and cell local space
            Array<Real,3> loc = {0.0,0.0,0.0};
            Array<int,3> iloc = {0,0,0};

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                loc[d] = (p.pos(d) - plo[d]) * dxi[d];
                iloc[d] = static_cast<int>(amrex::Math::floor(loc[d]));
            }

            // accumulate an independent counter for later use
            counter4(iloc[0], iloc[1], iloc[2]) += 1;

            if (get_count.first) plot4(iloc[0], iloc[1], iloc[2], get_count.second) += 1;
            if (get_u.first) plot4(iloc[0], iloc[1], iloc[2], get_u.second) += p.rdata(+ParticleIdxR::VX);
            if (get_v.first) plot4(iloc[0], iloc[1], iloc[2], get_v.second) += p.rdata(+ParticleIdxR::VY);
            if (get_w.first) plot4(iloc[0], iloc[1], iloc[2], get_w.second) += p.rdata(+ParticleIdxR::VZ);


        });

        // average
        ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {

            const int N = counter4(i,j,k);

            if (N > 0) {
                if (get_u.first) plot4(i,j,k, get_u.second) /= counter4(i,j,k);
                if (get_v.first) plot4(i,j,k, get_v.second) /= counter4(i,j,k);
                if (get_w.first) plot4(i,j,k, get_w.second) /= counter4(i,j,k);
            }

        });

//        plot_FAB_2d(plot,get_count.second,"count", false, true);
    }
}

void TracerParticle::write_info(nlohmann::json& js) const
{
    LagrangianState::write_info(js);

    js["state_idx"] = state_idx;
    js["type"] = tag;

}
#endif
