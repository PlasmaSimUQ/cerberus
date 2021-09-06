#ifdef AMREX_PARTICLES
#include "MFP_chargedparticle.H"
#include "MFP_state.H"

Vector<std::string> ChargedParticle::particle_real_names = {"mass", "charge", "x_vel", "y_vel", "z_vel"};
Vector<std::string> ChargedParticle::particle_int_names = {};

std::string ChargedParticle::tag = "charged_particle";
bool ChargedParticle::registered = GetStateFactory().Register(ChargedParticle::tag, StateBuilder<ChargedParticle>);


ChargedParticle::ChargedParticle()
{

}

ChargedParticle::ChargedParticle(const sol::table& def)
{
    name = def["name"];
    global_idx = def["global_idx"];
    verbosity = def.get_or("verbosity",0);

    // grab the initial positions
    sol::table particle_def = def["particles"];

    ClassFactory<ParticleDistribution> dfact = GetParticleDistributionFactory();

    std::string distribution_type = particle_def["type"].get<std::string>();

    initial = dfact.Build(distribution_type, particle_def);

    if (!initial)
        Abort("Failed to read particle distribution definition for '"+name+"', must be one of "+vec2str(dfact.getKeys()));

}

void ChargedParticle::init_data(MFP* mfp, const Real time)
{
    if (mfp->get_level() == 0)
        init(mfp);
}

// this should only be called at level 0
void ChargedParticle::init(MFP* mfp, bool make_particles)
{

    AmrCore* amr_core = mfp->get_parent();

    // generate the particle container
    particles = std::unique_ptr<AmrCParContType>(new AmrCParContType(amr_core));
    particles->SetVerbose(verbosity);

    if (make_particles) {
        // now make the particles

        constexpr int level = 0;

        const Geometry& geom = particles->Geom(level);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        IntVect pos_idx;

        switch (initial->get_type()) {
        case ParticleDistributionInit::Defined :
        {
            DefinedSpecies* d_species = static_cast<DefinedSpecies*>(initial.get());

            int np = d_species->num_particles();

            Vector<int> done(np, 0);

            // iterate over all of the boxes on this level and make particles if they fit into one
            for(MFIter mfi = particles->MakeMFIter(level); mfi.isValid(); ++mfi) {
                Box box = mfi.validbox();

                // get the tile of particles for the local box
                CParTileType& pc = particles->DefineAndReturnParticleTile(level, mfi.index(), mfi.LocalTileIndex());

                for (int pi=0; pi < np; ++pi) {
                    if (done[pi]) continue;

                    RealArray& pos = d_species->pos[pi];

                    // convert position to index
                    for (int dim=0; dim<AMREX_SPACEDIM; ++dim) {
                        pos_idx[dim] = std::floor((pos[dim] - plo[dim])*dxi[dim]);
                    }

                    if (box.contains(pos_idx)) {
                        Array<Real,3>& vel = d_species->vel[pi];
                        CParticleType p;
                        p.id()   = CParticleType::NextID();
                        p.cpu()  = ParallelDescriptor::MyProc();
                        AMREX_D_TERM(
                                    p.pos(0) = pos[0];,
                                p.pos(1) = pos[1];,
                        p.pos(2) = pos[2];
                        )

                        p.rdata(+ParticleIdxR::VX) = vel[0];
                        p.rdata(+ParticleIdxR::VY) = vel[1];
                        p.rdata(+ParticleIdxR::VZ) = vel[2];

                        p.rdata(+ParticleIdxR::Charge) = d_species->charge;
                        p.rdata(+ParticleIdxR::Mass) = d_species->mass;

                        pc.push_back(p);

                        done[pi] = 1;
                    }
                }
            }
        }

            break;
        case ParticleDistributionInit::Distribution :
        {
            DistributionSpecies* d_species = static_cast<DistributionSpecies*>(initial.get());

            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_real_distribution<double> dist(-0.5,0.5);

            const Real* prob_lo = geom.ProbLo();
            const Real* dx = geom.CellSize();
            const Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);

            for(amrex::MFIter mfi= particles->MakeMFIter(0); mfi.isValid();++mfi){

#ifdef AMREX_USE_EB
                const EBCellFlagFab& flag = mfp->get_eb_data(global_idx).flags[mfi];
                Array4<const EBCellFlag> const& flag4 = flag.array();
#endif

                // Each grid,tile has a their own local particle container
                CParTileType& pc = particles->DefineAndReturnParticleTile(level, mfi.index(), mfi.LocalTileIndex());

                auto box=mfi.validbox();

                Real x, y, z;
                const auto lo = amrex::lbound(box);
                const auto hi = amrex::ubound(box);

                for     (int k = lo.z; k <= hi.z; ++k) {
                    z = prob_lo[2] + (k + 0.5)*dx[2];
                    z = AMREX_D_PICK(0,0,z);
                    for   (int j = lo.y; j <= hi.y; ++j) {
                        y = prob_lo[1] + (j + 0.5)*dx[1];
                        y = AMREX_D_PICK(0,y,y);
                        AMREX_PRAGMA_SIMD
                                for (int i = lo.x; i <= hi.x; ++i) {
                            x = prob_lo[0] + (i + 0.5)*dx[0];

#ifdef AMREX_USE_EB
                            const EBCellFlag &cflag = flag4(i,j,k);
                            if (cflag.isCovered()) continue;
#endif

                            Real number_density = d_species->number_density(x,y,z);
                            int num_particles = std::floor(d_species->particles_per_cell(x,y,z));
                            if ((number_density <= 0.0) || (num_particles <= 0)) continue;

                            Real total_mass = number_density*d_species->mass*vol;
                            Real particle_mass = total_mass/num_particles;
                            Real particle_charge = (particle_mass/d_species->mass)*d_species->charge;

                            for (int n=0; n<num_particles; ++n) {
                                const Real x_ = x+dist(mt)*dx[0];
                                const Real y_ = y+dist(mt)*dx[1];
                                const Real z_ = z+dist(mt)*dx[2];

                                Array<Real,3> vel = d_species->velocity->operator ()(x_,y_,z_);

                                CParticleType p;
                                p.id()   = CParticleType::NextID();
                                p.cpu()  = ParallelDescriptor::MyProc();
                                AMREX_D_TERM(p.pos(0) = x_;, p.pos(1) = y_;, p.pos(2) = z_;)
                                p.rdata(+ParticleIdxR::Mass) = particle_mass;
                                p.rdata(+ParticleIdxR::Charge) = particle_charge;
                                p.rdata(+ParticleIdxR::VX) = vel[0];
                                p.rdata(+ParticleIdxR::VY) = vel[1];
                                p.rdata(+ParticleIdxR::VZ) = vel[2];

                                //                                p.idata(+ParticleIdxI::???) = ???;

                                pc.push_back(p);
                            }

                        }
                    }
                }
            }
        }
            break;
        default:
            Abort("How did we get here?");
        }

        particles->Redistribute();
    }
}

Real ChargedParticle::get_allowed_time_step(MFP* mfp) const
{
    BL_PROFILE("ChargedParticle::get_allowed_time_step");

    Real dt = std::numeric_limits<Real>::max();

//    const Real* dx = mfp->Geom().CellSize();

//    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);



//    for (CParIterType pti(*(particles.get()), mfp->get_level()); pti.isValid(); ++pti) {
//        Real wt = ParallelDescriptor::second();

//        // Each grid,tile has a their own local particle container
//        CParticleType *  AMREX_RESTRICT p = &(pti.GetArrayOfStructs()[0]);
//        const int np = pti.numParticles();

//        for (int i=0; i<np; ++i)
//        {
//            for (int d = 0; d<AMREX_SPACEDIM; ++d) {
//                const Real vel = particles[i].rdata(+ChargedParticle::ParticleIdxR::VX+d);
//                dt = std::min(dt, dx[d]/vel);
//            }
//        }

//        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
//        cost[pti].plus(wt, box);
//    }


    return dt;
}

void ChargedParticle::checkpoint(const std::string& dir)
{
    particles->Checkpoint(dir, "Particles_"+name, true, particle_real_names, particle_int_names);
}

void ChargedParticle::restart(const std::string& dir)
{
    particles->Restart(dir, "Particles_"+name);
}

void ChargedParticle::redistribute(int level, int finest_level, int ngrow, int local)
{
    particles->Redistribute(level, finest_level, ngrow, local);
}

void ChargedParticle::clear()
{
    particles.reset();
}

void ChargedParticle::push_particles(MFIter& mfi,
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
                                     )
{
    BL_PROFILE("TracerParticle::push_particles");

    // advance particle locations

    // grab the tile of particles

    CParTileType& ptile = particles->ParticlesAt(level, mfi);

    //    const Real          strttime = amrex::second();
    const auto          plo      = geom.ProbLoArray();
    const auto          dxi      = geom.InvCellSizeArray();


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif


    auto& aos  = ptile.GetArrayOfStructs();
    const int n  = aos.numParticles();

    // get the left and right cell values
    Array<Array4<const Real>, AMREX_SPACEDIM> rlo4;
    Array<Array4<const Real>, AMREX_SPACEDIM> rhi4;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        rlo4[d] = rlo[d].array();
        rhi4[d] = rhi[d].array();
    }

    const auto p4 = prim.array();
    auto  p_pbox = aos().data();

#ifdef AMREX_USE_EB

    const auto f4 = flag.array();
#endif

    Real scale_factor = 0.5*dt/MFP::Larmor;

    amrex::ParallelFor(n,
                       [=] AMREX_GPU_DEVICE (int i)
    {
        CParticleType& p  = p_pbox[i];
        if (p.id() <= 0) return;
        Real E[3] = {0.0,0.0,0.0};
        Real B[3] = {0.0,0.0,0.0};

        // implement particle boundary conditions here??

        //        int valid = 0;

        // calculate where we are in index space and cell local space
        Array<Real,3> loc = {0.0,0.0,0.0};
        Array<int,3> iloc = {0,0,0};

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            loc[d] = (p.pos(d) - plo[d]) * dxi[d] - 0.5;
            iloc[d] = static_cast<int>(amrex::Math::floor(loc[d]));
            loc[d] -= iloc[d];
        }
#ifdef AMREX_USE_EB
        if (f4(iloc[0], iloc[1], iloc[2]).isCovered()) {
            p.id() = -p.id();
            return;
        }
#endif

        // get the E and B fields at the particle according to the local slopes
        // obtained via the reconstructed cell face values

        for (int d1=0; d1<3; ++d1) {
            E[d1] = p4(iloc[0], iloc[1], iloc[2], E_idx + d1); // cell centre value
            B[d1] = p4(iloc[0], iloc[1], iloc[2], B_idx + d1); // cell centre value
            for (int d2=0; d2<AMREX_SPACEDIM; ++d2) {
                // adjustment due to slopes
                const Real E_lo_val = rlo4[d2](iloc[0], iloc[1], iloc[2],E_idx + d1);
                const Real E_hi_val = rhi4[d2](iloc[0], iloc[1], iloc[2],E_idx + d1);
                E[d1] += loc[d2]*(E_hi_val - E_lo_val);

                const Real B_lo_val = rlo4[d2](iloc[0], iloc[1], iloc[2],B_idx + d1);
                const Real B_hi_val = rhi4[d2](iloc[0], iloc[1], iloc[2],B_idx + d1);
                B[d1] += loc[d2]*(B_hi_val - B_lo_val);
            }
        }

        // apply updates to particle position and velocity using the Boris pusher

        Real q = p.rdata(+ParticleIdxR::Charge);
        Real m = p.rdata(+ParticleIdxR::Mass);

        Real E0[3], B0[3], Vn[3], V0[3], V1_[3], V1[3];
        for (int i=0; i<3; ++i) {
            E0[i] = scale_factor*q*E[i]/m;
            B0[i] = scale_factor*q*B[i]/m;
            Vn[i] = p.rdata(+ParticleIdxR::VX+i);
            V0[i] = Vn[i] + E0[i];
        }

        Real B02 = 1 + B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2];

        V1_[0] = 2*(V0[0] + V0[1]*B0[2] - V0[2]*B0[1])/B02;
        V1_[1] = 2*(V0[1] + V0[2]*B0[0] - V0[0]*B0[2])/B02;
        V1_[2] = 2*(V0[2] + V0[0]*B0[1] - V0[1]*B0[0])/B02;

        V1[0] = V0[0] + V1_[1]*B0[2] - V1_[2]*B0[1];
        V1[1] = V0[1] + V1_[2]*B0[0] - V1_[0]*B0[2];
        V1[2] = V0[2] + V1_[0]*B0[1] - V1_[1]*B0[0];

        for (int i=0; i<3; ++i) {
            p.rdata(+ParticleIdxR::VX+i) = V1[i] + E0[i];
        }

        for (int dim=0; dim<AMREX_SPACEDIM; ++dim) {
            p.pos(dim) += dt*p.rdata(+ParticleIdxR::VX+dim);
        }

    });

    return;
}

void ChargedParticle::calculate_source(MFIter& mfi, FArrayBox& S, Geometry& geom, int level) const
{

    const Real* dom_lo = geom.ProbLo();
    const Real* dxi = geom.InvCellSize();

    // zero out first
    S.setVal(0.0);

    const Real factor = MFP::Larmor/(MFP::lightspeed*MFP::Debye*MFP::Debye);

    Array4<Real> const& S4 = S.array();

    CParTileType& pc = particles->DefineAndReturnParticleTile(level, mfi.index(), mfi.LocalTileIndex());
    CParticleType *  AMREX_RESTRICT particle = &(pc.GetArrayOfStructs()[0]);
    const int np = pc.numParticles();

    for (int i=0; i<np; ++i) {

        const Real q = particle[i].rdata(+ParticleIdxR::Charge);

        Array<Real,3> vel;
        vel[0] = particle[i].rdata(+ParticleIdxR::VX);
        vel[1] = particle[i].rdata(+ParticleIdxR::VY);
        vel[2] = particle[i].rdata(+ParticleIdxR::VZ);

        // get the local coordinates of the particle
        Array<int,3> coord = {0,0,0};
        AMREX_D_TERM(
                    coord[0]=floor((particle[i].pos(0) - dom_lo[0])*dxi[0]);,
        coord[1]=floor((particle[i].pos(1) - dom_lo[1])*dxi[1]);,
        coord[2]=floor((particle[i].pos(2) - dom_lo[2])*dxi[2]);
        )

        // calculate the current
        for (int vi=0; vi<3; ++vi) {
            S4(coord[0], coord[1], coord[2], vi) -= factor*q*vel[vi];
        }
    }

    return;
}

#ifdef AMREX_USE_EB
void ChargedParticle::set_eb_bc(const sol::table &bc_def)
{

    std::string bc_type = bc_def.get<std::string>("type");


}
#endif


Vector<std::string> ChargedParticle::get_plot_output_names() const
{
    Vector<std::string> plot_names = {"count", "x_vel", "y_vel", "z_vel"};

    return plot_names;
}

void ChargedParticle::get_plot_output(MFP* mfp, MultiFab &plot_data, std::vector<std::string> &plot_names) const
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

        Array4<Real> const& plot4 = plot_data[mfi].array();

        // Each grid,tile has a their own local particle container
        CParTileType& pc = P->DefineAndReturnParticleTile(level, mfi.index(), mfi.LocalTileIndex());

        auto& aos  = pc.GetArrayOfStructs();
        const int n  = aos.numParticles();

        auto  p_pbox = aos().data();

        IArrayBox counter(mfi.tilebox());
        counter.setVal(0);

        Array4<int> const& counter4 = counter.array();

        amrex::ParallelFor(n,[=] AMREX_GPU_DEVICE (int i)
        {
            CParticleType& p  = p_pbox[i];
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
    }
}


void ChargedParticle::write_info(nlohmann::json& js) const
{
    LagrangianState::write_info(js);

    js["state_idx"] = state_idx;
    js["type"] = tag;

}

#endif
