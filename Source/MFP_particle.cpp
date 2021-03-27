#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include "MFP.H"

#ifdef AMREX_PARTICLES
#include "AMReX_AmrParticles.H"
#endif

using namespace amrex;

#ifdef AMREX_PARTICLES

Vector<std::string> MFP::particle_names;
Vector<int> MFP::particle_idx;
Vector<amrex::AmrTracerParticleContainer*> MFP::particles;


void
MFP::read_particle_params ()
{

    ParmParse pp("particles");
    pp.query("do_particles", gd.do_tracer_particles);
    pp.query("v", gd.particle_verbose);

    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

}

void
MFP::init_particles ()
{
    BL_PROFILE("MFP::init_particles()");

    if (level > 0)
        return;

    if (gd.do_tracer_particles) {
        for (int idx=0; idx<gd.num_solve_state; ++idx) {

            State &istate = gd.get_state(idx);

            if (!istate.particle_init.empty()) {
                AmrTracerParticleContainer* TracerPC = new AmrTracerParticleContainer(parent);
                TracerPC->SetVerbose(gd.particle_verbose);

                TracerPC->InitFromAsciiFile(istate.particle_init,0);

                istate.particle_index = particles.size();
                particles.push_back(TracerPC);
                particle_idx.push_back(idx);
                particle_names.push_back(gd.state_names[idx]);

            }




        }
    }
}

void MFP::writeParticles(const std::string& dir)
{

    if (gd.do_tracer_particles) {

        for (int idx=0; idx<particles.size(); ++idx) {
#ifdef AMREX_USE_HDF5
            particles[idx]->CheckpointHDF5 (dir, "Particles_"+particle_names[idx], true, {AMREX_D_DECL("vel-x", "vel-y", "vel-z")});
#else
            particles[idx]->Checkpoint(dir,"Particles_"+particle_names[idx], true, {AMREX_D_DECL("vel-x", "vel-y", "vel-z")});
#endif
        }
    }
}

void
MFP::ParticlePostRestart (const std::string& dir)
{
    if ((level == 0) && gd.do_tracer_particles) {

        // handle if we have archived level data in a restart folder
        Vector<std::string> to_remove;
        if (ParallelDescriptor::IOProcessor()) {

            // get the names of the particle folders
            for (int idx=0; idx<gd.num_solve_state; ++idx) {

                State &istate = gd.get_state(idx);
                if (!istate.particle_init.empty()) {
                    std::string &particle_name = gd.state_names[idx];

                    std::string particle_folder = "Particles_"+particle_name;
                    std::string FullPath = dir+"/"+particle_folder;

                    // check if we actually need to un-tar
                    if (FileExists(FullPath))
                        continue;

                    // check tar exists
                    if (!FileExists(FullPath + ".tar"))
                        continue;

                    // add to list of things to clean up later
                    to_remove.push_back(FullPath);


                    // perform un-tar
                    std::string cmd = "\\tar -C " + dir + " -xf " + FullPath + ".tar ";
                    const char * command = {cmd.c_str()};
                    int retVal = std::system(command);
                    if (retVal == -1 || WEXITSTATUS(retVal) != 0) {
                        Abort("Error: Unable to un-tar '"+FullPath+"'");
                    }
                }
            }
        }

        ParallelDescriptor::Barrier();

        // retrieve particle data
        for (int idx=0; idx<gd.num_solve_state; ++idx) {

            State &istate = gd.get_state(idx);
            if (!istate.particle_init.empty()) {

                AmrTracerParticleContainer* TracerPC = new AmrTracerParticleContainer(parent);
                particles.push_back(TracerPC);
                particle_idx.push_back(idx);
                particle_names.push_back(gd.state_names[idx]);

                TracerPC->SetVerbose(gd.particle_verbose);

#ifdef AMREX_USE_HDF5
                TracerPC->RestartHDF5(dir, "Particles_"+gd.state_names[idx]);
#else
                TracerPC->Restart(dir, "Particles_"+gd.state_names[idx]);
#endif
            }
        }

        // clean up after un-tar operation
        if (ParallelDescriptor::IOProcessor()) {
            for (const auto& path : to_remove) {
                FileSystem::RemoveAll(path);
            }
        }
    }
}

//
// Uses midpoint method to advance particles using cell-centered velocity
//

void
MFP::push_particles (ParTileType& ptile,
                     const FArrayBox& prim,
                     Array<FArrayBox, AMREX_SPACEDIM> &rlo,
                     Array<FArrayBox, AMREX_SPACEDIM> &rhi,
                     const int vel_idx,
                     const Real dt
                     EB_OPTIONAL(,const EBCellFlagFab& flag)
                     )
{


    // advance particle locations

    const Real          strttime = amrex::second();
    const auto          plo      = geom.ProbLoArray();
    const auto          dxi      = geom.InvCellSizeArray();

    for (int ipass = 0; ipass < 2; ipass++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif


        auto& aos  = ptile.GetArrayOfStructs();
        const int n          = aos.numParticles();

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
        amrex::ParallelFor(n,
                           [=] AMREX_GPU_DEVICE (int i)
        {
            ParticleType& p  = p_pbox[i];
            if (p.id() <= 0) return;
            Real v[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0,0.0,0.0)};

            // implement particle boundary conditions here??

            int valid = 0;

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

            // get the velocity at the particle according to the local slopes
            // obtained via the reconstructed cell face values

            for (int d1=0; d1<AMREX_SPACEDIM; ++d1) {
                v[d1] = p4(iloc[0], iloc[1], iloc[2],vel_idx + d1); // cell centre value
                for (int d2=0; d2<AMREX_SPACEDIM; ++d2) {
                    // adjustment due to slopes
                    const Real lo_val = rlo4[d2](iloc[0], iloc[1], iloc[2],vel_idx + d1);
                    const Real hi_val = rhi4[d2](iloc[0], iloc[1], iloc[2],vel_idx + d1);
                    v[d1] += loc[d2]*(hi_val - lo_val);
                }
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




#endif
