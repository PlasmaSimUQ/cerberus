

#include <AMReX_ParmParse.H>
#include <MFP.H>
#include <MFP_utility.H>
#include <algorithm>
#include <utility>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

#include "MFP_global.H"
#include "MFP_diagnostics.H"
#include "json.hpp"

#ifdef AMREX_USE_HDF5
#include <hdf5.h>
#endif



using namespace amrex;

std::string MFP::getVersion()
{
    return CERBERUS_GIT_NAME + std::string("_") + CERBERUS_GIT_VERSION;
}




/* load all of the data that we want to plot into a single MultiFAB
 *  and provide the name of that data.
 */
void MFP::getPlotData(MultiFab &plot_data,
                      std::vector<std::string> &plot_names)
{
    BL_PROFILE("MFP::getPlotData");
    int nv = 0; // specified variables
    int nf = gd.plot_functions.size(); // user defined functions
    int N = 0; // total outputs

    for (const auto& f : gd.plot_variables) {
        const int idx = f.second[0];
        if (idx > -1) {
            plot_names.resize(max((int)plot_names.size(), idx+1));
            plot_names[idx] = f.first;
        }
    }

    nv = plot_names.size();
    N = nv + nf;
    plot_names.resize(N);

    for (int fi=0; fi<nf; ++fi) {
        plot_names[nv+fi] = gd.plot_functions[fi].first;
    }

    plot_data.define(grids, dmap, N, 0);

    const Real* dx = geom.CellSize();
    const Real* prob_lo =  geom.ProbLo();
    const Real time = parent->cumTime();

    // get all the data
    Vector<MultiFab> U(gd.num_solve_state);
    for (int idx=0; idx<gd.num_solve_state; ++idx) {
        State &istate = gd.get_state(idx);
        // get a full array of data at this level
        int ns = desc_lst[idx].nComp();
        int ng = istate.num_grow;
        U[idx].define(grids, dmap, ns, ng, MFInfo(),Factory());

#ifdef AMREX_USE_EB
        EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate.eb2_index));
#endif
        FillPatch(*this, U[idx], ng, time, idx, 0, ns, 0);
    }


    std::map<std::string,Real> dat = {{"Larmor",gd.Larmor}, {"Debye",gd.Debye},
                                      {"c",gd.lightspeed}, {"skin_depth",gd.skin_depth},
                                      {"beta",gd.beta},{"t",time}};


    std::map<std::string,FArrayBox> dat_arrays;
    std::map<std::string,Array4<Real>> dat4_arrays;

    const Array<std::string, AMREX_SPACEDIM> grad = {AMREX_D_DECL("-dx","-dy","-dz")};
    Vector<std::string> updated;

    MultiFab& cost = get_new_data(gd.Cost_Idx);

    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {
        const Box& orig_box = mfi.tilebox();

        for (const auto& pair : gd.plot_variables) {
            dat_arrays[pair.first].resize(orig_box,1);
        }

        // cost
        const bool get_cost = gd.plot_variables.find("cost") != gd.plot_variables.end();
        if (get_cost) {
            dat_arrays["cost"].copy(cost[mfi]);
        }

        for (int idx=0; idx<gd.num_solve_state; ++idx) {
            State &istate = gd.get_state(idx);
            const Box box = grow(orig_box,istate.num_grow);

#ifdef AMREX_USE_EB
            const EBCellFlagFab& flag = getEBData(idx).flags[mfi];
            const FArrayBox& vfrac = getEBData(idx).volfrac[mfi];
#endif

            // raw values
            istate.get_state_values(box, U[idx][mfi], dat_arrays, updated EB_OPTIONAL(,vfrac));

            // shock detector
            std::string shock_name = "shock-"+istate.name;
            const bool get_shock = gd.plot_variables.find(shock_name) != gd.plot_variables.end();
            if ((istate.shock_idx > -1) && get_shock) {
                MultiFab& S = get_new_data(gd.Shock_Idx);
                dat_arrays[shock_name].copy(S[mfi], orig_box, SrcComp(istate.shock_idx), DestComp(0), NumComps(1));
            }

            // get any gradients

            for (const std::string& var_name : updated) {

                const Array<int,AMREX_SPACEDIM+1> var_grad = gd.plot_variables[var_name];

                for (int i=0; i<AMREX_SPACEDIM; ++i) {
                    if (!var_grad[i+1]) continue;

                    // this should really go in the startup routine in global_data.cpp
                    // but it would require a not insignificant amount of work, so here it is
                    if (istate.reconstruction->stencil_length < 2) {
                        PhysicsFactory<Reconstruction> rfact = GetReconstructionFactory();
                        std::string msg = "";
                        msg += "Gradient requested for "+var_name+" but state '";
                        msg += istate.name;
                        msg += "' has an invalid reconstruction '";
                        msg += istate.reconstruction->get_tag();
                        msg += "' (stencil length < 2). Options are ";
                        msg += vec2str(rfact.getKeys());
                        Abort(msg);
                    }

                    const std::string grad_name = var_name + grad[i];

                    dat_arrays[grad_name].resize(orig_box,1);

                    State::calc_slope(orig_box,
                                      dat_arrays[var_name],
                                      dat_arrays[grad_name],
                                      EB_OPTIONAL(flag,)
                                      dx, 0, i,
                                      *istate.reconstruction.get());
                }
            }
        }

//        plot_FABs_2d(dat_arrays,0,false,true);

        // now iterate over the data, load it into a map, and call the lua functions on it
        const Dim3 lo = amrex::lbound(orig_box);
        const Dim3 hi = amrex::ubound(orig_box);

        for (const auto& pair : dat_arrays) {
            dat4_arrays[pair.first] = dat_arrays[pair.first].array();
        }

        Array4<Real> const& pd4 = plot_data[mfi].array();

        for     (int k = lo.z; k <= hi.z; ++k) {
            dat["z"] = prob_lo[2] + (k + 0.5)*dx[2];
            for   (int j = lo.y; j <= hi.y; ++j) {
                dat["y"] = prob_lo[1] + (j + 0.5)*dx[1];
                for (int i = lo.x; i <= hi.x; ++i) {
                    dat["x"] = prob_lo[0] + (i + 0.5)*dx[0];

                    // load the data
                    for (const auto& var : gd.plot_variables) {
                        const std::string name = var.first;
                        const int index = var.second[0];

                        Real value = dat4_arrays[name](i,j,k);

                        dat[name] = value;
                        if (index > -1) pd4(i,j,k,index) = value;
                    }

                    // call the function
                    for (int fi=0; fi<nf; ++fi) {
                        pd4(i,j,k,fi+nv) = gd.plot_functions[fi].second(dat);
                    }
                }
            }
        }
    }

    return;
}


void
MFP::writePlotFile (const std::string& dir,
                    std::ostream&      os,
                    VisMF::How         how)
{
    BL_PROFILE("MFP::writePlotFile");
    MultiFab plot_data;
    Vector<std::string> plot_names;

    getPlotData(plot_data, plot_names);

    int n_data_items = plot_data.nComp();



    // get the time from the first State_Type
    // if the State_Type is ::Interval, this will get t^{n+1/2} instead of t^n
    Real cur_time = state[0].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            amrex::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

        //
        // Names of variables
        //
        for (int i =0; i < n_data_items; i++) {
            os << plot_names[i] << '\n';
        }

        os << AMREX_SPACEDIM << '\n';
        os << cur_time << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            os << Geom().ProbLo(i) << ' ';
        os << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            os << Geom().ProbHi(i) << ' ';
        os << '\n';
        for (int i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < AMREX_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geom().Coord() << '\n';
        os << "0\n"; // Write bndry data.

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string sLevel = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if ( ! FullPath.empty() && FullPath[FullPath.size()-1] != '/')
    {
        FullPath += '/';
    }
    FullPath += sLevel;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if ( ! levelDirectoryCreated) {
        if (ParallelDescriptor::IOProcessor()) {
            if ( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
                amrex::CreateDirectoryFailed(FullPath);
            }
        }
        // Force other processors to wait until directory is built.
        ParallelDescriptor::Barrier();
    }

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (int i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (int n = 0; n < AMREX_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = sLevel;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }

#ifdef AMREX_USE_EB
        // volfrac threshhold for amrvis
        if (level == parent->finestLevel()) {
            for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
                os << "1.0e-6\n";
            }
        }
#endif
    }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    if (AsyncOut::UseAsyncOut()) {
        VisMF::AsyncWrite(plot_data,TheFullPath);
    } else {
        VisMF::Write(plot_data,TheFullPath,how,true);
    }

    levelDirectoryCreated = false;  // ---- now that the plotfile is finished
}

void MFP::archive_folder(const std::string &dir)
{
    BL_PROFILE("MFP::archive_folder");
    // make sure everyone has finished writing to the directory
    ParallelDescriptor::Barrier();

    if (!ParallelDescriptor::IOProcessor())
        return;

    Vector<std::string> cmds;
    Vector<std::string> to_remove;

    // get the names of what we want to tar
    std::string LevelDir, FullPath;
    LevelDirectoryNames(dir, LevelDir, FullPath);

    cmds.push_back("\\cd "+dir+";\\tar -cf " + LevelDir + ".tar " + LevelDir);
    to_remove.push_back(FullPath);

#ifdef AMREX_PARTICLES
    if (gd.do_tracer_particles) {

        for (int idx=0; idx<particles.size(); ++idx) {
            std::string particle_folder = "Particles_"+particle_names[idx];
            cmds.push_back("\\cd "+dir+";\\tar -cf " + particle_folder + ".tar " + particle_folder);
            to_remove.push_back(dir+"/"+particle_folder);
        }
    }
#endif

    // perform the archiving operation
    for (const auto &cmd : cmds) {
        const char * command = {cmd.c_str()};
        int retVal = std::system(command);
        if (retVal == -1 || WEXITSTATUS(retVal) != 0) {
            Abort("Error: Unable to tar '"+FullPath+"'");
        }
    }

    // now delete the original data
    for (const auto &rm : to_remove) {
        FileSystem::RemoveAll(rm);
    }

    return;

}


void MFP::writePlotFilePost(const std::string &dir, std::ostream &os)
{
    BL_PROFILE("MFP::writePlotFilePost");
#ifdef AMREX_PARTICLES
    writeParticles(dir);
#endif

    if (level != parent->finestLevel() || !ParallelDescriptor::IOProcessor())
        return;



    // metadata for simulation

    nlohmann::json mfp;


    // write the global variables
    gd.write_info(mfp);

    for (int global_idx = 0; global_idx < gd.num_solve_state; ++global_idx) {
        State &istate = gd.get_state(global_idx);
        auto& grp = mfp["state_"+num2str(global_idx)];
        istate.write_info(grp);

        for (const auto &ode_ptr : gd.ode_source_terms) {
            const ODESystem& ode = *ode_ptr;
            for (const std::unique_ptr<SourceTerm> &src: ode.sources) {
                for (const auto &idx : src->offsets) {
                    if (idx.global == global_idx) {
                        src->write_info(grp);
                    }
                }
            }
        }

        grp["cons_names"] = istate.get_cons_names();

#ifdef AMREX_PARTICLES
        if (gd.do_tracer_particles) {
            for (int pidx=0; pidx<particles.size(); ++pidx) {
                if (particle_idx[pidx] != global_idx)
                    continue;

                grp["particle_data"] = "Particles_"+particle_names[pidx];
            }
        }
#endif
    }

    mfp["version"] = getVersion();
    mfp["AMReX_version"] = amrex::Version();

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());


    // write out the json info

    std::ofstream InfoFile;
    InfoFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    InfoFile.open(dir+"/info", std::ios::out | std::ios::trunc | std::ios::binary);
    if ( ! InfoFile.good()) {
        amrex::FileOpenFailed("info");
    }

    InfoFile << std::setw(4) << mfp << std::endl;

    InfoFile.close();

    // write out the config script

    std::ofstream InputsFile;
    InputsFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    InputsFile.open(dir+"/inputs", std::ios::out | std::ios::trunc | std::ios::binary);
    if ( ! InputsFile.good()) {
        amrex::FileOpenFailed("inputs");
    }

    ParmParse::dumpTable(InputsFile, false);

    InputsFile.close();

}


void MFP::checkPointPost (const std::string& dir,std::ostream& os)
{
    BL_PROFILE("MFP::checkPointPost");
#ifdef AMREX_PARTICLES
    writeParticles(dir);
#endif

    writePlotFilePost(dir, os);

    if (archive_checkpoint)
        archive_folder(dir);
}
