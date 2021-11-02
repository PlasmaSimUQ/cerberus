
#include <AMReX_ParmParse.H>
#include <algorithm>
#include <utility>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

#include "MFP.H"
#include "MFP_state.H"
#include "MFP_eulerian.H"
#include "MFP_lagrangian.H"
#include "MFP_utility.H"
#include "MFP_diagnostics.H"
#include "json.hpp"




using namespace amrex;

std::string MFP::getVersion()
{
    return CERBERUS_GIT_NAME + std::string("_") + CERBERUS_GIT_VERSION;
}

#ifdef AMREX_USE_EB
bool MFP::check_covered_stencil(Array4<const EBCellFlag> const& flag, int i, int j, int k, int d, int stencil_length)
{
    BL_PROFILE("State::check_covered_stencil");
    Array<int,3> stencil_index;
    int offset = stencil_length/2;
    // cell that references a covered cell doesn't need calculating
    stencil_index.fill(0);
    for (int s=0; s<stencil_length; ++s) {
        stencil_index[d] = s - offset;
        // check if any of the stencil values are from a covered cell
        if (flag(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2]).isCovered()) {
            return true;
        }
    }
    return false;
}
#endif

void MFP::calc_slope(const Box& box,
                     const FArrayBox& src,
                     FArrayBox &slope,
                     #ifdef AMREX_USE_EB
                     const EBCellFlagFab &flag,
                     #endif
                     const Real *dx,
                     int index,
                     int dim,
                     Reconstruction &reco)
{

    BL_PROFILE("State::calc_slope");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const Real> const& src4 = src.array();

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& f4 = flag.array();
    bool check_eb = flag.getType() != FabType::regular;

#endif

    Vector<Real> stencil(reco.stencil_length);
    int offset = reco.stencil_length/2;
    Array<int,3> stencil_index;



    // make sure our array is the corect size
    slope.resize(box);

    Array4<Real> const& s4 = slope.array();



    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {


                stencil_index.fill(0);

#ifdef AMREX_USE_EB
                if (check_eb) {
                    // covered cell doesn't need calculating
                    if (f4(i,j,k).isCovered() || check_covered_stencil(f4, i, j, k, dim, reco.stencil_length)) {
                        s4(i,j,k) = 0.0;
                        continue;
                    }
                }
#endif

                for (int s=0; s<reco.stencil_length; ++s) {
                    stencil_index[dim] = s - offset;
                    stencil[s] = src4(i+stencil_index[0], j+stencil_index[1], k+stencil_index[2], index);
                }

                // perform reconstruction
                s4(i,j,k) = reco.get_slope(stencil)/dx[dim]; // account for local cell size
            }
        }
    }


    return;
}



/* load all of the data that we want to plot into a single MultiFAB
 *  and provide the name of that data.
 */
void MFP::getPlotData(MultiFab &plot_data,
                      std::vector<std::string> &plot_names)
{
    BL_PROFILE("MFP::getPlotData");
    int nv = 0; // specified variables
    int nf = plot_functions.size(); // user defined functions
    int N = 0; // total outputs

    for (const auto& f : plot_variables) {
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
        plot_names[nv+fi] = plot_functions[fi].first;
    }

    plot_data.define(grids, dmap, N, 0);
    plot_data.setVal(0.0);

    const Real* dx = geom.CellSize();
    const Real* prob_lo =  geom.ProbLo();
    const Real time = parent->cumTime();

    // get all the data
    Vector<MultiFab> U(states.size());

    for (int data_idx=0; data_idx<eulerian_states.size(); ++data_idx) {

        EulerianState &istate = EulerianState::get_state(data_idx);

        // get a full array of data at this level
        int ns = desc_lst[data_idx].nComp();
        int ng = istate.get_num_grow();
        U[data_idx].define(grids, dmap, ns, ng, MFInfo(),Factory());

#ifdef AMREX_USE_EB
        EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate.eb2_index));
#endif
        FillPatch(*this, U[data_idx], ng, time, data_idx, 0, ns, 0);
    }


    std::map<std::string,Real> dat = {{"Larmor",Larmor}, {"Debye",Debye},
                                      {"c",lightspeed}, {"skin_depth",skin_depth},
                                      {"beta",beta},{"t",time}};


    std::map<std::string,FArrayBox> dat_arrays;
    std::map<std::string,Array4<Real>> dat4_arrays;

    const Array<std::string, AMREX_SPACEDIM> grad = {AMREX_D_DECL("-dx","-dy","-dz")};
    Vector<std::string> updated;

    MultiFab& cost = get_data(Cost_Idx, time);


    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {
        const Box& orig_box = mfi.tilebox();

        for (const auto& pair : plot_variables) {
            dat_arrays[pair.first].resize(orig_box,1);
        }

        // cost
        const bool get_cost = plot_variables.find("cost") != plot_variables.end();
        if (get_cost) {
            dat_arrays["cost"].copy(cost[mfi]);
        }

        for (int data_idx=0; data_idx<eulerian_states.size(); ++data_idx) {
            EulerianState &istate = EulerianState::get_state(data_idx);

            const Box box = grow(orig_box,istate.get_num_grow());

#ifdef AMREX_USE_EB
            EBData& eb = get_eb_data(istate.global_idx);
            const EBCellFlagFab& flag = eb.flags[mfi];
            const FArrayBox& vfrac = eb.volfrac[mfi];
#endif

            // raw values
            istate.get_plot_output(box,
                                   U[data_idx][mfi],
                                   dat_arrays,
                                   updated
                       #ifdef AMREX_USE_EB
                                   ,vfrac
                       #endif
                                   );

            // get any gradients

            for (const std::string& var_name : updated) {

                const Array<int,AMREX_SPACEDIM+1> var_grad = plot_variables[var_name];

                for (int i=0; i<AMREX_SPACEDIM; ++i) {
                    if (!var_grad[i+1]) continue;

                    // this should really go in the startup routine in global_data.cpp
                    // but it would require a not insignificant amount of work, so here it is
                    if (istate.reconstructor->stencil_length < 2) {
                        ClassFactory<Reconstruction> rfact = GetReconstructionFactory();
                        std::string msg = "";
                        msg += "Gradient requested for "+var_name+" but state '";
                        msg += istate.name;
                        msg += "' has an invalid reconstruction '";
                        msg += istate.reconstructor->get_tag();
                        msg += "' (stencil length < 2). Options are ";
                        msg += vec2str(rfact.getKeys());
                        Abort(msg);
                    }

                    const std::string grad_name = var_name + grad[i];

                    dat_arrays[grad_name].resize(orig_box,1);

                    calc_slope(orig_box,
                               dat_arrays[var_name],
                               dat_arrays[grad_name],
           #ifdef AMREX_USE_EB
                               flag,
           #endif
                               dx, 0, i,
                               *istate.reconstructor.get());
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
                    for (const auto& var : plot_variables) {
                        const std::string name = var.first;
                        const int index = var.second[0];

                        Real value = dat4_arrays[name](i,j,k);

                        dat[name] = value;
                        if (index > -1) pd4(i,j,k,index) = value;
                    }

                    // call the function
                    for (int fi=0; fi<nf; ++fi) {
                        pd4(i,j,k,fi+nv) = plot_functions[fi].second(dat);
                    }
                }
            }
        }

    }

    //    for (int i=0; i<plot_data.nComp(); ++i) {
    //        plot_FAB_2d(plot_data, i, 0, plot_names[i], false, true);
    //    }

#ifdef AMREX_PARTICLES
    // now do the lagrangian states
    for (const int& global_idx : lagrangian_states) {
        LagrangianState& istate = LagrangianState::get_state_global(global_idx);

        istate.get_plot_output(this, plot_data, plot_names);
    }
#endif

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
    for (const auto& i : lagrangian_states) {
        LagrangianState& istate = LagrangianState::get_state_global(i);
        std::string particle_folder = "Particles_"+istate.name;
        cmds.push_back("\\cd "+dir+";\\tar -cf " + particle_folder + ".tar " + particle_folder);
        to_remove.push_back(dir+"/"+particle_folder);
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
    write_info(mfp);

    for (const auto& state : states) {
        auto& grp = mfp["state_"+num2str(state->global_idx)];
        state->write_info(grp);
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

    writePlotFilePost(dir, os);

    if (archive_checkpoint)
        archive_folder(dir);
}

void MFP::write_info(nlohmann::json &js) const
{
    BL_PROFILE("GlobalData::write_info");
    // write out globally defined data

    js["num_state"] = states.size();
    js["x_ref"] = x_ref;
    js["n_ref"] = n_ref;
    js["m_ref"] = m_ref;
    js["rho_ref"] = rho_ref;
    js["T_ref"] = T_ref;
    js["u_ref"] = u_ref;
    js["n0"] = n0;
    js["lightspeed"] = lightspeed;
    js["beta"] = beta;
    js["skin_depth"] = skin_depth;
    js["Larmor"] = Larmor;
    js["Debye"] = Debye;

}
