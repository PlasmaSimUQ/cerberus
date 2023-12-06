#ifdef AMREX_USE_EB
    #include "MFP.H"
    #include "MFP_eb_sdf.H"
    #include "MFP_state.H"

    #include <AMReX_EB2.H>
    #include <AMReX_EB2_IF.H>

using namespace amrex;

void initialize_EB2(Amr& amr, int ngrow)
{
    // we are running an explicit time stepper where coarse cell values are over-written
    // by fine cells and thus the coarse cell does not need to be calculated exactly
    // see AMReX issue : https://github.com/AMReX-Codes/amrex/issues/608
    constexpr bool build_coarse_level_by_coarsening = false;

    Geometry geom = amr.Geom(amr.maxLevel());
    const int required_coarsening_level = amr.maxLevel();
    const int max_coarsening_level = amr.maxLevel();

    Vector<FunctionIF> funcs(MFP::states.size());

    if (!MFP::eb_def.empty()) {
        // make all of the individual geometry and collate the geometry functions
        // per state

        for (auto& eb : MFP::eb_def) {
            FunctionIF lua_geom(&eb.geom_func, eb.inside);
            EB2::GeometryShop<FunctionIF> gshop(lua_geom);
            EB2::Build(gshop,
                       geom,
                       required_coarsening_level,
                       max_coarsening_level,
                       ngrow,
                       build_coarse_level_by_coarsening);
            eb.index_space = EB2::TopIndexSpace();

            for (const auto& si : eb.states) {
                if (eb.insertion_type > 0) {
                    funcs[si.first].plus(&eb.geom_func, eb.inside);
                } else {
                    funcs[si.first].minus(&eb.geom_func, eb.inside);
                }
            }
        }
    }

    // now make a per-state geometry
    const EB2::IndexSpace* all_regular = nullptr;
    for (int idx = 0; idx < MFP::states.size(); ++idx) {
        State& istate = MFP::get_state(idx);

        if (funcs[idx].empty()) {
            if (!all_regular) {
                EB2::AllRegularIF rif;
                EB2::GeometryShop<EB2::AllRegularIF> gshop(rif);
                EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level, ngrow);
                all_regular = EB2::TopIndexSpace();
            }
            istate.eb2_index = all_regular;
            istate.eb_all_regular = true;
        } else {
            EB2::GeometryShop<FunctionIF> gshop(funcs[idx]);
            EB2::Build(gshop,
                       geom,
                       required_coarsening_level,
                       max_coarsening_level,
                       ngrow,
                       build_coarse_level_by_coarsening);
            istate.eb2_index = EB2::TopIndexSpace();
            istate.eb_all_regular = false;
        }
    }
}
#endif
