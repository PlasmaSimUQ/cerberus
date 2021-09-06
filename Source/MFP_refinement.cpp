#include "MFP.H"
#include "MFP_eulerian.H"

void MFP::errorEst(TagBoxArray& tags, int, int, Real time, int, int) {
    BL_PROFILE("MFP::errorEst()");

    const int tagval = TagBox::SET;
    const int clearval = TagBox::CLEAR;


    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    for (const auto& istate : states) {
        istate->get_refinement_tags(this, tags);
    }

    if (!refine_boxes.empty()) {
        const Real* problo = geom.ProbLo();
        const Real* dx = geom.CellSize();

        const MultiFab& cost = get_new_data(Cost_Idx);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(cost); mfi.isValid(); ++mfi) {
            auto& fab = tags[mfi];
            const Box& bx = mfi.tilebox();
            for (BoxIterator bi(bx); bi.ok(); ++bi) {
                const IntVect& cell = bi();
                RealVect pos{AMREX_D_DECL((cell[0] + 0.5) * dx[0] + problo[0],
                            (cell[1] + 0.5) * dx[1] + problo[1],
                            (cell[2] + 0.5) * dx[2] + problo[2])};

                int set_type = -1;

                if (only_refine_in_box)
                    set_type = TagBox::CLEAR;

                for (int bidx=0; bidx < refine_boxes.size(); ++bidx) {

                    const auto& rbx = refine_boxes[bidx];

                    if (rbx.contains(pos)) {
                        if (refine_box_type[bidx] == RefineBoxType::OnlyRefine) {
                            set_type = fab(cell); // do nothing
                        } else if (refine_box_type[bidx] == RefineBoxType::NoRefine) {
                            set_type = std::max((int)TagBox::CLEAR, set_type); // clear, but only if we haven't already tagged for refinement
                        } else if (refine_box_type[bidx] == RefineBoxType::ForceRefine) {
                            set_type = std::max((int)TagBox::SET, set_type); // tag for refinement
                        }
                    }
                }

                if (set_type > -1)
                    fab(cell) = set_type;
            }
        }
    }
}
