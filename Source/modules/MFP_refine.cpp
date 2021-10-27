#include "MFP_refine.H"
#include "MFP_state.H"

Refinement::Refinement(){}

#ifdef AMREX_USE_EB
void Refinement::tag_cut_cells (MFP* mfp, TagBoxArray& tags, const int idx)
{
    BL_PROFILE("Refinement::tag_cut_cells");

    const char   tagval = TagBox::SET;

    auto const& flags = mfp->get_eb_data(idx).flags;

    const auto &cost = mfp->get_new_data(mfp->Cost_Idx);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cost, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const auto& flag = flags[mfi];

        const FabType typ = flag.getType(bx);
        if (typ != FabType::regular && typ != FabType::covered)
        {
            Array4<char> const& tagarr = tags.array(mfi);
            Array4<EBCellFlag const> const& flagarr = flag.const_array();
            AMREX_HOST_DEVICE_FOR_3D (bx, i, j, k,
            {
                                          if (flagarr(i,j,k).isSingleValued()) {
                                              tagarr(i,j,k) = tagval;
                                          }
                                      });
        }
    }
}
#endif
