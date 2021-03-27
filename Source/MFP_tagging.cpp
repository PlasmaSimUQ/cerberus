#include "MFP.H"

#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#endif

Real refine_criteria(Array<Real,3> S, Real low_val=1e-12)
{
  Real a = S[2] - 2*S[1] + S[0];
  Real b = std::abs(S[2] - S[1]) + std::abs(S[1] - S[0]) + 0.01*(S[2] + 2*S[1] + S[0]);
  if (std::abs(b) <= low_val) {
    return 0.0;
  } else {
    return std::abs(a/b);
  }
}

void MFP::tag_refinement(const Box& box,
                         const FArrayBox& src,
                         EB_OPTIONAL(const EBCellFlagFab& flags,)
                         TagBox& tags,
                         const int n,
                         const Real threshold,
                         const Real min_val,
                         const char tagval,
                         const char clearval)
{

  const Dim3 lo = amrex::lbound(box);
  const Dim3 hi = amrex::ubound(box);

  Array4<Real const> const& src4 = src.array();
  Array4<char> const& tag = tags.array();

  Real val;
  Array<Real,3> S;
  Array<int,3> index;

#ifdef AMREX_USE_EB
  Array4<const EBCellFlag> const& flag4 = flags.array();
#endif

  for     (int k = lo.z; k <= hi.z; ++k) {
    for   (int j = lo.y; j <= hi.y; ++j) {
      AMREX_PRAGMA_SIMD
      for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
          if (flag4(i,j,k).isCovered())
              continue;
#endif
        val = 0.0;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
          // fill our stencil
          index.fill(0);
          index[d] = 1;
          S[0] = src4(i-index[0], j-index[1], k-index[2], n);
          S[1] = src4(i, j, k, n);
          S[2] = src4(i+index[0], j+index[1], k+index[2], n);
          // calculate the refinement criteria and compare it to previous value
          val = std::max(val, refine_criteria(S, min_val));
        }

        // check against threshold and mark for refinement if necessary
        if (val >= threshold) {
          tag(i,j,k) = tagval;
        }
      }
    }
  }
}

#ifdef AMREX_USE_EB
void MFP::tag_cut_cells (TagBoxArray& tags, const int idx)
{
    const char   tagval = TagBox::SET;

    EBData &eb_data = getEBData(idx);

    auto const& flags = eb_data.flags;

    const auto &state = get_new_data(idx);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
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

void MFP::errorEst(TagBoxArray& tags, int, int, Real time, int, int) {
    BL_PROFILE("MFP::errorEst()");

    const int tagval = TagBox::SET;
    const int clearval = TagBox::CLEAR;


    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    FArrayBox Q;
    FArrayBox *fbox;
    bool has_prim;
    for (int idx = 0; idx < gd.num_solve_state; ++idx) {

        State &istate = *gd.states[idx];

        int n_prim = istate.n_prim();
        int n_cons = istate.n_cons();
        if ((level < istate.refine_grad_max_lev) || (istate.refine_grad_max_lev < 0)) {

            // check if we need the primitive values for this state
            has_prim = false;
            for (const auto& info : istate.refine_grad_threshold) {
                if (std::get<0>(info)){
                    has_prim = true;
                }
            }

            // grab the conservative state
            const int num_grow = 1;
            MultiFab grab(grids, dmap, n_cons, num_grow, MFInfo(),Factory());

#ifdef AMREX_USE_EB
            EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate.eb2_index));
#endif

            FillPatch(*this, grab, num_grow, time, idx, 0, n_cons);

#ifdef AMREX_USE_EB
            auto const& flags = getEBData(idx).flags;
            auto const& volfrac = getEBData(idx).volfrac;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(grab); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
                const FArrayBox& vfrac = volfrac[mfi];
                const EBCellFlagFab& flag = flags[mfi];
                FabType flag_type = flag.getType(bx);
#else
                FabType flag_type = FabType::regular;
#endif

                // check if we can skip an empty box
                if (is_vacuum(bx, grab[mfi], idx)) {
                    flag_type = FabType::covered;
                }

                if (flag_type != FabType::covered) {

                    // only calculate the primitives if we need them
                    if (has_prim) {

                        Box pbx = grow(bx, 1);

                        Q.resize(pbx, n_prim);

                        // get primitives
                        istate.calc_primitives(pbx,
                                               grab[mfi],
                                               Q,
                                               dx,
                                               parent->cumTime(),
                                               prob_lo
                                               EB_OPTIONAL(,vfrac)
                                               );
                    }

                    // now go through the list of things to check for high gradients
                    for (const auto& info : istate.refine_grad_threshold) {
                        const int& icomp = std::get<1>(info); // index
                        const Real& refine_grad = std::get<2>(info); // threshold value
                        const Real& refine_min_val = std::get<3>(info); // minimum value value

                        // select primitive or conserved data
                        if (std::get<0>(info)){
                            fbox = &Q;
                        } else {
                            fbox = &grab[mfi];
                        }

                        tag_refinement(bx, *fbox,
                                       EB_OPTIONAL(flags[mfi],)
                                       tags[mfi],
                                       icomp,
                                       refine_grad,
                                       refine_min_val,
                                       tagval,
                                       clearval);

                    }
                }
            }
        }
    }


#ifdef AMREX_USE_EB
    if (gd.refine_cutcells) {
        for (int idx = 0; idx < gd.num_solve_state; ++idx) {
            tag_cut_cells (tags, idx);
        }
    }
#endif

    if (!gd.refine_boxes.empty()) {
        const Real* problo = geom.ProbLo();
        const Real* dx = geom.CellSize();

        const MultiFab& S_new = get_new_data(0);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            auto& fab = tags[mfi];
            const Box& bx = mfi.tilebox();
            for (BoxIterator bi(bx); bi.ok(); ++bi) {
                const IntVect& cell = bi();
                RealVect pos{AMREX_D_DECL((cell[0] + 0.5) * dx[0] + problo[0],
                            (cell[1] + 0.5) * dx[1] + problo[1],
                            (cell[2] + 0.5) * dx[2] + problo[2])};

                int set_type = -1;

                if (gd.only_refine_in_box)
                    set_type = TagBox::CLEAR;

                for (int bidx=0; bidx < gd.refine_boxes.size(); ++bidx) {

                    const auto& rbx = gd.refine_boxes[bidx];

                    if (rbx.contains(pos)) {
                        if (gd.refine_box_type[bidx] == RefineBoxType::OnlyRefine) {
                            set_type = fab(cell); // do nothing
                        } else if (gd.refine_box_type[bidx] == RefineBoxType::NoRefine) {
                            set_type = std::max((int)TagBox::CLEAR, set_type); // clear, but only if we haven't already tagged for refinement
                        } else if (gd.refine_box_type[bidx] == RefineBoxType::ForceRefine) {
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
