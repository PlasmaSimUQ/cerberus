#ifndef REFINE_H
#define REFINE_H

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#endif

#include "MFP_factory.H"
#include "MFP.H"

using namespace amrex;

class Refinement
{
public:
    Refinement();

    virtual void get_tags(MFP* mfp, TagBoxArray& tags) const = 0;
#ifdef AMREX_USE_EB
    static void tag_cut_cells (MFP* mfp, TagBoxArray& tags, const int idx);
#endif

    int idx;
};

#endif // REFINEMENT_H
