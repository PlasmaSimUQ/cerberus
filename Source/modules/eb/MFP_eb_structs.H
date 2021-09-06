#ifndef MFP_EB_STRUCTS_H
#define MFP_EB_STRUCTS_H

#ifdef AMREX_USE_EB

#include <AMReX_EB2.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>

using namespace amrex;

struct EBData {
    std::unique_ptr<EBFArrayBoxFactory> ebfactory;
    FabArray<EBCellFlagFab> flags;
    MultiFab volfrac;
    const amrex::MultiCutFab* bndrycent;
    const amrex::MultiCutFab* bndrynorm;
    Array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
    Array<const MultiCutFab*, AMREX_SPACEDIM> facecent;
    MultiCutFab bndryidx;
};

struct DataEB {
    int index;
    const EB2::IndexSpace* index_space;
    std::function<Real (AMREX_D_DECL(Real, Real, Real))> geom_func;
    int insertion_type; // add = 1, subtract = -1
    bool inside;
    Vector<std::pair<int,int>> states; // [state global index, bc index withn state]
    bool regular;
};
#endif

#endif // MFP_EB_STRUCTS_H
