#ifndef MFP_EB_DIVERGENCE_H
#define MFP_EB_DIVERGENCE_H

#ifdef AMREX_USE_EB

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_IArrayBox.H>

#include "MFP_factory.H"

#include "sol.hpp"


using namespace amrex;

class DivergenceEB
{
public:
    DivergenceEB();
    ~DivergenceEB();

    virtual void calc_eb_divergence(const Box& box,
                                    const FArrayBox &cons,
                                    Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                    Array<FArrayBox, AMREX_SPACEDIM> &wall_fluxes,
                                    FArrayBox& du,
                                    const EBCellFlagFab& flag,
                                    const FArrayBox& vfrac,
                                    const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                                    const Array<const FArrayBox *, AMREX_SPACEDIM> &fcent,
                                    int as_crse,
                                    int as_fine,
                                    const IArrayBox* rrflag_as_crse,
                                    const IArrayBox& levmsk,
                                    FArrayBox* rr_drho_crse,
                                    FArrayBox &dm_as_fine,
                                    const Real *dx,
                                    const Real dt) const;

    virtual void merge_cells(const Box& box,
                            FArrayBox &cons,
                            FArrayBox& du,
                            const EBCellFlagFab& flag,
                            const FArrayBox& vfrac,
                            const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                            int as_fine,
                            FArrayBox &dm_as_fine,
                            const IArrayBox& levmsk) const;

    static bool is_inside(const int i,const int j, const int k, const Dim3 &lo, const Dim3 &hi);

    virtual const std::string& get_tag() const = 0;
    virtual std::string str() const = 0;

    int global_idx;
};

template <typename D>
std::unique_ptr<DivergenceEB> DivergenceEBBuilder(const sol::table& def)
{

    sol::optional<sol::table> div_def = def["eb_divergence"];

    if (!div_def) {
        return nullptr;
    }

    std::string type = div_def->get<std::string>("type");
    if (type == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

PhysicsFactory<DivergenceEB>& GetDivergenceEBBuilder();

//=============================================================================

class RedistributeEB : public DivergenceEB
{
public:
    RedistributeEB();
    RedistributeEB(const sol::table& def);
    ~RedistributeEB();

    enum class RedistributionEB : int {
            Uniform=0,
            VolumeFraction,
            Density,
            Energy
        };

    virtual void calc_eb_divergence(const Box& box,
                                    const FArrayBox &cons,
                                    Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                    Array<FArrayBox, AMREX_SPACEDIM> &wall_fluxes,
                                    FArrayBox& du,
                                    const EBCellFlagFab& flag,
                                    const FArrayBox& vfrac,
                                    const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                                    const Array<const FArrayBox *, AMREX_SPACEDIM> &fcent,
                                    int as_crse,
                                    int as_fine,
                                    const IArrayBox* rrflag_as_crse,
                                    const IArrayBox& levmsk,
                                    FArrayBox* rr_drho_crse,
                                    FArrayBox &dm_as_fine,
                                    const Real *dx,
                                    const Real dt) const override;

    virtual const std::string& get_tag() const override {return tag;}
    virtual std::string str() const override;

    static std::string tag;
    static bool registered;
    static Vector<std::string> options;

    RedistributionEB redistribution_strategy;
    Real reredistribution_threshold;
};

//=============================================================================

class MergeEB : public DivergenceEB
{
public:
    MergeEB();
    MergeEB(const sol::table& def);
    ~MergeEB();

    virtual void calc_eb_divergence(const Box& box,
                                    const FArrayBox &cons,
                                    Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                    Array<FArrayBox, AMREX_SPACEDIM> &wall_fluxes,
                                    FArrayBox& du,
                                    const EBCellFlagFab& flag,
                                    const FArrayBox& vfrac,
                                    const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                                    const Array<const FArrayBox *, AMREX_SPACEDIM> &fcent,
                                    int as_crse,
                                    int as_fine,
                                    const IArrayBox* rrflag_as_crse,
                                    const IArrayBox& levmsk,
                                    FArrayBox* rr_drho_crse,
                                    FArrayBox &dm_as_fine,
                                    const Real *dx,
                                    const Real dt) const override;

    virtual void merge_cells(const Box& box,
                            FArrayBox &cons,
                            FArrayBox& du,
                            const EBCellFlagFab& flag,
                            const FArrayBox& vfrac,
                            const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                            int as_fine,
                            FArrayBox &dm_as_fine,
                            const IArrayBox& levmsk) const override;

    virtual const std::string& get_tag() const override {return tag;}
    virtual std::string str() const override;
    static std::string tag;
    static bool registered;

    Real merge_threshold = 0.5;
};

#endif // AMREX_USE_EB
#endif // MFP_EB_DIVERGENCE_H
