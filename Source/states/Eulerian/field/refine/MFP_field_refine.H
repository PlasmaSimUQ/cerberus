#ifndef FIELD_REFINE_H
#define FIELD_REFINE_H

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_TagBox.H>
#include <AMReX_Vector.H>

#ifdef AMREX_USE_EB
    #include <AMReX_EB2.H>
    #include <AMReX_EBAmrUtil.H>
    #include <AMReX_EBFArrayBox.H>
    #include <AMReX_EBMultiFabUtil.H>
#endif

#include "MFP_field_defs.H"
#include "MFP_refine.H"

class FieldGradientRefinement : public Refinement
{
  public:
    FieldGradientRefinement();
    FieldGradientRefinement(const int global_idx, const sol::table& def);

    virtual void get_tags(MFP* mfp, TagBoxArray& tags) const override;

    Real refine_criteria(Array<Real, 3> S, Real low_val = 1e-12) const;

    void tag_refinement(const Box& box,
                        const FArrayBox& src,
#ifdef AMREX_USE_EB
                        const EBCellFlagFab& flags,
#endif
                        TagBox& tags,
                        const int n,
                        const Real threshold) const;

    // the index,threshold pairs variables
    Vector<std::pair<int, Real>> cons;
    int max_level;
    Real min_value;

    static bool registered;
    static std::string tag;
};

class FieldValueRefinement : public Refinement
{
  public:
    FieldValueRefinement();
    FieldValueRefinement(const int global_idx, const sol::table& def);

    virtual void get_tags(MFP* mfp, TagBoxArray& tags) const override;

    Real refine_criteria(Real S, Real low_val = 1e-12) const;

    void tag_refinement(const Box& box,
                        const FArrayBox& src,
#ifdef AMREX_USE_EB
                        const EBCellFlagFab& flags,
#endif
                        TagBox& tags,
                        const int n,
                        const Real threshold) const;

    // the index,threshold pairs variables
    Vector<std::pair<int, Real>> cons;
    int max_level;
    Real min_value;

    static bool registered;
    static std::string tag;
};

template <typename D> std::unique_ptr<Refinement> FieldRefinementBuilder(const sol::table& def)
{
    std::unique_ptr<Refinement> refine;

    const int global_idx = def["global_idx"];

    if (def.valid()) {
        if (def["name"] == D::tag) {
            return std::unique_ptr<D>(new D(global_idx, def));
        } else {
            return nullptr;
        }
    } else {
        return nullptr;
    }
}

ClassFactory<Refinement>& GetFieldRefinementFactory();

#endif  // HYDRO_REFINE_H
