#ifndef MFP_FILLBC_H
#define MFP_FILLBC_H

#include <AMReX_StateDescriptor.H>
#include <MFP_utility.H>
#include "MFP_optional_func.H"

using namespace amrex;

struct StateBC
{
    std::map<std::string,Optional3D1VFunction> functions;
    bool has_function = false;
    bool is_valid = false;
    void add_function(const std::string name, const Optional3D1VFunction &fun);
};

class BoundaryState
{
public:
    void post_init();
    std::string str(const std::string &prefix="") const;

    void set_lo(int d, const std::string name, const Optional3D1VFunction &f);
    void set_hi(int d, const std::string name, const Optional3D1VFunction &f);
    void set(int d, const std::string name, int s, const Optional3D1VFunction &f);

    const Optional3D1VFunction& get(int hl, int dir, const std::string name) const;

    const Optional3D1VFunction& get_lo(int dir, const std::string name) const;
    Optional3D1VFunction& get_lo(int dir, const std::string name);
    const Optional3D1VFunction& get_hi(int dir, const std::string name) const;
    Optional3D1VFunction& get_hi(int dir, const std::string name);

    const StateBC& get_lo(int dir) const;
    const StateBC& get_hi(int dir) const;

    Vector<BCRec> phys_fill_bc;
    Vector<BCRec> fill_bc; // what to run before performing an inflow
    BCRec where_is_inflow;
#ifdef AMREX_USE_EB
    BCRec eb_bc;
#endif
    bool has_valid_ext_dir = false;
private:

    // [dimension, state]
    Array<Array<StateBC, AMREX_SPACEDIM>, 2> data; // the function to run for an inflow condition

};

class FillBC : public StateDescriptor::BndryFunc
{
public:

    FillBC() noexcept {}
    FillBC(BoundaryState& inflow, const int idx);

    virtual FillBC* clone () const override;

    virtual bool hasFabVersion () const noexcept override;

    virtual void operator() (Box const& bx, FArrayBox& dest,
                             const int dcomp, const int numcomp,
                             Geometry const& geom, const Real time,
                             const Vector<BCRec>& bcr, const int bcomp,
                             const int orig_comp) const override;

};

class NullFillBC : public StateDescriptor::BndryFunc
{
public:

    NullFillBC() noexcept {}

    virtual NullFillBC* clone () const override;

    virtual bool hasFabVersion () const noexcept override;

    virtual void operator() (Box const& bx, FArrayBox& dest,
                             const int dcomp, const int numcomp,
                             Geometry const& geom, const Real time,
                             const Vector<BCRec>& bcr, const int bcomp,
                             const int orig_comp) const override;

};



#endif // MFP_FILLBC_H
