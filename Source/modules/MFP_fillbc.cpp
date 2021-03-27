#include "MFP_fillbc.H"
#include <AMReX_filcc_f.H>

#include "MFP_global.H"
#include "MFP_state.H"

#include <sstream>

using GD = GlobalData;

void StateBC::add_function(const std::string name, const Optional3D1VFunction &fun)
{
    functions[name] = fun;

    if (fun.has_func()) {
        has_function = true;
    }
    if (fun.is_valid()) {
        is_valid = true;
    }
}

void BoundaryState::set(int d, const std::string name, int s, const Optional3D1VFunction &f)
{
    data[s][d].add_function(name, f);
}

void BoundaryState::set_lo(int d, const std::string name, const Optional3D1VFunction &f)
{
    data[0][d].add_function(name, f);
}

void BoundaryState::set_hi(int d, const std::string name, const Optional3D1VFunction &f)
{
    data[1][d].add_function(name, f);
}

const Optional3D1VFunction& BoundaryState::get(int hl, int dir, const std::string name) const {
    return data[hl][dir].functions.at(name);
}

const Optional3D1VFunction& BoundaryState::get_lo(int dir, const std::string name) const {
    return data[0][dir].functions.at(name);
}

Optional3D1VFunction& BoundaryState::get_lo(int dir, const std::string name) {
    return data[0][dir].functions.at(name);
}


const Optional3D1VFunction& BoundaryState::get_hi(int dir, const std::string name) const {
    return data[1][dir].functions.at(name);
}

Optional3D1VFunction& BoundaryState::get_hi(int dir, const std::string name) {
    return data[1][dir].functions[name];
}

const StateBC& BoundaryState::get_lo(int dir) const
{
    return data[0][dir];
}

const StateBC& BoundaryState::get_hi(int dir) const
{
    return data[1][dir];
}

void BoundaryState::post_init()
{
    bool has_inflow;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // lo side
        has_inflow = false;
        if (data[0][d].is_valid) {
            where_is_inflow.setLo(d, BCType::ext_dir);
            has_valid_ext_dir = true;
        }

        // hi side
        has_inflow = false;
        if (data[1][d].is_valid) {
            where_is_inflow.setHi(d, BCType::ext_dir);
            has_valid_ext_dir = true;
        }
    }
}

std::string BoundaryState::str(const std::string &prefix) const {
    std::stringstream ss;

    const Array<std::string, 3> ax = {"x","y","z"};
    const Array<std::string, 2> hl = {"lo", "hi"};

    // see AMReX_BC_TYPES.H for definitions
    const std::map<int,std::string> bc_types = {
        {-666,"bogus"},
        {-1,"reflect_odd"},
        {0, "int_dir"},
        {1, "reflect_even"},
        {2, "foextrap"},
        {3, "ext_dir"},
        {4, "hoextrap"},
        {5, "hoextrapcc"},
    };

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        for (int i=0; i<2; ++i) {

            Vector<std::string> bc;
            for (const auto &bcr : fill_bc) {
                bc.push_back(bc_types.at(bcr.data()[i*AMREX_SPACEDIM + d]));
            }

            ss << prefix << ax[d] << "(" << hl[i] << "): \n"
               << prefix << "  " << vec2str(bc,"type=[");
            if (where_is_inflow.data()[i*AMREX_SPACEDIM + d] == BCType::ext_dir) {

                // get all of the keys
                const std::map<std::string, Optional3D1VFunction> &m = data[i][d].functions;
                Vector<std::string> keys;
                for(std::map<std::string,Optional3D1VFunction>::const_iterator it = m.begin(); it != m.end(); ++it) {
                  keys.push_back(it->first);
                }


                ss << "\n" << prefix << "  " << vec2str(keys, "UDF =[");
            }
            ss << std::endl;
        }
    }
    return ss.str();
}

FillBC* FillBC::clone () const
{
    return new FillBC(*this);
}

bool FillBC::hasFabVersion () const noexcept
{
    return true;
}

void FillBC::operator() (Box const& bx, FArrayBox& dest,
                         const int dcomp, const int numcomp,
                         Geometry const& geom, const Real time,
                         const Vector<BCRec>& bcr, const int bcomp,
                         const int orig_comp) const
{
    const int* lo = dest.loVect();
    const Box& domain = geom.Domain();
    const int* dom_lo = domain.loVect();
    const Real* dx = geom.CellSize();
    const Real* problo = geom.ProbLo();
    //    const int* bc = bcr[bcomp].vect();
    Real xlo[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        xlo[i] = problo[i] + dx[i]*(lo[i]-dom_lo[i]);
    }

    fab_filcc(bx, dest.array(dcomp), numcomp, domain, dx, xlo, &bcr[0]);
    //    amrex_fab_filcc(BL_TO_FORTRAN_N_ANYD(dest,dcomp), &numcomp,
    //                      BL_TO_FORTRAN_BOX(domain),
    //                      dx, xlo, bc);

}

//==================================================================================

NullFillBC* NullFillBC::clone () const
{
    return new NullFillBC(*this);
}

bool NullFillBC::hasFabVersion () const noexcept
{
    return true;
}

void NullFillBC::operator() (Box const& bx, FArrayBox& dest,
                             const int dcomp, const int numcomp,
                             Geometry const& geom, const Real time,
                             const Vector<BCRec>& bcr, const int bcomp,
                             const int orig_comp) const
{
    amrex::Abort("How did we get here?");
}
