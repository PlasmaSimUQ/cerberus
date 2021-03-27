#include "MFP_reconstruction.H"

#include <algorithm>
#include <math.h>
#include <MFP_utility.H>

Reconstruction::Reconstruction()
{
    // do nothing
}

Reconstruction::~Reconstruction()
{
    // do nothing
}

Real Reconstruction::get_slope(Vector<Real>& stencil)
{
    return 0;
}

Real Reconstruction::get_slope(Vector<Real>& stencil, const int side)
{
    if (side == 0) {
        return get_slope(stencil);
    } else if (side == 1) {
        return stencil[stencil_length/2+1] - stencil[stencil_length/2];
    } else {
        return stencil[stencil_length/2] - stencil[stencil_length/2-1];
    }
}

void Reconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi)
{
    // do nothing
}

void Reconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi, const int side)
{
    if (side == 0) {
        get_face_values(stencil, lo, hi);
    } else {

        Real g = get_slope(stencil, side);

        lo = stencil[stencil_length/2] - 0.5*g;
        hi = stencil[stencil_length/2] + 0.5*g;
    }
}

PhysicsFactory<Reconstruction>& GetReconstructionFactory()
{
    static PhysicsFactory<Reconstruction> F;
    return F;
}

//================================================================================

std::string NullReconstruction::tag = "null";
bool NullReconstruction::registered = GetReconstructionFactory().Register(NullReconstruction::tag, ReconstructionBuilder<NullReconstruction>);

NullReconstruction::NullReconstruction()
{
    stencil_length = 0;
    num_grow = 0;
}

//================================================================================

std::string ConstantReconstruction::tag = "constant";
bool ConstantReconstruction::registered = GetReconstructionFactory().Register(ConstantReconstruction::tag, ReconstructionBuilder<ConstantReconstruction>);

ConstantReconstruction::ConstantReconstruction()
{
    stencil_length = 1;
    num_grow = 1; 
}

void ConstantReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi)
{
    lo = stencil[0];
    hi = stencil[0];
}


//================================================================================

std::string MinModReconstruction::tag = "minmod";
bool MinModReconstruction::registered = GetReconstructionFactory().Register(MinModReconstruction::tag, ReconstructionBuilder<MinModReconstruction>);

MinModReconstruction::MinModReconstruction()
{
    stencil_length = 3;
    num_grow = 2;
}

Real MinModReconstruction::get_slope(Vector<Real>& stencil)
{
    Real a = stencil[1] - stencil[0];
    Real b = stencil[2] - stencil[1];

    Real g = 0.0;

    if ((a >= 0.0) && (b >= 0.0)) {
        g = std::min(a, b);
    } else if ((a <= 0.0) && (b <= 0.0)) {
        g = std::max(a,b);
    }

    return g;
}

void MinModReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi)
{

    Real g = get_slope(stencil);

    lo = stencil[1] - 0.5*g;
    hi = stencil[1] + 0.5*g;
}

//================================================================================

std::string VanLeerReconstruction::tag = "vanLeer";
bool VanLeerReconstruction::registered = GetReconstructionFactory().Register(VanLeerReconstruction::tag, ReconstructionBuilder<VanLeerReconstruction>);

VanLeerReconstruction::VanLeerReconstruction()
{
    stencil_length = 3;
    num_grow = 2;
}

Real VanLeerReconstruction::get_slope(Vector<Real>& stencil)
{
    Real a = stencil[1] - stencil[0];
    Real b = stencil[2] - stencil[1];

    Real sa = sgn(a);
    Real sb = sgn(b);
    Real aa = std::abs(a);
    Real bb = std::abs(b);
    Real g = (sa + sb)*((aa*bb)/(aa + bb + 1e-20));

    return g;
}

void VanLeerReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi)
{

    Real g = get_slope(stencil);

    lo = stencil[1] - 0.5*g;
    hi = stencil[1] + 0.5*g;
}

//================================================================================

std::string MCReconstruction::tag = "MC";
bool MCReconstruction::registered = GetReconstructionFactory().Register(MCReconstruction::tag, ReconstructionBuilder<MCReconstruction>);

MCReconstruction::MCReconstruction()
{
    stencil_length = 3;
    num_grow = 2;
}

Real MCReconstruction::get_slope(Vector<Real>& stencil)
{
    Real x1 = stencil[0];
    Real x2 = stencil[1];
    Real x3 = stencil[2];

    Real dqm = (x2 - x1);
    Real dqp = (x3 - x2);
    Real dqc = 0.5*(x3 - x1);

    Real s = dqm*dqp;

    if (s <= 0.0) {
        s = 0.0;
    } else {
        if ((std::abs(dqm) < std::abs(dqp)) && (std::abs(dqm) < std::abs(dqc))) {
            s = dqm;
        } else if (std::abs(dqp) < std::abs(dqc)) {
            s = dqp;
        } else {
            s = dqc;
        }
    }

    return s;
}

void MCReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi)
{

    Real g = get_slope(stencil);

    lo = stencil[1] - 0.5*g;
    hi = stencil[1] + 0.5*g;
}

//================================================================================

std::string CentredReconstruction::tag = "centre";
bool CentredReconstruction::registered = GetReconstructionFactory().Register(CentredReconstruction::tag, ReconstructionBuilder<CentredReconstruction>);

CentredReconstruction::CentredReconstruction()
{
    stencil_length = 3;
    num_grow = 2;
}

Real CentredReconstruction::get_slope(Vector<Real>& stencil)
{
    return 0.5*(stencil[2] - stencil[0]);
}

void CentredReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi)
{
    Real g =get_slope(stencil);

    lo = stencil[1] - 0.5*g;
    hi = stencil[1] + 0.5*g;
}

//================================================================================

std::string SixthOrderReconstruction::tag = "O6";
bool SixthOrderReconstruction::registered = GetReconstructionFactory().Register(SixthOrderReconstruction::tag, ReconstructionBuilder<SixthOrderReconstruction>);

SixthOrderReconstruction::SixthOrderReconstruction()
{
    stencil_length = 7;
    num_grow = 4;
}

Real SixthOrderReconstruction::get_slope(Vector<Real>& stencil)
{
    // just use an approximation based on the cell edge values
    Real lo, hi;
    get_face_values(stencil, lo, hi);
    return hi - lo;
}

void SixthOrderReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi)
{
    hi = (37.0/60.0)*(stencil[3] + stencil[4]) - (2.0/15.0)*(stencil[2] + stencil[5]) + (1.0/60.0)*(stencil[1] + stencil[6]);
    lo = (37.0/60.0)*(stencil[3] + stencil[2]) - (2.0/15.0)*(stencil[4] + stencil[1]) + (1.0/60.0)*(stencil[5] + stencil[0]);
}

//================================================================================
