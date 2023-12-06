#include "MFP_reconstruction.H"

#include <MFP_utility.H>
#include <algorithm>
#include <math.h>

Reconstruction::Reconstruction()
{
    // do nothing
}

Reconstruction::~Reconstruction()
{
    // do nothing
}

Real Reconstruction::get_slope(Vector<Real>& stencil) const
{
    BL_PROFILE("Reconstruction::get_slope");

    return 0;
}

void Reconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const
{
    // do nothing
}

ClassFactory<Reconstruction>& GetReconstructionFactory()
{
    static ClassFactory<Reconstruction> F;
    return F;
}

//================================================================================

std::string NullReconstruction::tag = "null";
bool NullReconstruction::registered = GetReconstructionFactory().Register(
  NullReconstruction::tag, ReconstructionBuilder<NullReconstruction>);

NullReconstruction::NullReconstruction()
{
    BL_PROFILE("NullReconstruction::NullReconstruction");

    stencil_length = 0;
    num_grow = 0;
}

//================================================================================

std::string ConstantReconstruction::tag = "constant";
bool ConstantReconstruction::registered = GetReconstructionFactory().Register(
  ConstantReconstruction::tag, ReconstructionBuilder<ConstantReconstruction>);

ConstantReconstruction::ConstantReconstruction()
{
    BL_PROFILE("ConstantReconstruction::ConstantReconstruction");

    stencil_length = 1;
    num_grow = 1;
}

void ConstantReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const
{
    BL_PROFILE("ConstantReconstruction::get_face_values");

    lo = stencil[0];
    hi = stencil[0];
}

//================================================================================

std::string MinModReconstruction::tag = "minmod";
bool MinModReconstruction::registered = GetReconstructionFactory().Register(
  MinModReconstruction::tag, ReconstructionBuilder<MinModReconstruction>);

MinModReconstruction::MinModReconstruction()
{
    BL_PROFILE("MinModReconstruction::MinModReconstruction");

    stencil_length = 3;
    num_grow = 2;
}

Real MinModReconstruction::get_slope(Vector<Real>& stencil) const
{
    BL_PROFILE("MinModReconstruction::get_slope");

    Real a = stencil[1] - stencil[0];
    Real b = stencil[2] - stencil[1];

    Real g = 0.0;

    if ((a >= 0.0) && (b >= 0.0)) {
        g = std::min(a, b);
    } else if ((a <= 0.0) && (b <= 0.0)) {
        g = std::max(a, b);
    }

    return g;
}

void MinModReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const
{
    BL_PROFILE("MinModReconstruction::get_face_values");

    Real g = get_slope(stencil);

    lo = stencil[1] - 0.5 * g;
    hi = stencil[1] + 0.5 * g;
}

//================================================================================

std::string VanLeerReconstruction::tag = "vanLeer";
bool VanLeerReconstruction::registered = GetReconstructionFactory().Register(
  VanLeerReconstruction::tag, ReconstructionBuilder<VanLeerReconstruction>);

VanLeerReconstruction::VanLeerReconstruction()
{
    BL_PROFILE("VanLeerReconstruction::VanLeerReconstruction");

    stencil_length = 3;
    num_grow = 2;
}

Real VanLeerReconstruction::get_slope(Vector<Real>& stencil) const
{
    BL_PROFILE("VanLeerReconstruction::get_slope");

    Real a = stencil[1] - stencil[0];
    Real b = stencil[2] - stencil[1];

    Real sa = sgn(a);
    Real sb = sgn(b);
    Real aa = std::abs(a);
    Real bb = std::abs(b);
    Real g = (sa + sb) * ((aa * bb) / (aa + bb + 1e-20));

    return g;
}

void VanLeerReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const
{
    BL_PROFILE("VanLeerReconstruction::get_face_values");

    Real g = get_slope(stencil);

    lo = stencil[1] - 0.5 * g;
    hi = stencil[1] + 0.5 * g;
}

//================================================================================

std::string MCReconstruction::tag = "MC";
bool MCReconstruction::registered = GetReconstructionFactory().Register(
  MCReconstruction::tag, ReconstructionBuilder<MCReconstruction>);

MCReconstruction::MCReconstruction()
{
    BL_PROFILE("MCReconstruction::MCReconstruction");

    stencil_length = 3;
    num_grow = 2;
}

Real MCReconstruction::get_slope(Vector<Real>& stencil) const
{
    BL_PROFILE("MCReconstruction::get_slope");

    Real x1 = stencil[0];
    Real x2 = stencil[1];
    Real x3 = stencil[2];

    Real dqm = (x2 - x1);
    Real dqp = (x3 - x2);
    Real dqc = 0.5 * (x3 - x1);

    Real s = dqm * dqp;

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

void MCReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const
{
    BL_PROFILE("MCReconstruction::get_face_values");

    Real g = get_slope(stencil);

    lo = stencil[1] - 0.5 * g;
    hi = stencil[1] + 0.5 * g;
}

//================================================================================

std::string CentredReconstruction::tag = "centre";
bool CentredReconstruction::registered = GetReconstructionFactory().Register(
  CentredReconstruction::tag, ReconstructionBuilder<CentredReconstruction>);

CentredReconstruction::CentredReconstruction()
{
    BL_PROFILE("CentredReconstruction::CentredReconstruction");

    stencil_length = 3;
    num_grow = 2;
}

Real CentredReconstruction::get_slope(Vector<Real>& stencil) const
{
    BL_PROFILE("CentredReconstruction::get_slope");

    return 0.5 * (stencil[2] - stencil[0]);
}

void CentredReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const
{
    BL_PROFILE("CentredReconstruction::get_face_values");

    Real g = get_slope(stencil);

    lo = stencil[1] - 0.5 * g;
    hi = stencil[1] + 0.5 * g;
}

//================================================================================

std::string SixthOrderReconstruction::tag = "O6";
bool SixthOrderReconstruction::registered = GetReconstructionFactory().Register(
  SixthOrderReconstruction::tag, ReconstructionBuilder<SixthOrderReconstruction>);

SixthOrderReconstruction::SixthOrderReconstruction()
{
    BL_PROFILE("SixthOrderReconstruction::SixthOrderReconstruction");

    stencil_length = 7;
    num_grow = 4;
}

Real SixthOrderReconstruction::get_slope(Vector<Real>& stencil) const
{
    BL_PROFILE("SixthOrderReconstruction::get_slope");

    // just use an approximation based on the cell edge values
    Real lo, hi;
    get_face_values(stencil, lo, hi);
    return hi - lo;
}

void SixthOrderReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const
{
    BL_PROFILE("SixthOrderReconstruction::get_face_values");

    hi = (37.0 / 60.0) * (stencil[3] + stencil[4]) - (2.0 / 15.0) * (stencil[2] + stencil[5]) +
         (1.0 / 60.0) * (stencil[1] + stencil[6]);
    lo = (37.0 / 60.0) * (stencil[3] + stencil[2]) - (2.0 / 15.0) * (stencil[4] + stencil[1]) +
         (1.0 / 60.0) * (stencil[5] + stencil[0]);
}

//================================================================================

// M.P. Martin et al. Journal of Computational Physics 220 (2006) 270-289

std::string WENOReconstruction::tag = "WENO";
bool WENOReconstruction::registered = GetReconstructionFactory().Register(
  WENOReconstruction::tag, ReconstructionBuilder<WENOReconstruction>);

WENOReconstruction::WENOReconstruction()
{
    BL_PROFILE("WENOReconstruction::WENOReconstruction");

    stencil_length = 7;
    num_grow = 4;
}

Real WENOReconstruction::get_slope(Vector<Real>& stencil) const
{
    BL_PROFILE("WENOReconstruction::get_slope");

    // just use an approximation based on the cell edge values
    Real lo, hi;
    get_face_values(stencil, lo, hi);
    return hi - lo;
}

Real WENOReconstruction::WENOSYMBO(Vector<Real>& stencil, int upwind) const
// i     : centre of stencil
// upwind: direction of upwinding (1 = -->, -1 = <--)
{
    BL_PROFILE("WENOReconstruction::WENOSYMBO");

    size_t i = stencil_length / 2;

    constexpr size_t r = 3;
    constexpr Real C[4] = {0.094647545896, 0.428074212384, 0.408289331408, 0.068988910311};
    constexpr Real akl[4][3] = {{2 / 6.0, -7 / 6.0, 11 / 6.0},
                                {-1 / 6.0, 5 / 6.0, 2 / 6.0},
                                {2 / 6.0, 5 / 6.0, -1 / 6.0},
                                {11 / 6.0, -7 / 6.0, 2 / 6.0}};
    constexpr Real b2 = std::sqrt(13 / 12.0);
    constexpr Real dkml[4][2][3] = {{{1 / 2.0, -4 / 2.0, 3 / 2.0}, {b2, -2 * b2, b2}},
                                    {{-1 / 2.0, 0, 1 / 2.0}, {b2, -2 * b2, b2}},
                                    {{-3 / 2.0, 4 / 2.0, -1 / 2.0}, {b2, -2 * b2, b2}},
                                    {{-5 / 2.0, 8 / 2.0, -3 / 2.0}, {b2, -2 * b2, b2}}};
    const Real epsilon = 1e-10;

    // calculate smoothness indicators
    Real IS[r + 1] = {0, 0, 0, 0};
    for (size_t k = 0; k <= r; ++k) {
        for (size_t m = 0; m < r - 1; ++m) {
            Real is = 0.0;
            for (size_t l = 0; l < r; ++l) {
                is += dkml[k][m][l] * stencil[i - upwind * (r - 1 - k - l)];
            }
            IS[k] += is * is;
        }
    }

    if (upwind != 0) {
        // limit the downwind indicator
        Real max_IS = IS[0];
        for (size_t k = 1; k <= r; ++k) { max_IS = std::max(IS[k], max_IS); }
        IS[r] = max_IS;
    }

    // calculate alpha
    Real alpha[r + 1];
    for (size_t k = 0; k <= r; ++k) { alpha[k] = C[k] / (epsilon + IS[k]); }

    // sum of alpha (k < r)
    Real sum_alpha = alpha[0] + alpha[1] + alpha[2];

    // calculate omega
    Real omega[r + 1];
    for (size_t k = 0; k <= r; ++k) { omega[k] = alpha[k] / sum_alpha; }

    // calculate face value
    Real face_value = 0.0;
    Real q;
    for (size_t k = 0; k < r; ++k) {
        q = 0.0;
        for (size_t l = 0; l < r; ++l) { q += akl[k][l] * stencil[i - upwind * (r - 1 - k - l)]; }
        face_value += omega[k] * q;
    }
    return face_value;
}

void WENOReconstruction::get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const
{
    BL_PROFILE("WENOReconstruction::get_face_values");

    hi = WENOSYMBO(stencil, 1);
    lo = WENOSYMBO(stencil, -1);
}

//================================================================================
