#ifndef MFP_RECONSTRUCTION_H
#define MFP_RECONSTRUCTION_H

#include "MFP_factory.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

using namespace amrex;

class Reconstruction
{
  public:
    Reconstruction();
    virtual ~Reconstruction();

    virtual Real get_slope(Vector<Real>& stencil) const;
    virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const;
    virtual int get_num_grow() const { return 0; }

    virtual std::string get_tag() const { return ""; }

    int stencil_length;
    int num_grow;
};

template <typename D> std::unique_ptr<Reconstruction> ReconstructionBuilder(const sol::table& def)
{
    if (def["reconstruction"] == D::tag) {
        return std::unique_ptr<D>(new D());
    } else {
        return nullptr;
    }
}

ClassFactory<Reconstruction>& GetReconstructionFactory();

//================================================================================

class NullReconstruction : public Reconstruction
{
  public:
    NullReconstruction();
    virtual std::string get_tag() const override { return tag; }
    virtual int get_num_grow() const override { return num_grow; }
    static std::string tag;
    static bool registered;
};

//================================================================================

class ConstantReconstruction : public Reconstruction
{
  public:
    ConstantReconstruction();
    virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const override;
    virtual int get_num_grow() const override { return num_grow; }

    virtual std::string get_tag() const override { return tag; }
    static std::string tag;
    static bool registered;
};

//================================================================================

class MinModReconstruction : public Reconstruction
{
  public:
    MinModReconstruction();
    virtual Real get_slope(Vector<Real>& stencil) const override;
    virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const override;
    virtual int get_num_grow() const override { return num_grow; }

    virtual std::string get_tag() const override { return tag; }
    static std::string tag;
    static bool registered;
};

//================================================================================

class VanLeerReconstruction : public Reconstruction
{
  public:
    VanLeerReconstruction();
    virtual Real get_slope(Vector<Real>& stencil) const override;
    virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const override;
    virtual int get_num_grow() const override { return num_grow; }

    virtual std::string get_tag() const override { return tag; }
    static std::string tag;
    static bool registered;
};

//================================================================================

class MCReconstruction : public Reconstruction
{
  public:
    MCReconstruction();
    virtual Real get_slope(Vector<Real>& stencil) const override;
    virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const override;
    virtual int get_num_grow() const override { return num_grow; }

    virtual std::string get_tag() const override { return tag; }
    static std::string tag;
    static bool registered;
};

//================================================================================

class CentredReconstruction : public Reconstruction
{
  public:
    CentredReconstruction();
    virtual Real get_slope(Vector<Real>& stencil) const override;
    virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const override;
    virtual int get_num_grow() const override { return num_grow; }

    virtual std::string get_tag() const override { return tag; }
    static std::string tag;
    static bool registered;
};

//================================================================================

class SixthOrderReconstruction : public Reconstruction
{
  public:
    SixthOrderReconstruction();
    virtual Real get_slope(Vector<Real>& stencil) const override;
    virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const override;
    virtual int get_num_grow() const override { return num_grow; }

    virtual std::string get_tag() const override { return tag; }
    static std::string tag;
    static bool registered;
};

//================================================================================

class WENOReconstruction : public Reconstruction
{
  public:
    WENOReconstruction();
    Real WENOSYMBO(Vector<Real>& stencil, int upwind) const;
    virtual Real get_slope(Vector<Real>& stencil) const override;
    virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) const override;
    virtual int get_num_grow() const override { return num_grow; }

    virtual std::string get_tag() const override { return tag; }
    static std::string tag;
    static bool registered;
};

//================================================================================

#endif  // MFP_RECONSTRUCTION_H
