#ifndef MFP_RECONSTRUCTION_H
#define MFP_RECONSTRUCTION_H

#include <AMReX_Vector.H>
#include <AMReX_REAL.H>

#include "MFP_factory.H"

using namespace amrex;

class Reconstruction
{
public:
  Reconstruction();
  virtual ~Reconstruction();

  virtual Real get_slope(Vector<Real>& stencil);
  virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi);
  virtual int get_num_grow()=0;

  virtual const std::string& get_tag() const = 0;

  int stencil_length;
  int num_grow;
};

template <typename D>
std::unique_ptr<Reconstruction> ReconstructionBuilder(const sol::table& def)
{
    if (def["reconstruction"] == D::tag) {
        return std::unique_ptr<D>(new D());
    } else {
        return nullptr;
    }
}

PhysicsFactory<Reconstruction>& GetReconstructionFactory();

//================================================================================

class NullReconstruction : public Reconstruction
{
public:
  NullReconstruction();
  virtual const std::string& get_tag() const override {return tag;}
  virtual int get_num_grow() override {return num_grow;}
  static std::string tag;
  static bool registered;
};

//================================================================================

class ConstantReconstruction : public Reconstruction
{
public:
  ConstantReconstruction();
  virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) override;
  virtual int get_num_grow() override {return num_grow;}

  virtual const std::string& get_tag() const override {return tag;}
  static std::string tag;
  static bool registered;
};

//================================================================================

class MinModReconstruction : public Reconstruction
{
public:
  MinModReconstruction();
  virtual Real get_slope(Vector<Real>& stencil) override;
  virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) override;
  virtual int get_num_grow() override {return num_grow;}

  virtual const std::string& get_tag() const override {return tag;}
  static std::string tag;
  static bool registered;
};

//================================================================================

class VanLeerReconstruction : public Reconstruction
{
public:
  VanLeerReconstruction();
  virtual Real get_slope(Vector<Real>& stencil) override;
  virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) override;
  virtual int get_num_grow() override {return num_grow;}

  virtual const std::string& get_tag() const override {return tag;}
  static std::string tag;
  static bool registered;
};

//================================================================================

class MCReconstruction : public Reconstruction
{
public:
  MCReconstruction();
  virtual Real get_slope(Vector<Real>& stencil) override;
  virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) override;
  virtual int get_num_grow() override {return num_grow;}

  virtual const std::string& get_tag() const override {return tag;}
  static std::string tag;
  static bool registered;
};

//================================================================================

class CentredReconstruction : public Reconstruction
{
public:
  CentredReconstruction();
  virtual Real get_slope(Vector<Real>& stencil) override;
  virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) override;
  virtual int get_num_grow() override {return num_grow;}

  virtual const std::string& get_tag() const override {return tag;}
  static std::string tag;
  static bool registered;
};

//================================================================================

class SixthOrderReconstruction : public Reconstruction
{
public:
  SixthOrderReconstruction();
  virtual Real get_slope(Vector<Real>& stencil) override;
  virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) override;
  virtual int get_num_grow() override {return num_grow;}

  virtual const std::string& get_tag() const override {return tag;}
  static std::string tag;
  static bool registered;
};

//================================================================================

class WENOReconstruction : public Reconstruction
{
public:
  WENOReconstruction();
  Real WENOSYMBO(Vector<Real>& stencil, int upwind);
  virtual Real get_slope(Vector<Real>& stencil) override;
  virtual void get_face_values(Vector<Real>& stencil, Real& lo, Real& hi) override;
  virtual int get_num_grow() override {return num_grow;}

  virtual const std::string& get_tag() const override {return tag;}
  static std::string tag;
  static bool registered;
};

//================================================================================

#endif // MFP_RECONSTRUCTION_H
