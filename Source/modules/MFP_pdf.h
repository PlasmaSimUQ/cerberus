#ifndef PDF_H
#define PDF_H

#include <AMReX_REAL.H>
#include <string>
#include <random>

#include "MFP_factory.H"
#include "MFP_optional_func.H"

using namespace amrex;

class PDF
{
public:
    PDF();
    virtual ~PDF();

    virtual const std::string& get_tag() const = 0;
    virtual Array<Real,3> operator()(Real x=0, Real y=0, Real z=0, Real t=0) const = 0;
};

template <typename D>
std::unique_ptr<PDF> PDFBuilder(const sol::table& def)
{

    if (def["distribution"] == D::tag) {
        return std::unique_ptr<D>(new D(def));
    } else {
        return nullptr;
    }
}

ClassFactory<PDF>& GetPDFFactory();

//-----------------------------------------------------------------------------

class Maxwellian : public PDF
{
public:
    Maxwellian();
    Maxwellian(const sol::table& def);

    virtual const std::string& get_tag() const override {return tag;}
    virtual Array<Real,3> operator()(Real x=0, Real y=0, Real z=0, Real t=0) const override;

    static std::string tag;
    static bool registered;

    std::unique_ptr<std::mt19937> mt;

private:
    Optional3D1VFunction temperature;
    Optional3D1VFunction x_vel, y_vel, z_vel;
};

#endif // PDF_H
