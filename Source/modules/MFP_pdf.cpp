
#include "MFP_pdf.h"

PDF::PDF()
{

}

PDF::~PDF()
{

}

ClassFactory<PDF>& GetPDFFactory()
{
    static ClassFactory<PDF> F;
    return F;
}

//-----------------------------------------------------------------------------

std::string Maxwellian::tag = "Maxwell";
bool Maxwellian::registered = GetPDFFactory().Register(Maxwellian::tag, PDFBuilder<Maxwellian>);

Maxwellian::Maxwellian(){}

Maxwellian::Maxwellian(const sol::table& def)
{
    temperature = get_udf(def["T"]);
    x_vel = get_udf(def["x_vel"]);
    y_vel = get_udf(def["y_vel"]);
    z_vel = get_udf(def["z_vel"]);

    std::random_device rd;
    mt = std::unique_ptr<std::mt19937> (new std::mt19937(rd()));


}

Array<Real,3> Maxwellian::operator()(Real x, Real y, Real z, Real t) const
{

    Real vx = x_vel(x,y,z,t);
    Real vy = y_vel(x,y,z,t);
    Real vz = z_vel(x,y,z,t);
    Real T = temperature(x,y,z,t);

    std::normal_distribution<Real> fu(vx,T);
    std::normal_distribution<Real> fv(vy,T);
    std::normal_distribution<Real> fw(vz,T);

    return {fu(*mt), fv(*mt), fw(*mt)};
}
