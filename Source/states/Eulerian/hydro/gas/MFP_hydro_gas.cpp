#include "MFP_hydro_gas.H"

#include "MFP.H"
#include "MFP_state.H"
#include "MFP_hydro.H"

// ====================================================================================

HydroGas::HydroGas()
{
}

HydroGas::~HydroGas()
{
}

Real HydroGas::get_density_from_number_density(const Real nd, Vector<Real>& Q) const
{

    Real mass_inv = 0.0;
    if (mass_const)
        return nd * mass[0];

    Real alphai = 0.0;
    Real S_alphas = 0.0;
    for (int i=0; i < n_tracers(); ++i) {
        alphai = Q[+HydroDef::PrimIdx::NUM + i];
        alphai = std::clamp(alphai, 0.0, 1.0);
        S_alphas += alphai;
        mass_inv += alphai / mass[i];
    }
    mass_inv += (1.0-S_alphas) / mass[n_tracers()]; // deduce alpha value for final component

    return nd / mass_inv;
}

void HydroGas::write_info(nlohmann::json &js) const
{
    js["mass"] = mass;
    js["charge"] = charge;
    js["comp_names"] = comp_names;
}

ClassFactory<HydroGas>& GetHydroGasFactory() {
    static ClassFactory<HydroGas> F;
    return F;
}
