#ifndef AMREX_USE_EB
#include "MFP_CTU_Braginskii.H"
#include "MFP.H"
#include "MFP_state.H"
#include "MFP_hydro.H"
#include "MFP_field.H"
#include "sol.hpp"
#include "Eigen"
#include "Dense"
#include "MFP_diagnostics.H"
#include "MFP_rk4.H"

std::string BraginskiiCTU::tag = "BraginskiiCTU";
bool BraginskiiCTU::registered = GetActionFactory().Register(BraginskiiCTU::tag, ActionBuilder<BraginskiiCTU>);

bool BraginskiiCTU::braginskii_anisotropic = true;
bool BraginskiiCTU::srin_switch = false;

BraginskiiCTU::BraginskiiCTU(){}
BraginskiiCTU::~BraginskiiCTU(){}

BraginskiiCTU::BraginskiiCTU(const int idx, const sol::table &def)
{

    Abort("You have selected an action of type 'BraginskiiCTU' which is currently incomplete, comment out this Abort() only if you know what you are doing.");

    action_idx = idx;
    name = def["name"];

    do_CTU = def.get_or("corner_transport", true);

    const sol::table state_names = def["states"];

    ion_state = &HydroState::get_state(state_names["ion"]);
    states[+BraginskiiStateIdx::Ion] = ion_state;
    state_indexes.push_back(ion_state->global_idx);

    electron_state = &HydroState::get_state(state_names["electron"]);
    states[+BraginskiiStateIdx::Electron] = electron_state;
    state_indexes.push_back(electron_state->global_idx);

    field_state = &FieldState::get_state(state_names["field"]);
    states[+BraginskiiStateIdx::Field] = field_state;
    state_indexes.push_back(field_state->global_idx);

    srin_switch = def.get_or("srin_switch",false);
    braginskii_anisotropic = def.get_or("anisotropic",false);
    cfl = def.get_or("cfl",1.0);

    ion_coeffs.forceViscosityValue = def.get_or("force_ion_viscosity",0.0);
    ion_coeffs.forceViscosity = ion_coeffs.forceViscosityValue > 0.0;

    electron_coeffs.forceViscosityValue = def.get_or("force_electron_viscosity",0.0);
    electron_coeffs.forceViscosity = electron_coeffs.forceViscosityValue > 0.0;

    time_refinement_factor = def.get_or("time_refinement_factor",10);
    max_time_refinement = def.get_or("max_time_refinement_levels",10);


    return;
}

void BraginskiiCTU::get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const
{
    BL_PROFILE("BraginskiiCTU::get_data");

    Vector<Array<int,2>> options = {
        {ion_state->global_idx, 1},{electron_state->global_idx, 1}, {field_state->global_idx, 1}
    };

    Action::get_data(mfp, options, update, time);

}

Real BraginskiiCTU::get_coulomb_logarithm(const Real& T_i, const Real& T_e, const Real& nd_e)
{
    BL_PROFILE("BraginskiiSource::get_coulomb_logarithm");
    // Coulomb logairthm as reported in BRaginskii OG paper.
    // Alternative is to use the formulation from
    // "Ionic transport in high-energy-density matter" Stanton 2016
    return 10.;
    Real T_ref = MFP::T_ref;
    Real n_ref = MFP::n_ref;

    if (T_e < 50*11600/T_ref) {//where refernce value is in K and conversion 1 eV = 11600K
        return 23.4 - 1.15 * log10( nd_e*n_ref ) + 3.45 * log10( T_e*T_ref/11600 );

    }
    else {
        Real val = 25.3 - 1.15 * log10( nd_e*n_ref ) + 2.3*log10( T_e*T_ref/11600 );
        return val;
        //return 25.3 - 1.15 * log10( nd_e*n_ref ) + 2.3*log10( T_e*T_ref/11600 );
    }
}


void BraginskiiCTU::get_ion_coeffs(const Vector<Real>& Q_i,
                                   const Vector<Real>& Q_e,
                                   const Array<Real,3>& B_xyz,
                                   Real& T_i,
                                   Real& eta0,
                                   Real& eta1,
                                   Real& eta2,
                                   Real& eta3,
                                   Real& eta4,
                                   Real& kappa1,
                                   Real& kappa2,
                                   Real& kappa3,
                                   int& truncatedTau)
{
    BL_PROFILE("BraginskiiCTU::get_ion_coeffs");
    truncatedTau = 0;
    Real mass_i,mass_e,charge_i,charge_e,T_e,nd_i,nd_e;
    HydroState &istate = *ion_state;
    HydroState &estate = *electron_state;

    //Extract and assign parameters from Q_i and Q_e
    //--- electron state and properties required for calcs -------Note move this
    charge_e= estate.gas->get_charge_from_prim(Q_e); // electron propertis
    mass_e  = estate.gas->get_mass_from_prim(Q_e);
    T_e     = estate.gas->get_temperature_from_prim(Q_e);
    nd_e    = Q_e[+HydroDef::PrimIdx::Density]/mass_e;
    //--- ion state and properties required for calcs
    charge_i= istate.gas->get_charge_from_prim(Q_i);
    mass_i  = istate.gas->get_mass_from_prim(Q_i);
    T_i     = istate.gas->get_temperature_from_prim(Q_i); //
    nd_i    = Q_i[+HydroDef::PrimIdx::Density]/mass_i;
    //Magnetic field
    Real Bx=B_xyz[0], By=B_xyz[1], Bz=B_xyz[2];

    // See page 215 (document numbering) of Braginskii's original transport paper
    Real t_collision_ion, p_lambda, omega_ci, omega_p;
    p_lambda = get_coulomb_logarithm(T_i,T_e,nd_e);
    // absence of the boltzmann constant due to usage of nondimensional variables.
    // Note that here charge_i = e^4*Z^4

    Real Debye = MFP::Debye, Larmor = MFP::Larmor;

    Real n0_ref=MFP::n0;


    t_collision_ion = std::pow(Debye,4)*n0_ref
            *(12*std::sqrt(mass_i)*std::pow(3.14159*T_i, 3./2.)) /
            (p_lambda * std::pow(charge_i,4) * nd_i);

    omega_ci = charge_i * std::sqrt( Bx*Bx + By*By + Bz*Bz ) / mass_i / Larmor;
    omega_p  = std::sqrt(nd_i*charge_i*charge_i/mass_i/Debye/Debye) ;

    if (1/t_collision_ion < effective_zero) {

        t_collision_ion = 1/effective_zero;
    }

    if (srin_switch && (1/t_collision_ion < omega_ci/10/2/3.14159) && (1/t_collision_ion < omega_p/10/2/3.14159)) {
        t_collision_ion = 1/std::min(omega_ci/2/3.14159, omega_p/2/3.14159); truncatedTau = 1;
    }

    // coefficients used exclusively in the braginskii transport
    Real delta_kappa, delta_eta, delta_eta2, x_coef ;

    x_coef = omega_ci*t_collision_ion;

    //TODO fix up coefficients here also with tabled depending atomic number

    delta_kappa = x_coef*x_coef*x_coef*x_coef + 2.700*x_coef*x_coef + 0.677;
    delta_eta   = x_coef*x_coef*x_coef*x_coef + 4.030*x_coef*x_coef + 2.330;
    delta_eta2  = 16*x_coef*x_coef*x_coef*x_coef + 4*4.030*x_coef*x_coef + 2.330;

    //assign viscosity value
    if (ion_coeffs.forceViscosity) eta0 = ion_coeffs.forceViscosityValue;
    else eta0 = 0.96*nd_i*T_i*t_collision_ion ;//* n0_ref;
    //Print() << "Ion eta0:\t" << eta0 << "\n";
    if (braginskii_anisotropic) {
        eta2 = nd_i*T_i*t_collision_ion*(6./5.*x_coef*x_coef+2.23)/delta_eta;
        eta1 = nd_i*T_i*t_collision_ion*(6./5.*(2*x_coef)*(2*x_coef)+2.23)/delta_eta2;
        eta4 = nd_i*T_i*t_collision_ion*x_coef*(x_coef*x_coef + 2.38)/delta_eta;
        eta3 = nd_i*T_i*t_collision_ion*(2*x_coef)*((2*x_coef)*(2*x_coef) + 2.38)/delta_eta2;
        if ((x_coef < 1e-8 ) && (omega_ci > effective_zero)) {
            eta2 = eta2/x_coef/x_coef;
            eta1 = eta1/x_coef/x_coef;
            eta4 = eta4/x_coef;
            eta3 = eta3/x_coef; }
    } else {
        eta2 = 0;
        eta1 = 0;
        eta4 = 0;
        eta3 = 0; }


    kappa1 = 3.906*nd_i*T_i*t_collision_ion/mass_i;
    if (braginskii_anisotropic) {
        kappa2 = (2.*x_coef*x_coef + 2.645)/delta_kappa*nd_i*T_i*t_collision_ion/mass_i;
        kappa3 = (5./2.*x_coef*x_coef + 4.65)*x_coef*nd_i*T_i*t_collision_ion/mass_i/delta_kappa;
    } else {kappa2 = 0; kappa3 = 0; }


    if ((kappa1<0) || (kappa2<0) || (kappa3 < 0)) {
        amrex::Abort("Braginski Ion coefficients are non-physical");
    }

    if (kappa1 < kappa2) {
        kappa2 = kappa1;
    }
    if (kappa1 < kappa3) {
        kappa3 = kappa1;
    }

    return;
}

Real BraginskiiCTU::get_max_speed_ions(const Vector<amrex::Real>& U_i,
                                       const Vector<amrex::Real>& U_e,
                                       const Vector<amrex::Real>& U_f
                                       )
{
    BL_PROFILE("BraginskiiIon::get_max_speed");

    HydroState &estate = *electron_state;
    HydroState &istate = *ion_state;

    //---Calculate the coefficients from scratch
    Real mass_i,mass_e,charge_i,charge_e,T_e,T_i,nd_i,nd_e;
    Real eta0, eta1, eta2, eta3, eta4, kappa1, kappa2, kappa3;


    //Extract and assign parameters from Q_i and Q_e
    //--- electron state and properties required for calcs -------Note move this
    charge_e= estate.gas->get_charge_from_cons(U_e); // electron propertis
    mass_e  = estate.gas->get_mass_from_cons(U_e);
    T_e     = estate.gas->get_temperature_from_cons(U_e);
    nd_e    = U_e[+HydroDef::ConsIdx::Density]/mass_e;

    //--- ion state and properties required for calcs
    charge_i= istate.gas->get_charge_from_cons(U_i);
    mass_i  = istate.gas->get_mass_from_cons(U_i);
    T_i     = istate.gas->get_temperature_from_cons(U_i); //
    nd_i    = U_i[+HydroDef::ConsIdx::Density]/mass_i;

    //Magnetic field
    Real Bx = U_f[+FieldDef::ConsIdx::Bx];
    Real By = U_f[+FieldDef::ConsIdx::By];
    Real Bz = U_f[+FieldDef::ConsIdx::Bz];


    // See page 215 (document numbering) of Braginskii's original transport paper
    Real t_collision_ion, p_lambda, omega_ci, omega_p;
    p_lambda = get_coulomb_logarithm(T_i,T_e,nd_e);
    // absence of the boltzmann constant due to usage of nondimensional variables.
    // Note that here charge_i = e^4*Z^4

    Real Debye = MFP::Debye, Larmor = MFP::Larmor, n0_ref=MFP::n0;

    t_collision_ion = std::pow(Debye,4)*n0_ref
            *(12*std::sqrt(mass_i)*std::pow(3.14159*T_i, 3./2.)) /
            (p_lambda * std::pow(charge_i,4) * nd_i);

    omega_ci = charge_i * std::sqrt( Bx*Bx + By*By + Bz*Bz ) / mass_i / Larmor;
    omega_p  = std::sqrt(nd_i*charge_i*charge_i/mass_i/Debye/Debye) ;

    if (1/t_collision_ion < effective_zero) t_collision_ion = 1/effective_zero;

    if (srin_switch && (1/t_collision_ion < omega_ci/10/2/3.14159) && (1/t_collision_ion < omega_p/10/2/3.14159)) {
        t_collision_ion = 1/std::min(omega_ci/2/3.14159, omega_p/2/3.14159) ;
    }

    //-- simple plasma with magnetic field
    Real delta_kappa, delta_eta, delta_eta2, x_coef ;

    x_coef = omega_ci*t_collision_ion;

    //TODO fix up coefficients here also with tabled depending atomic number
    delta_kappa = x_coef*x_coef*x_coef*x_coef + 2.700*x_coef*x_coef + 0.677;
    delta_eta   = x_coef*x_coef*x_coef*x_coef + 4.030*x_coef*x_coef + 2.330;
    delta_eta2  = 16*x_coef*x_coef*x_coef*x_coef + 4*4.030*x_coef*x_coef + 2.330;

    //assign viscosity value
    if (ion_coeffs.forceViscosity) eta0 = ion_coeffs.forceViscosityValue;
    else eta0 = 0.96*nd_i*T_i*t_collision_ion ;//* n0_ref;

    if (braginskii_anisotropic) {
        eta2 = nd_i*T_i*t_collision_ion*(6./5.*x_coef*x_coef+2.23)/delta_eta;
        eta1 = nd_i*T_i*t_collision_ion*(6./5.*(2*x_coef)*(2*x_coef)+2.23)/delta_eta2;
        eta4 = nd_i*T_i*t_collision_ion*x_coef*(x_coef*x_coef + 2.38)/delta_eta;
        eta3 = nd_i*T_i*t_collision_ion*(2*x_coef)*((2*x_coef)*(2*x_coef) + 2.38)/delta_eta2;
        if ((x_coef < 1e-8 ) && (omega_ci > effective_zero)) {
            eta2 = eta2/x_coef/x_coef;
            eta1 = eta1/x_coef/x_coef;
            eta4 = eta4/x_coef;
            eta3 = eta3/x_coef; }
    } else {
        eta2 = 0;
        eta1 = 0;
        eta4 = 0;
        eta3 = 0; }


    //From Braginskii OG paper page 250 of paper in journal heading 4
    // Kinetics of a simple plasma (Quantitative Analyis)
    kappa1 = 3.906*nd_i*T_i*t_collision_ion/mass_i;
    if (braginskii_anisotropic) {
        kappa2 = (2.*x_coef*x_coef + 2.645)/delta_kappa*nd_i*T_i*t_collision_ion/mass_i;
        kappa3 = (5./2.*x_coef*x_coef + 4.65)*x_coef*nd_i*T_i*t_collision_ion/mass_i/delta_kappa;
    } else {kappa2 = 0; kappa3 = 0; }


    if ((kappa1<0) || (kappa2<0) || (kappa3 < 0)) {
        amrex::Abort("Braginski Ion coefficients are non-physical");
    }

    if (kappa1 < kappa2) {

        kappa2 = kappa1;
    }
    if (kappa1 < kappa3) {
        kappa3 = kappa1;
    }

    Real rho = U_i[+HydroDef::ConsIdx::Density];
    Real cp_ion = istate.gas->get_cp_from_cons(U_i);
    Real nu_thermal = kappa1/rho/cp_ion/cfl; //thermal diffusivity
    Real nu_visc = (eta0/rho)/cfl;

    Real nu;
    if (nu_thermal> nu_visc) {
        nu = nu_thermal ;
    } else if (nu_thermal <= nu_visc) {
        nu = nu_visc ;
    }

    return 1.0;

}

void BraginskiiCTU::get_electron_coeffs(
        const Vector<Real>& Q_i,
        const Vector<Real>& Q_e,
        const Array<Real,3>& B_xyz,
        Real& T_e,
        Real& eta0,
        Real& eta1,
        Real& eta2,
        Real& eta3,
        Real& eta4,
        Real& kappa1,
        Real& kappa2,
        Real& kappa3,
        Real& beta1,
        Real& beta2,
        Real& beta3,
        int& truncatedTau
        )
{
    BL_PROFILE("BraginskiiEle::get_electron_coeffs");

    HydroState &estate = *electron_state;
    HydroState &istate = *ion_state;

    truncatedTau = 0;
    Real mass_i,mass_e,charge_i,charge_e,T_i,nd_i,nd_e;

    //Extract and assign parameters from Q_i and Q_e
    //--- electron state and properties required for calcs -------Note move this
    charge_i= istate.gas->get_charge_from_prim(Q_i); // electron propertis
    mass_i  = istate.gas->get_mass_from_prim(Q_i);
    T_i     = istate.gas->get_temperature_from_prim(Q_i);
    nd_i    = Q_i[+HydroDef::PrimIdx::Density]/mass_i;
    //--- ion state and properties required for calcs
    charge_e= estate.gas->get_charge_from_prim(Q_e);
    mass_e  = estate.gas->get_mass_from_prim(Q_e);
    T_e     = estate.gas->get_temperature_from_prim(Q_e); //
    nd_e    = Q_e[+HydroDef::PrimIdx::Density]/mass_e;
    //Magnetic field
    Real Bx=B_xyz[0], By=B_xyz[1], Bz=B_xyz[2];

    // See page 215 (document numbering) of Braginskii's original transport paper
    Real t_collision_ele, p_lambda, omega_ce, omega_p;
    p_lambda = get_coulomb_logarithm(T_i,T_e,nd_e);

    Real Debye = MFP::Debye, Larmor = MFP::Larmor, n0_ref=MFP::n0;

    t_collision_ele = std::pow(Debye,4)*n0_ref
            *(6*std::sqrt(2*mass_e)*std::pow(3.14159*T_e, 3./2.)) /
            (p_lambda*std::pow((charge_i/-charge_e), 2)*nd_i);

    omega_ce = -charge_e * std::sqrt( Bx*Bx + By*By + Bz*Bz ) / mass_e/ Larmor;
    omega_p  = std::sqrt(nd_e*charge_e*charge_e/mass_e/Debye/Debye) ;

    if (1/t_collision_ele < effective_zero) {
        t_collision_ele = 1/effective_zero;
    }

    if (srin_switch && (1/t_collision_ele < omega_ce/10/2/3.14159) && (1/t_collision_ele < omega_p/10/2/3.14159)) {
        t_collision_ele = 1/std::min(omega_ce/2/3.14159, omega_p/2/3.14159); truncatedTau = 1;
    }

    Real delta_kappa, delta_eta, delta_eta2, x_coef;// coefficients used exclusively in the braginskii
    x_coef = omega_ce*t_collision_ele;
    // TODO fix up these table 2 page 251 BT
    Real delta_0=3.7703, delta_1=14.79;
    delta_kappa= x_coef*x_coef*x_coef*x_coef+delta_1*x_coef*x_coef + delta_0;
    delta_eta  = x_coef*x_coef*x_coef*x_coef+13.8*x_coef*x_coef + 11.6;
    delta_eta2 = 16*x_coef*x_coef*x_coef*x_coef+4*13.8*x_coef*x_coef + 11.6;


    if (electron_coeffs.forceViscosity) eta0 = electron_coeffs.forceViscosityValue;
    else eta0 = 0.733*nd_e *T_e * t_collision_ele;
    if (braginskii_anisotropic) {
        eta2 = nd_e *T_e*t_collision_ele*(2.05*x_coef*x_coef+8.5)/delta_eta;
        eta1 = nd_e *T_e*t_collision_ele*(2.05*(2*x_coef)*(2*x_coef)+8.5)/delta_eta2;
        eta4 = -nd_e*T_e*t_collision_ele*x_coef*(x_coef*x_coef+7.91)/delta_eta;
        eta3 = -nd_e*T_e*t_collision_ele*(2*x_coef)*((2*x_coef)*(2*x_coef)+7.91)/delta_eta2;
        if ((x_coef < 1e-8 ) && (omega_ce > effective_zero)) {
            eta2 = eta2/x_coef/x_coef;
            eta1 = eta1/x_coef/x_coef;
            eta4 = eta4/x_coef;
            eta3 = eta3/x_coef; }
    } else {
        eta2 = 0;
        eta1 = 0;
        eta4 = 0;
        eta3 = 0; }



    //From Braginskii OG paper page 250 of paper in journal heading 4
    // Kinetics of a simple plasma (Quantitative Analyis)
    //TODO change coefficient values for different Z values
    // Currently set for a hydrogen  plasma
    Real BT_gamma_0=11.92/3.7703,BT_gamma_0_p=11.92,BT_gamma_1_p=4.664,BT_gamma_1_pp=5./2.;
    Real BT_gamma_0_pp=21.67;

    kappa1=nd_e*T_e*t_collision_ele/mass_e*BT_gamma_0;
    if (braginskii_anisotropic) {
        kappa2=(BT_gamma_1_p*x_coef*x_coef+BT_gamma_0_p)/delta_kappa*nd_e*T_e*t_collision_ele/mass_e;
        kappa3=(BT_gamma_1_pp*x_coef*x_coef+BT_gamma_0_pp)*x_coef*nd_e*T_e*t_collision_ele
                /mass_e/delta_kappa;
    } else {kappa2 = 0; kappa3 = 0;}


    if ((kappa1<0.) || (kappa2<0.) || (kappa3 < 0.)) {
        //amrex::Abort("mfp_viscous.cpp ln: 673 - Braginski Ion coefficients are non-physical");
        Print() << "\nmfp_viscous.cpp ln: 673 - Braginski Ion coefficients are non-physical\n";
        Print() << "\n" << kappa1 << "\n" << kappa2 << "\n" << kappa3 << "\nomega_ce = " << omega_ce;
    }

    if (kappa1 < kappa2) {
        kappa2 = kappa1;
    }
    if (kappa1 < kappa3) {
        kappa3 = kappa1;
    }

    //--- beta terms for the thermal component of thermal heat flux of the electrons.

    Real b_0 = 0.711, b_0_pp = 3.053, b_0_p=2.681, b_1_p=5.101, b_1_pp=3./2.;
    beta1 = nd_e*b_0*T_e;
    if (braginskii_anisotropic) {
        beta2 = nd_e*(b_1_p*x_coef*x_coef+b_0_p)/delta_kappa*T_e;
        beta3 = nd_e*x_coef*(b_1_pp*x_coef*x_coef+b_0_pp)/delta_kappa*T_e;
    } else {beta2 = 0; beta3 = 0;}

    return;
}

Real BraginskiiCTU::get_max_speed_electrons(const Vector<Real>& U_e,
                                            const Vector<Real>& U_i,
                                            const Vector<Real>& U_f
                                            )
{
    BL_PROFILE("BraginskiiEle::get_max_speed");

    HydroState &estate = *electron_state;
    HydroState &istate = *ion_state;

    //TODO Fix get_max_speed to be compatible with electrona and ion
    //---Calculate the coefficients from scratch
    Real mass_i,mass_e,charge_i,charge_e,T_e,T_i,nd_i,nd_e;

    //Extract and assign parameters from Q_i and Q_e

    //--- electron state and properties required for calcs -------Note move this
    charge_e= estate.gas->get_charge_from_cons(U_e); // electron propertis
    mass_e  = estate.gas->get_mass_from_cons(U_e);

    //Print() << "\nmass_i" << mass_i << "\n";
    T_e     = estate.gas->get_temperature_from_cons(U_e);
    const Real rho_e = U_e[+HydroDef::ConsIdx::Density];
    nd_e    = rho_e/mass_e;

    //--- ion state and properties required for calcs
    charge_i= istate.gas->get_charge_from_cons(U_i);
    mass_i  = istate.gas->get_mass_from_cons(U_i);
    //Print() << "\nmass_i" << mass_i << "\n";
    T_i     = istate.gas->get_temperature_from_cons(U_i); //
    const Real rho_i = U_i[+HydroDef::ConsIdx::Density];
    nd_i    = rho_i/mass_i;

    //Magnetic field
    Real Bx = U_f[+FieldDef::ConsIdx::Bx];
    Real By = U_f[+FieldDef::ConsIdx::By];
    Real Bz = U_f[+FieldDef::ConsIdx::Bz];

    Real t_collision_ele, p_lambda, omega_ce, omega_p;
    p_lambda = get_coulomb_logarithm(T_i,T_e,nd_e);

    //collision time nondimensional
    Real Debye = MFP::Debye, Larmor = MFP::Larmor, n0_ref=MFP::n0;

    t_collision_ele = std::pow(Debye,4)*n0_ref
            *(6*std::sqrt(2*mass_e)*std::pow(3.14159*T_e, 3./2.)) /
            (p_lambda*std::pow((charge_i/-charge_e), 2)*nd_i);
    // the issue with the 'c' value in the

    omega_ce = -charge_e * std::sqrt( Bx*Bx + By*By + Bz*Bz ) / mass_e/Larmor;
    omega_p  = std::sqrt(nd_e*charge_e*charge_e/mass_e/Debye/Debye) ;

    if (1/t_collision_ele < effective_zero) t_collision_ele = 1/effective_zero;

    if (srin_switch && (1/t_collision_ele < omega_ce/10/2/3.14159) && (1/t_collision_ele < omega_p/10/2/3.14159)) {
        t_collision_ele = 1/std::min(omega_ce/2/3.14159, omega_p/2/3.14159) ;
    }

    Real delta_kappa, delta_eta, delta_eta2,x_coef;// coefficients used exclusively in the braginskii
    x_coef = omega_ce*t_collision_ele;

    // TODO fix up these table 2 page 251 BT
    Real delta_0=3.7703, delta_1=14.79;
    delta_kappa= x_coef*x_coef*x_coef*x_coef+delta_1*x_coef*x_coef + delta_0;
    delta_eta  = x_coef*x_coef*x_coef*x_coef+13.8*x_coef*x_coef + 11.6;
    delta_eta2 = 16*x_coef*x_coef*x_coef*x_coef+4*13.8*x_coef*x_coef + 11.6;

    if (electron_coeffs.forceViscosity) electron_coeffs.eta0 = electron_coeffs.forceViscosityValue;
    else electron_coeffs.eta0 = 0.733*nd_e *T_e * t_collision_ele;

    if (braginskii_anisotropic) {
        electron_coeffs.eta2 = nd_e *T_e*t_collision_ele*(2.05*x_coef*x_coef+8.5)/delta_eta;
        electron_coeffs.eta1 = nd_e *T_e*t_collision_ele*(2.05*(2*x_coef)*(2*x_coef)+8.5)/delta_eta2;
        electron_coeffs.eta4 = -nd_e*T_e*t_collision_ele*x_coef*(x_coef*x_coef+7.91)/delta_eta;
        electron_coeffs.eta3 = -nd_e*T_e*t_collision_ele*(2*x_coef)*((2*x_coef)*(2*x_coef)+7.91)/delta_eta2;
        if ((x_coef < 1e-8 ) && (omega_ce > effective_zero)) {
            electron_coeffs.eta2 /= x_coef/x_coef;
            electron_coeffs.eta1 /= x_coef/x_coef;
            electron_coeffs.eta4 /= x_coef;
            electron_coeffs.eta3 /= x_coef; }
    } else {
        electron_coeffs.eta2 = 0;
        electron_coeffs.eta1 = 0;
        electron_coeffs.eta4 = 0;
        electron_coeffs.eta3 = 0; }

    //From Braginskii OG paper page 250 of paper in journal heading 4
    // Kinetics of a simple plasma (Quantitative Analyis)
    //TODO change coefficient values for different Z values
    // Currently set for a hydrogen  plasma
    Real BT_gamma_0=11.92/3.7703,BT_gamma_0_p=11.92,BT_gamma_1_p=4.664,BT_gamma_1_pp=5./2.;
    Real BT_gamma_0_pp=21.67;

    electron_coeffs.kappa1=nd_e*T_e*t_collision_ele/mass_e*BT_gamma_0;
    if (braginskii_anisotropic) {
        electron_coeffs.kappa2=(BT_gamma_1_p*x_coef*x_coef+BT_gamma_0_p)/delta_kappa*nd_e*T_e*t_collision_ele/mass_e;
        electron_coeffs.kappa3=(BT_gamma_1_pp*x_coef*x_coef+BT_gamma_0_pp)*x_coef*nd_e*T_e*t_collision_ele
                /mass_e/delta_kappa;}

    if ((electron_coeffs.kappa1<0.) || (electron_coeffs.kappa2<0.) || (electron_coeffs.kappa3 < 0.)) {
        //amrex::Abort("mfp_viscous.cpp ln: 673 - Braginski Ion coefficients are non-physical");
        Print() << "\nmfp_viscous.cpp ln: 673 - Braginski Ion coefficients are non-physical\n";
        Print() << "\n" << electron_coeffs.kappa1 << "\n" << electron_coeffs.kappa2 << "\n" << electron_coeffs.kappa3 << "\nomega_ce = " << omega_ce;
    }

    if (electron_coeffs.kappa1 < electron_coeffs.kappa2) {
        electron_coeffs.kappa2 = electron_coeffs.kappa1;
    }
    if (electron_coeffs.kappa1 < electron_coeffs.kappa3) {
        electron_coeffs.kappa3 = electron_coeffs.kappa1;
    }

    // add in the kappa and beta, i have tio figure out how to take these into account,
    // the viscous characteristic speed is not properly taken account of.
    Real cp_ele = estate.gas->get_cp_from_cons(U_e);
    Real nu_thermal = electron_coeffs.kappa1/rho_e/cp_ele/cfl; //thermal diffusivity
    Real nu_visc = (electron_coeffs.eta0/rho_e)/cfl;


    Real nu;
    if (nu_thermal> nu_visc) {
        nu = nu_thermal ;
    } else if (nu_thermal <= nu_visc) {
        nu = nu_visc ;
    }


    return 1.0;

    return 2*nu;
}


// ====================================================================================
void BraginskiiCTU::calc_ion_diffusion_terms(const Box& box,
                                             Array4<const Real> const& prim_i4,
                                             Array4<const Real> const& prim_e4,
                                             Array4<const Real> const& prim_f4,
                                             FArrayBox& diff
                                             ) {
    BL_PROFILE("BraginskiiCTU::calc_ion_diffusion_terms");

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<Real> const& d4 = diff.array();


    const int np_i = ion_state->n_prim();
    const int np_e = electron_state->n_prim();

    Vector<Real> Q_i(np_i), Q_e(np_e);
    Array<Real,3> B_xyz;
    //prefix of p_ denotes particle characteristic
    Real T_i, eta_0, eta_1, eta_2, eta_3, eta_4, kappa_1, kappa_2, kappa_3;
    int truncatedTau;
    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                for (int n=0; n<np_i; ++n) {
                    Q_i[n] = prim_i4(i,j,k,n);
                }
                for (int n=0; n<np_e; ++n) {
                    Q_e[n] = prim_e4(i,j,k,n);
                }

                B_xyz[0] = prim_f4(i,j,k,+FieldDef::ConsIdx::Bx);
                B_xyz[1] = prim_f4(i,j,k,+FieldDef::ConsIdx::By);
                B_xyz[2] = prim_f4(i,j,k,+FieldDef::ConsIdx::Bz);

                get_ion_coeffs(Q_i,Q_e,B_xyz,T_i,eta_0,eta_1,eta_2,eta_3,
                               eta_4, kappa_1,kappa_2, kappa_3, truncatedTau);
                //if the switch was used, notify
                //if (truncatedTau) Print() <<"\nion srin_switch used cell i, j, k: " << i << " " << j << " " << k << "\n";
                //assign values to the diff (d4) matrix for usage in the superior function

                d4(i,j,k,+IonDiffusionCoeffs::IonTemp) = T_i;
                d4(i,j,k,+IonDiffusionCoeffs::IonKappa1) = kappa_1;
                d4(i,j,k,+IonDiffusionCoeffs::IonKappa2) = kappa_2;
                d4(i,j,k,+IonDiffusionCoeffs::IonKappa3) = kappa_3;
                d4(i,j,k,+IonDiffusionCoeffs::IonEta0) = eta_0;
                d4(i,j,k,+IonDiffusionCoeffs::IonEta1) = eta_1;
                d4(i,j,k,+IonDiffusionCoeffs::IonEta2) = eta_2;
                d4(i,j,k,+IonDiffusionCoeffs::IonEta3) = eta_3;
                d4(i,j,k,+IonDiffusionCoeffs::IonEta4) = eta_4; // Note we could store the magnetic field but then we are doubling up on their storafe, perhpas better to just tolerate the access penalty
            }
        }
    }
    return;
}

void BraginskiiCTU::calc_ion_viscous_fluxes(const Box& box,
                                            Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                            const Box& pbox,
                                            Array4<const Real> const& prim_i4,
                                            Array4<const Real> const& prim_e4,
                                            Array4<const Real> const& prim_f4,
                                            const Real* dx) {
    BL_PROFILE("BraginskiiCTU::calc_ion_viscous_fluxes");
    //data strucutre for the diffusion coefficients.
    FArrayBox diff_ion(pbox, +IonDiffusionCoeffs::NUM_ION_DIFF_COEFFS);


    //--- diffusion coefficients for each cell to be used

    calc_ion_diffusion_terms(pbox,
                             prim_i4,
                             prim_e4,
                             prim_f4,
                             diff_ion);

    calc_charged_viscous_fluxes(FluxSpecies::IonFlux, box, fluxes,
                                prim_i4, prim_i4, prim_e4, prim_f4,
                                dx, diff_ion);
    return;
}

// ====================================================================================

void BraginskiiCTU::calc_electron_diffusion_terms(const Box& box,
                                                  Array4<const Real> const& prim_i4,
                                                  Array4<const Real> const& prim_e4,
                                                  Array4<const Real> const& prim_f4,
                                                  FArrayBox& diff
                                                  ) {
    BL_PROFILE("HydroState::calc_electron_diffusion_terms");

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<Real> const& d4 = diff.array();


    const int np_i = ion_state->n_prim();
    const int np_e = electron_state->n_prim();

    Vector<Real> Q_i(np_i), Q_e(np_e);
    Array<Real,3> B_xyz;
    //prefix of p_ denotes particle characteristic
    Real T_e, eta_0, eta_1, eta_2, eta_3, eta_4, kappa_1, kappa_2, kappa_3, beta1, beta2, beta3;
    int truncatedTau;
    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                for (int n=0; n<np_e; ++n) {
                    Q_e[n] = prim_e4(i,j,k,n);
                }
                for (int n=0; n<np_i; ++n) {
                    Q_i[n] = prim_i4(i,j,k,n);
                }

                B_xyz[0] = prim_f4(i,j,k,+FieldDef::ConsIdx::Bx);
                B_xyz[1] = prim_f4(i,j,k,+FieldDef::ConsIdx::By);
                B_xyz[2] = prim_f4(i,j,k,+FieldDef::ConsIdx::Bz);


                get_electron_coeffs(Q_i,Q_e,B_xyz,T_e,eta_0,eta_1,eta_2,eta_3,
                                    eta_4, kappa_1,kappa_2, kappa_3, beta1, beta2, beta3, truncatedTau);

                d4(i,j,k,+ElectronDiffusionCoeffs::EleTemp) = T_e;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleKappa1) = kappa_1;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleKappa2) = kappa_2;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleKappa3) = kappa_3;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleEta0) = eta_0;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleEta1) = eta_1;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleEta2) = eta_2;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleEta3) = eta_3;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleEta4) = eta_4;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleBeta1) = beta1;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleBeta2) = beta2;
                d4(i,j,k,+ElectronDiffusionCoeffs::EleBeta3) = beta3;
            }
        }
    }

    return;
}

void BraginskiiCTU::calc_electron_viscous_fluxes(const Box& box,
                                                 Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                                 const Box& pbox,
                                                 Array4<const Real> const& prim_i4,
                                                 Array4<const Real> const& prim_e4,
                                                 Array4<const Real> const& prim_f4,
                                                 const Real* dx) {
    BL_PROFILE("HydroState::calc_electron_viscous_fluxes");
    FArrayBox diff_ele(pbox, +ElectronDiffusionCoeffs::NUM_ELE_DIFF_COEFFS);


    //--- diffusion coefficients for each cell to be used
    calc_electron_diffusion_terms(pbox,
                                  prim_i4,
                                  prim_e4,
                                  prim_f4,
                                  diff_ele
                                  );
    //handle all the generic flux calculations

    calc_charged_viscous_fluxes(FluxSpecies::ElectronFlux,
                                box, fluxes,
                                prim_e4,
                                prim_i4,
                                prim_e4,
                                prim_f4,
                                dx, diff_ele);
    return;
}

// ====================================================================================

void BraginskiiCTU::calc_charged_viscous_fluxes(FluxSpecies flux_type,
                                                const Box& box,
                                                Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                                Array4<const Real> const& p4,
                                                Array4<const Real> const& prim_i4,
                                                Array4<const Real> const& prim_e4,
                                                Array4<const Real> const& prim_f4,
                                                const Real* dx, FArrayBox& diff)
{
    BL_PROFILE("HydroState::calc_charged_viscous_fluxes");

    //create a box for the viscous sress tensor where we only store the 6 unique
    // elements in order of 0:tauxx, 1:tauyy, 2:tauzz, 3:tauxy, 4:tauxz, 5:tauyz
    Array<Real,6> ViscTens;
    Array<Real,3> q_flux;
    //--- em state info

    //---Sorting out indexing and storage access
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<Real, AMREX_SPACEDIM> dxinv;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        dxinv[d] = 1/dx[d];
    }

    Array4<const Real> const& d4 = diff.array();

    //Array4<Real> const& d4 = diff.array();

    //Need to stage these properly once the algebraric expressions for the
    //viscous stresses are in
    Real dudx=0., dudy=0., dudz=0., dvdx=0., dvdy=0., dvdz=0., dwdx=0., dwdy=0.,
            dwdz=0., divu=0.;
    Real dTdx, dTdy, dTdz;

    //--------------Indexes for the trasnport coefficients
    int iTemp=-1, iEta0=-1, iEta1=-1, iEta2=-1, iEta3=-1, iEta4=-1, iKappa1=-1,
            iKappa2=-1, iKappa3=-1, iBeta1=-1, iBeta2=-1, iBeta3=-1;

    HydroState* state;

    if (flux_type == FluxSpecies::IonFlux) {
        iTemp = +IonDiffusionCoeffs::IonTemp;
        iEta0 = +IonDiffusionCoeffs::IonEta0;
        iEta1 = +IonDiffusionCoeffs::IonEta1;
        iEta2 = +IonDiffusionCoeffs::IonEta2;
        iEta3 = +IonDiffusionCoeffs::IonEta3;
        iEta4 = +IonDiffusionCoeffs::IonEta4;
        iKappa1=+IonDiffusionCoeffs::IonKappa1;
        iKappa2=+IonDiffusionCoeffs::IonKappa2;
        iKappa3 = +IonDiffusionCoeffs::IonKappa3;
        iTemp = +IonDiffusionCoeffs::IonTemp;
        state = ion_state;
    } else {
        iTemp = +ElectronDiffusionCoeffs::EleTemp;
        iEta0 = +ElectronDiffusionCoeffs::EleEta0;
        iEta1 = +ElectronDiffusionCoeffs::EleEta1;
        iEta2 = +ElectronDiffusionCoeffs::EleEta2;
        iEta3 = +ElectronDiffusionCoeffs::EleEta3;
        iEta4 = +ElectronDiffusionCoeffs::EleEta4;
        iKappa1= +ElectronDiffusionCoeffs::EleKappa1;
        iKappa2=+ElectronDiffusionCoeffs::EleKappa2;
        iKappa3= +ElectronDiffusionCoeffs::EleKappa3;
        iBeta1 =+ElectronDiffusionCoeffs::EleBeta1;
        iBeta2= +ElectronDiffusionCoeffs::EleBeta2;
        iBeta3 = +ElectronDiffusionCoeffs::EleBeta3;
        iTemp = +ElectronDiffusionCoeffs::EleTemp;
        state = electron_state;
    }




    constexpr int Xvel = +HydroDef::PrimIdx::Xvel;
    constexpr int Yvel = +HydroDef::PrimIdx::Yvel;
    constexpr int Zvel = +HydroDef::PrimIdx::Zvel;

    constexpr int Xmom = +HydroDef::ConsIdx::Xmom;
    constexpr int Ymom = +HydroDef::ConsIdx::Ymom;
    constexpr int Zmom = +HydroDef::ConsIdx::Zmom;

    constexpr int Eden = +HydroDef::ConsIdx::Eden;

    //Magnetic field components used for Braginskii terms, `p' represents prime.
    //Real bx_pp=0.,by_pp=0.,bz_pp=0.,bx_p=0.,by_p=0.,B=0.,B_pp=0.,B_p=0.;
    Real xB, yB, zB;

    Array<Real, +ElectronDiffusionCoeffs::NUM_ELE_DIFF_COEFFS> faceCoefficients; // use the beta places and leave them uninitiated if ion

    Vector<Real> u_rel(3);
    Array4<Real> const& fluxX = fluxes[0].array();
    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x + 1; ++i) {

                if (flux_type == FluxSpecies::ElectronFlux) {
                    faceCoefficients[iBeta1] = 0.5*(d4(i,j,k,iBeta1)+d4(i-1,j,k,iBeta1));
                    faceCoefficients[iBeta2] = 0.5*(d4(i,j,k,iBeta2)+d4(i-1,j,k,iBeta2));
                    faceCoefficients[iBeta3] = 0.5*(d4(i,j,k,iBeta3)+d4(i-1,j,k,iBeta3));
                }

                faceCoefficients[iKappa1] = 0.5*(d4(i,j,k,iKappa1)+d4(i-1,j,k,iKappa1));
                faceCoefficients[iKappa2] = 0.5*(d4(i,j,k,iKappa2)+d4(i-1,j,k,iKappa2));
                faceCoefficients[iKappa3] = 0.5*(d4(i,j,k,iKappa3)+d4(i-1,j,k,iKappa3));

                faceCoefficients[iEta0] = 0.5*(d4(i,j,k,iEta0)+d4(i-1,j,k,iEta0));
                faceCoefficients[iEta1] = 0.5*(d4(i,j,k,iEta1)+d4(i-1,j,k,iEta1));
                faceCoefficients[iEta2] = 0.5*(d4(i,j,k,iEta2)+d4(i-1,j,k,iEta2));
                faceCoefficients[iEta3] = 0.5*(d4(i,j,k,iEta3)+d4(i-1,j,k,iEta3));
                faceCoefficients[iEta4] = 0.5*(d4(i,j,k,iEta4)+d4(i-1,j,k,iEta4));

                xB = 0.5*(prim_f4(i,j,k,+FieldDef::ConsIdx::Bx) + prim_f4(i-1,j,k,+FieldDef::ConsIdx::Bx)); // using i j k  means you are taking the magnetic field in the cell i, not on the interface
                yB = 0.5*(prim_f4(i,j,k,+FieldDef::ConsIdx::By) + prim_f4(i-1,j,k,+FieldDef::ConsIdx::By));
                zB = 0.5*(prim_f4(i,j,k,+FieldDef::ConsIdx::Bz) + prim_f4(i-1,j,k,+FieldDef::ConsIdx::Bz));
                //if (global_idx == ion_idx)
                u_rel[0] = 0.5*(prim_e4(i,j,k,Xvel) + prim_e4(i-1,j,k,Xvel) - prim_i4(i,j,k,Xvel) - prim_i4(i-1,j,k,Xvel)); //TODO fix up flux
                u_rel[1] = 0.5*(prim_e4(i,j,k,Yvel) + prim_e4(i-1,j,k,Yvel) - prim_i4(i,j,k,Yvel) - prim_i4(i-1,j,k,Yvel));
                u_rel[2] = 0.5*(prim_e4(i,j,k,Zvel) + prim_e4(i-1,j,k,Zvel) - prim_i4(i,j,k,Zvel) - prim_i4(i-1,j,k,Zvel));

                dTdx = (d4(i,j,k,iTemp) - d4(i-1,j,k,iTemp))*dxinv[0];

                dudx = (p4(i,j,k,Xvel) - p4(i-1,j,k,Xvel))*dxinv[0];
                dvdx = (p4(i,j,k,Yvel) - p4(i-1,j,k,Yvel))*dxinv[0];
                dwdx = (p4(i,j,k,Zvel) - p4(i-1,j,k,Zvel))*dxinv[0];

#if AMREX_SPACEDIM >= 2
                dTdy = (d4(i,j+1,k,iTemp)+d4(i-1,j+1,k,iTemp)-d4(i,j-1,k,iTemp)-d4(i-1,j-1,k,iTemp))*(0.25*dxinv[1]);
                dudy = (p4(i,j+1,k,Xvel)+p4(i-1,j+1,k,Xvel)-p4(i,j-1,k,Xvel)-p4(i-1,j-1,k,Xvel))*(0.25*dxinv[1]);
                dvdy = (p4(i,j+1,k,Yvel)+p4(i-1,j+1,k,Yvel)-p4(i,j-1,k,Yvel)-p4(i-1,j-1,k,Yvel))*(0.25*dxinv[1]);

                // Put in to facilitate the matrix operations that will one day be
                // (//TODO) replaced with the explicit algebraic expression of the
                // viscous stress tensor entries hack without having correct
                // dimensional staging
                dwdy = (p4(i,j+1,k,Zvel)+p4(i-1,j+1,k,Zvel)-p4(i,j-1,k,Zvel)-p4(i-1,j-1,k,Zvel))*(0.25*dxinv[1]);
#endif

#if AMREX_SPACEDIM == 3
                dTdz = (d4(i,j,k+1,iTemp)+d4(i-1,j,k+1,iTemp)-d4(i,j,k-1,iTemp)-d4(i-1,j,k-1,iTemp))*(0.25*dxinv[1]);
                //(d4(i,j,k,iTemp)-d4(i,j,k-1,iTemp))*dxinv[2];

                dudz = (p4(i,j,k+1,Xvel)+p4(i-1,j,k+1,Xvel)-p4(i,j,k-1,Xvel)-p4(i-1,j,k-1,Xvel))*(0.25*dxinv[2]);

                dwdz = (p4(i,j,k+1,Zvel)+p4(i-1,j,k+1,Zvel)-p4(i,j,k-1,Zvel)-p4(i-1,j,k-1,Zvel))*(0.25*dxinv[2]);
                //Put in to hack without having correct dimensional staging
                dvdz = (p4(i,j,k+1,Yvel)+p4(i,j-1,k+1,Yvel)-p4(i,j,k-1,Yvel)-p4(i,j-1,k-1,Yvel))*(0.25*dxinv[2]);
#endif
                divu = dudx + dvdy + dwdz;

                ///TODO hacks because of transform which needs to be turned into algebra...

                //--- retrive the viscous stress tensor and heat flux vector on this face

                //TODO Print() << "Check changes to BRaginskiiViscousTensor... and calculation of xB, .. on interfaces are correct.";
                //Print() << "braginskii_anisotropic\t" << braginskii_anisotropic << "\n";
                if (braginskii_anisotropic) {
                    BraginskiiViscousTensorHeatFlux(flux_type,
                                                    xB, yB, zB, u_rel, dTdx, dTdy, dTdz,
                                                    dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz,
                                                    faceCoefficients, ViscTens, q_flux);
                } else {
                    IsotropicBraginskiiViscousTensorHeatFlux(flux_type,
                                                             u_rel, dTdx, dTdy, dTdz,
                                                             dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz,
                                                             faceCoefficients, ViscTens, q_flux);
                }


                fluxX(i,j,k,Xmom) += ViscTens[0];
                fluxX(i,j,k,Ymom) += ViscTens[3];
                fluxX(i,j,k,Zmom) += ViscTens[5];
                //assume typo in livescue formulation
                fluxX(i,j,k,Eden) += 0.5*((p4(i,j,k,Xvel) + p4(i-1,j,k,Xvel))*ViscTens[0]+
                        (p4(i,j,k,Yvel) + p4(i-1,j,k,Yvel))*ViscTens[3]+
                        (p4(i,j,k,Zvel) + p4(i-1,j,k,Zvel))*ViscTens[5])
                        + q_flux[0];

                AMREX_ASSERT(std::isfinite(fluxX(i,j,k,Xmom)));
                AMREX_ASSERT(std::isfinite(fluxX(i,j,k,Ymom)));
                AMREX_ASSERT(std::isfinite(fluxX(i,j,k,Zmom)));
                AMREX_ASSERT(std::isfinite(fluxX(i,j,k,Eden)));
            }
        }
    }



#if AMREX_SPACEDIM >= 2

    Array4<Real> const& fluxY = fluxes[1].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                if (flux_type == FluxSpecies::ElectronFlux) {

                    faceCoefficients[iBeta1] = 0.5*(d4(i,j,k,iBeta1)+d4(i,j-1,k,iBeta1));
                    faceCoefficients[iBeta2] = 0.5*(d4(i,j,k,iBeta2)+d4(i,j-1,k,iBeta2));
                    faceCoefficients[iBeta3] = 0.5*(d4(i,j,k,iBeta3)+d4(i,j-1,k,iBeta3));
                }


                faceCoefficients[iKappa1] = 0.5*(d4(i,j,k,iKappa1)+d4(i,j-1,k,iKappa1));
                faceCoefficients[iKappa2] = 0.5*(d4(i,j,k,iKappa2)+d4(i,j-1,k,iKappa2));
                faceCoefficients[iKappa3] = 0.5*(d4(i,j,k,iKappa3)+d4(i,j-1,k,iKappa3));

                faceCoefficients[iEta0] = 0.5*(d4(i,j,k,iEta0)+d4(i,j-1,k,iEta0));
                faceCoefficients[iEta1] = 0.5*(d4(i,j,k,iEta1)+d4(i,j-1,k,iEta1));
                faceCoefficients[iEta2] = 0.5*(d4(i,j,k,iEta2)+d4(i,j-1,k,iEta2));
                faceCoefficients[iEta3] = 0.5*(d4(i,j,k,iEta3)+d4(i,j-1,k,iEta3));
                faceCoefficients[iEta4] = 0.5*(d4(i,j,k,iEta4)+d4(i,j-1,k,iEta4));

                xB = 0.5*(prim_f4(i,j,k,+FieldDef::ConsIdx::Bx) + prim_f4(i,j-1,k,+FieldDef::ConsIdx::Bx)); // using i j k  means you are taking the magnetic field in the cell i, not on the interface
                yB = 0.5*(prim_f4(i,j,k,+FieldDef::ConsIdx::By) + prim_f4(i,j-1,k,+FieldDef::ConsIdx::By));
                zB = 0.5*(prim_f4(i,j,k,+FieldDef::ConsIdx::Bz) + prim_f4(i,j-1,k,+FieldDef::ConsIdx::Bz));
                u_rel[0] = 0.5*(prim_e4(i,j,k,Xvel) + prim_e4(i,j-1,k,Xvel) - prim_i4(i,j,k,Xvel) - prim_i4(i,j-1,k,Xvel)); //TODO fix up flux
                u_rel[1] = 0.5*(prim_e4(i,j,k,Yvel) + prim_e4(i,j-1,k,Yvel) - prim_i4(i,j,k,Yvel) - prim_i4(i,j-1,k,Yvel));
                u_rel[2] = 0.5*(prim_e4(i,j,k,Zvel) + prim_e4(i,j-1,k,Zvel) - prim_i4(i,j,k,Zvel) - prim_i4(i,j-1,k,Zvel));

                dTdy = (d4(i,j,k,iTemp)-d4(i,j-1,k,iTemp))*dxinv[1];

                dudy = (p4(i,j,k,Xvel)-p4(i,j-1,k,Xvel))*dxinv[1];
                dvdy = (p4(i,j,k,Yvel)-p4(i,j-1,k,Yvel))*dxinv[1];
                dwdy = (p4(i,j,k,Zvel)-p4(i,j-1,k,Zvel))*dxinv[1];

                dudx = (p4(i+1,j,k,Xvel)+p4(i+1,j-1,k,Xvel)-p4(i-1,j,k,Xvel)-p4(i-1,j-1,k,Xvel))*(0.25*dxinv[0]);
                dvdx = (p4(i+1,j,k,Yvel)+p4(i+1,j-1,k,Yvel)-p4(i-1,j,k,Yvel)-p4(i-1,j-1,k,Yvel))*(0.25*dxinv[0]);
                //--- retrive the viscous stress tensor and heat flux vector on this face
                ///TODO hacks because of transform which needs to be turned into algebra...
                dwdx = (p4(i+1,j,k,Zvel)+p4(i+1,j-1,k,Zvel)-p4(i-1,j,k,Zvel)-p4(i-1,j-1,k,Zvel))*(0.25*dxinv[0]);

#if AMREX_SPACEDIM == 3
                dvdz = (p4(i,j,k+1,Yvel)+p4(i,j-1,k+1,Yvel)-p4(i,j,k-1,Yvel)-p4(i,j-1,k-1,Yvel))*(0.25*dxinv[2]);
                dwdz = (p4(i,j,k+1,Zvel)+p4(i,j-1,k+1,Zvel)-p4(i,j,k-1,Zvel)-p4(i,j-1,k-1,Zvel))*(0.25*dxinv[2]);

                //--- retrive the viscous stress tensor and heat flux vector on this face
                ///TODO hacks because of transform which needs to be turned into algebra...
                dudz = (p4(i,j,k+1,Xvel)+p4(i,j-1,k+1,Xvel)-p4(i,j,k-1,Xvel)-p4(i,j-1,k-1,Xvel))*(0.25*dxinv[2]);

#endif
                divu = dudx + dvdy + dwdz;


                BraginskiiViscousTensorHeatFlux(flux_type,
                                                xB, yB, zB, u_rel, dTdx, dTdy, dTdz,
                                                dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz,
                                                faceCoefficients, ViscTens, q_flux);




                fluxY(i,j,k,Xmom) += ViscTens[3];//tauxy;
                fluxY(i,j,k,Ymom) += ViscTens[1];//tauyy;
                fluxY(i,j,k,Zmom) += ViscTens[4];//tauyz;
                fluxY(i,j,k,Eden) += +0.5*((p4(i,j,k,Xvel)+p4(i,j-1,k,Xvel))*ViscTens[3]+
                        (p4(i,j,k,Yvel)+p4(i,j-1,k,Yvel))*ViscTens[1]+
                        (p4(i,j,k,Zvel)+p4(i,j-1,k,Zvel))*ViscTens[4])
                        + q_flux[1];

                AMREX_ASSERT(std::isfinite(fluxY(i,j,k,Xmom)));
                AMREX_ASSERT(std::isfinite(fluxY(i,j,k,Ymom)));
                AMREX_ASSERT(std::isfinite(fluxY(i,j,k,Zmom)));
                AMREX_ASSERT(std::isfinite(fluxY(i,j,k,Eden)));
            }
        }
    }
#endif

#if AMREX_SPACEDIM == 3
    Array4<Real> const& fluxZ = fluxes[2].array();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y + 1; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                if (flux_type == ElectronFlux) {
                    faceCoefficients[iBeta1] = 0.5*(d4(i,j,k,iBeta1)+d4(i,j,k-1,iBeta1));
                    faceCoefficients[iBeta2] = 0.5*(d4(i,j,k,iBeta2)+d4(i,j,k-1,iBeta2));
                    faceCoefficients[iBeta3] = 0.5*(d4(i,j,k,iBeta3)+d4(i,j,k-1,iBeta3));
                }

                faceCoefficients[iKappa1] = 0.5*(d4(i,j,k,iKappa1)+d4(i,j,k-1,iKappa1));
                faceCoefficients[iKappa2] = 0.5*(d4(i,j,k,iKappa2)+d4(i,j,k-1,iKappa2));
                faceCoefficients[iKappa3] = 0.5*(d4(i,j,k,iKappa3)+d4(i,j,k-1,iKappa3));

                faceCoefficients[iEta0] = 0.5*(d4(i,j,k,iEta0)+d4(i,j,k-1,iEta0));
                faceCoefficients[iEta1] = 0.5*(d4(i,j,k,iEta1)+d4(i,j,k-1,iEta1));
                faceCoefficients[iEta2] = 0.5*(d4(i,j,k,iEta2)+d4(i,j,k-1,iEta2));
                faceCoefficients[iEta3] = 0.5*(d4(i,j,k,iEta3)+d4(i,j,k-1,iEta3));
                faceCoefficients[iEta4] = 0.5*(d4(i,j,k,iEta4)+d4(i,j,k-1,iEta4));

                xB = 0.5*(prim_f4(i,j,k,+FieldDef::ConsIdx::Bx) + prim_f4(i,j,k-1,+FieldDef::ConsIdx::Bx)); // using i j k  means you are taking the magnetic field in the cell i, not on the interface
                yB = 0.5*(prim_f4(i,j,k,+FieldDef::ConsIdx::By) + prim_f4(i,j,k-1,+FieldDef::ConsIdx::By));
                zB = 0.5*(prim_f4(i,j,k,+FieldDef::ConsIdx::Bz) + prim_f4(i,j,k-1,+FieldDef::ConsIdx::Bz));
                u_rel[0] = 0.5*(prim_e4(i,j,k,Xvel) + prim_e4(i,j,k-1,Xvel) - prim_i4(i,j,k,Xvel) - prim_i4(i,j,k-1,Xvel)); //TODO fix up flux
                u_rel[1] = 0.5*(prim_e4(i,j,k,Yvel) + prim_e4(i,j,k-1,Yvel) - prim_i4(i,j,k,Yvel) - prim_i4(i,j,k-1,Yvel));
                u_rel[2] = 0.5*(prim_e4(i,j,k,Zvel) + prim_e4(i,j,k-1,Zvel) - prim_i4(i,j,k,Zvel) - prim_i4(i,j,k-1,Zvel));

                dTdz = (d4(i,j,k,iTemp)-d4(i,j,k-1,iTemp))*dxinv[2];
                dudz = (p4(i,j,k,Xvel)-p4(i,j,k-1,Xvel))*dxinv[2];
                dvdz = (p4(i,j,k,Yvel)-p4(i,j,k-1,Yvel))*dxinv[2];
                dwdz = (p4(i,j,k,Zvel)-p4(i,j,k-1,Zvel))*dxinv[2];

                dudx = (p4(i+1,j,k,Xvel)+p4(i+1,j,k-1,Xvel)-p4(i-1,j,k,Xvel)-p4(i-1,j,k-1,Xvel))*(0.25*dxinv[0]);
                dwdx = (p4(i+1,j,k,Zvel)+p4(i+1,j,k-1,Zvel)-p4(i-1,j,k,Zvel)-p4(i-1,j,k-1,Zvel))*(0.25*dxinv[0]);
                dvdy = (p4(i,j+1,k,Yvel)+p4(i,j+1,k-1,Yvel)-p4(i,j-1,k,Yvel)-p4(i,j-1,k-1,Yvel))*(0.25*dxinv[1]);
                dwdy = (p4(i,j+1,k,Zvel)+p4(i,j+1,k-1,Zvel)-p4(i,j-1,k,Zvel)-p4(i,j-1,k-1,Zvel))*(0.25*dxinv[1]);
                divu = dudx + dvdy + dwdz;

                ///TODO hacks because of transform which needs to be turned into algebra...
                dvdx = (p4(i+1,j,k,Yvel)+p4(i+1,j,k-1,Yvel)-p4(i-1,j,k,Yvel)-p4(i-1,j,k-1,Yvel))*(0.25*dxinv[0]);
                dudy = (p4(i,j+1,k,Xvel)+p4(i,j+1,k-1,Xvel)-p4(i,j-1,k,Xvel)-p4(i,j-1,k-1,Xvel))*(0.25*dxinv[1]);

                //--- retrive the viscous stress tensor and heat flux vector on this face
                BraginskiiViscousTensorHeatFlux(flux_type,xB, yB, zB, u_rel, dTdx, dTdy, dTdz,
                                                dudx, dudy, dudz, dvdx, dvdy, dvdz,
                                                dwdx, dwdy, dwdz, faceCoefficients, ViscTens, q_flux);

                fluxZ(i,j,k,Xmom) += ViscTens[5];//tauxz;
                fluxZ(i,j,k,Ymom) += ViscTens[4];//tauyz;
                fluxZ(i,j,k,Zmom) += ViscTens[2];//tauzz;
                fluxZ(i,j,k,Eden) += +0.5*((p4(i,j,k,Xvel)+p4(i,j,k-1,Xvel))*ViscTens[5]+
                        (p4(i,j,k,Yvel)+p4(i,j,k-1,Yvel))*ViscTens[4]+
                        (p4(i,j,k,Zvel)+p4(i,j,k-1,Zvel))*ViscTens[2])
                        + q_flux[2];

                AMREX_ASSERT(std::isfinite(fluxZ(i,j,k,Xmom)));
                AMREX_ASSERT(std::isfinite(fluxZ(i,j,k,Ymom)));
                AMREX_ASSERT(std::isfinite(fluxZ(i,j,k,Zmom)));
                AMREX_ASSERT(std::isfinite(fluxZ(i,j,k,Eden)));

            }
        }
    }

#endif

    return;
}

// ===================================================================================
// function for calculating the viscous stress tensor and heat flux vector
// for the braginskii transport
void BraginskiiCTU::BraginskiiViscousTensorHeatFlux(FluxSpecies flux_type,
                                                    Real xB, Real yB, Real zB,
                                                    Vector<Real> u_rel,
                                                    Real dTdx, Real dTdy, Real dTdz,
                                                    Real dudx, Real dudy, Real dudz,
                                                    Real dvdx, Real dvdy, Real dvdz,
                                                    Real dwdx, Real dwdy, Real dwdz,
                                                    Array<Real,+ElectronDiffusionCoeffs::NUM_ELE_DIFF_COEFFS> faceCoefficients,
                                                    Array<Real, 6> &ViscTens,
                                                    Array<Real, 3> &q_flux) {
    BL_PROFILE("HydroState::BraginskiiViscousTensorHeatFlux");
    //Print() << "braginskii anisotropic function\n";
    //Note all the properties used in here need to be for the interface, not just the cell i!!!

    //---Sorting out indexing and storage access

    Real bx_pp=0.,by_pp=0.,bz_pp=0.,bx_p=0.,by_p=0.,B=0.,B_pp=0.,B_p=0.;
    Real qu_e_temp_max=0., qt_e_temp_max=0., qt_i_temp_max=0.;

    Real divu = dudx + dvdy + dwdz ;
    int i_disp, j_disp, k_disp;
    //Using the formulation of Li 2018 "High-order two-fluid plasma
    // solver for direct numerical simulations of plasma flowswith full
    // transport phenomena"  --- this is outdated, i think i was eefering
    // to the coulomb loagrithm which i just ended up taking fro braginskii

    //---get magnetic field orientated unit vectors
    B = xB*xB + yB*yB + zB*zB;

    if (B < effective_zero) {
        if (MFP::verbosity >= 3) {
            Print() << "\nZero magnetic field \n";
        }
        //amrex::Warning("zero magnetic field check the viscous stress matrix is not transformed at all.");
        B_pp = 0.;
        B_p  = 0.;
    } else if ( (std::abs(xB) < effective_zero) && (std::abs(yB) < effective_zero) &&
                (std::abs(zB) > effective_zero) ) {
        if (MFP::verbosity >= 4) {
            Print() << "\nZero x and y magnetic field \n";
            amrex::Warning("fixed frame aligned mag field, check the viscous stress matrix is not transformed at all.");
        }
        B_pp = 1/sqrt(B); // B prime prime
        B_p  = 0.;
    } else {
        B_pp = 1/sqrt(B); // B prime prime
        B_p  = 1/sqrt(xB*xB + yB*yB); // B prime
    }


    bx_pp = xB*B_pp; bx_p = xB*B_p;
    by_pp = yB*B_pp; by_p = yB*B_p;
    bz_pp = zB*B_pp;

    Array<Real,3> B_unit; //Magnetic field unit vector
    Array<Real,3> u_para; //Velocity parallel to B_unit
    Array<Real,3> u_perp; //Velocity perpendicular to B_unit
    Array<Real,3> u_chev; //unit vector perp to u and B_unit
    Array<Real,3> TG_para;//Temp grad parallel to B_unit
    Array<Real,3> TG_perp;//Temp grad perpendicular to B_unit
    Array<Real,3> TG_chev;//unit vector perp to gradT and B_unit

    B_unit[0] = bx_pp; B_unit[1] = by_pp; B_unit[2] = bz_pp;

    int iTemp=-1, iEta0=-1, iEta1=-1, iEta2=-1, iEta3=-1, iEta4=-1, iKappa1=-1,
            iKappa2=-1, iKappa3=-1, iBeta1=-1, iBeta2=-1, iBeta3=-1;

    if (flux_type == FluxSpecies::IonFlux) {
        iTemp = +IonDiffusionCoeffs::IonTemp;
        iEta0 = +IonDiffusionCoeffs::IonEta0;
        iEta1 = +IonDiffusionCoeffs::IonEta1;
        iEta2 = +IonDiffusionCoeffs::IonEta2;
        iEta3 = +IonDiffusionCoeffs::IonEta3;
        iEta4 = +IonDiffusionCoeffs::IonEta4;
        iKappa1=+IonDiffusionCoeffs::IonKappa1;
        iKappa2=+IonDiffusionCoeffs::IonKappa2;
        iKappa3 = +IonDiffusionCoeffs::IonKappa3;

    } else {
        iTemp = +ElectronDiffusionCoeffs::EleTemp;
        iEta0 = +ElectronDiffusionCoeffs::EleEta0;
        iEta1 = +ElectronDiffusionCoeffs::EleEta1;
        iEta2 = +ElectronDiffusionCoeffs::EleEta2;
        iEta3 = +ElectronDiffusionCoeffs::EleEta3;
        iEta4 = +ElectronDiffusionCoeffs::EleEta4;
        iKappa1=+ElectronDiffusionCoeffs::EleKappa1;
        iKappa2=+ElectronDiffusionCoeffs::EleKappa2;
        iKappa3= +ElectronDiffusionCoeffs::EleKappa3;
        iBeta1 =+ElectronDiffusionCoeffs::EleBeta1;
        iBeta2= +ElectronDiffusionCoeffs::EleBeta2;
        iBeta3 = +ElectronDiffusionCoeffs::EleBeta3;

    }

    Real dot_B_unit_TG, dot_B_unit_U ;//temp variables

    dot_B_unit_U = bx_pp*u_rel[0] + by_pp*u_rel[1] + bz_pp*u_rel[2];
    dot_B_unit_TG= bx_pp*dTdx + by_pp*dTdy + bz_pp*dTdz;
    // Prepare the x-components components
    u_para[0] = B_unit[0]*dot_B_unit_U ;
    TG_para[0]= B_unit[0]*dot_B_unit_TG ;
    u_perp[0] = u_rel[0] - u_para[0];
    u_chev[0] = B_unit[1]*u_rel[2] - B_unit[2]*u_rel[1];
    TG_perp[0]= dTdx - TG_para[0];
    TG_chev[0]= B_unit[1]*dTdz-B_unit[2]*dTdy;
    // Prepare the y-components components
    u_para[1] = B_unit[1]*dot_B_unit_U ;
    TG_para[1]= B_unit[1]*dot_B_unit_TG ;
    u_perp[1] = u_rel[1] - u_para[1];
    u_chev[1] = -(B_unit[0]*u_rel[2]-B_unit[2]*u_rel[0]);
    TG_perp[1]= dTdy - TG_para[1];
    TG_chev[1]= -(B_unit[0]*dTdz-B_unit[2]*dTdx);
    // Prepare the z-components components
    u_para[2] = B_unit[2]*dot_B_unit_U ;
    TG_para[2]= B_unit[2]*dot_B_unit_TG ;
    u_perp[2] = u_rel[2] - u_para[2];
    u_chev[2] = B_unit[0]*u_rel[1]-B_unit[1]*u_rel[0];
    TG_perp[2]= dTdz - TG_para[2];
    TG_chev[2]= B_unit[0]*dTdy-B_unit[1]*dTdx;

    int rowFirst=3,columnFirst=3,columnSecond=3;
    /* Matrix representations of viscous stree tensor and associated
    quantities are defined as:
    Trans       - the transformation matrix Q
    Strain      - Strain rate matrix mate W
    StrainTrans - Strain rate matrix transformed into the B aligned frame
    ViscStress  - Viscous stress tensor in lab frame, PI
    ViscStressTrans - Viscous stress tensor in B aligned frame, PI'
    TransT      - Transpose of the transformation matrix, Q'
    */
    Array<Array<Real,3>,3> Trans, Strain, StrainTrans, ViscStress, ViscStressTrans, TransT, WorkingMatrix;

    //if (global_idx == electron_idx)
    if (flux_type == FluxSpecies::ElectronFlux) {
        const Real qu_e_temp = faceCoefficients[iBeta1]*u_para[0] + faceCoefficients[iBeta2]*u_perp[0] + faceCoefficients[iBeta3]*u_chev[0] ;
        const Real qt_e_temp = -faceCoefficients[iKappa1]*TG_para[0] - faceCoefficients[iKappa2]*TG_perp[0] - faceCoefficients[iKappa3]*TG_chev[0];

        qu_e_temp_max = std::max(std::fabs(qu_e_temp_max), std::fabs(qu_e_temp));
        qt_e_temp_max = std::max(std::fabs(qt_e_temp_max), std::fabs(qt_e_temp));

        q_flux[0] = qt_e_temp + qu_e_temp ;

        q_flux[1] = faceCoefficients[iBeta1]*u_para[1] + faceCoefficients[iBeta2]*u_perp[1] + faceCoefficients[iBeta3]*u_chev[1]
                -faceCoefficients[iKappa1]*TG_para[1] - faceCoefficients[iKappa2]*TG_perp[1] - faceCoefficients[iKappa3]*TG_chev[1];

        q_flux[2] = faceCoefficients[iBeta1]*u_para[2] + faceCoefficients[iBeta2]*u_perp[2] + faceCoefficients[iBeta3]*u_chev[2]
                -faceCoefficients[iKappa1]*TG_para[2] - faceCoefficients[iKappa2]*TG_perp[2] - faceCoefficients[iKappa3]*TG_chev[2];

    } else {
        q_flux[0] = -faceCoefficients[iKappa1]*TG_para[0] - faceCoefficients[iKappa2]*TG_perp[0] + faceCoefficients[iKappa3]*TG_chev[0];

        qt_i_temp_max=std::max(std::fabs(qt_i_temp_max),std::fabs(q_flux[0]));

        q_flux[1] = -faceCoefficients[iKappa1]*TG_para[1] - faceCoefficients[iKappa2]*TG_perp[1] + faceCoefficients[iKappa3]*TG_chev[1];



        q_flux[2] = -faceCoefficients[iKappa1]*TG_para[2] - faceCoefficients[iKappa2]*TG_perp[2] + faceCoefficients[iKappa3]*TG_chev[2];

    }

    //Calculate the viscous stress tensor
    //Populate strain rate tensor in B unit aligned cartesian frame
    //This is braginskii's  - the strain rate tensor multplied by negative one later to Li
    // Livescue formulation
    Strain[0][0] = 2*dudx - 2./3.*divu;
    Strain[0][1] = dudy + dvdx;
    Strain[0][2] = dwdx + dudz;
    Strain[1][0] = Strain[0][1];
    Strain[1][1] = 2*dvdy - 2./3.*divu;
    Strain[1][2] = dvdz + dwdy;
    Strain[2][0] = Strain[0][2];
    Strain[2][1] = Strain[1][2];
    Strain[2][2] = 2*dwdz - 2./3.*divu;

    for (i_disp=0; i_disp<3; ++i_disp) { // set to zero
        for (j_disp=0; j_disp<3; ++j_disp) {
            WorkingMatrix[i_disp][j_disp] = 0.;
            StrainTrans[i_disp][j_disp] = 0.;
        }
    }

    // Do we have a special case of B=0 or Bx=By=0 and Bz = B?
    if (B < effective_zero) { // B = 0

        ViscStress[0][0] = -faceCoefficients[iEta0]*Strain[0][0];

        ViscStress[0][1] = -faceCoefficients[iEta0]*Strain[0][1];
        ViscStress[1][0] = ViscStress[0][1];

        ViscStress[0][2] = -faceCoefficients[iEta0]*Strain[0][2];
        ViscStress[2][0] = ViscStress[0][2];

        ViscStress[1][1] = -faceCoefficients[iEta0]*Strain[1][1];
        ViscStress[1][2] = -faceCoefficients[iEta0]*Strain[1][2];
        ViscStress[2][1] = ViscStress[1][2];

        ViscStress[2][2] = -faceCoefficients[iEta0]*Strain[2][2];

    } else if ( (std::abs(xB) < effective_zero) && (std::abs(yB) < effective_zero) &&
                (std::abs(zB) > effective_zero) ) { //B aligned with z direction
        ViscStress[0][0]=-1/2*faceCoefficients[iEta0]*
                (Strain[0][0] + Strain[1][1])
                -1/2*faceCoefficients[iEta1]*
                (Strain[0][0] - Strain[1][1])
                -faceCoefficients[iEta3]*(Strain[0][1]);

        ViscStress[0][1]=-faceCoefficients[iEta1]*Strain[0][1]
                +1/2*faceCoefficients[iEta3]*
                (Strain[0][0] - Strain[1][1]);

        ViscStress[1][0]= ViscStress[0][1];

        ViscStress[0][2]=-faceCoefficients[iEta2]*Strain[0][2]
                - faceCoefficients[iEta4]*Strain[1][2];

        ViscStress[2][0]= ViscStress[0][2];

        ViscStress[1][1]=-1/2*faceCoefficients[iEta0]*
                (Strain[0][0] + Strain[1][1])
                -1/2*faceCoefficients[iEta1]*
                (Strain[1][1] - Strain[0][0])
                +faceCoefficients[iEta3]*Strain[0][1];

        ViscStress[1][2]=-faceCoefficients[iEta2]*Strain[1][2] +
                faceCoefficients[iEta4]*Strain[0][2];

        ViscStress[2][1]=ViscStress[1][2];

        ViscStress[2][2]=-faceCoefficients[iEta0]*Strain[2][2];

    } else { //the generic case with non trivial magnetic field requiring a transform
        //Populate the transformation matrix from cartesian normal to B unit
        // aligned cartesian - Li 2018

        //Print() << "braginskii anisotropic matrix transforms executed\n";
        Trans[0][0]=-by_p; Trans[0][1]=-bx_p*bz_pp; Trans[0][2]=bx_pp;

        Trans[1][0]= bx_p; Trans[1][1]=-by_p*bz_pp; Trans[1][2] = by_pp;

        Trans[2][0]= 0;    Trans[2][1]= bx_p*bx_pp
                +by_p*by_pp; Trans[2][2] =bz_pp;
        //Populate the transpose of the transformation matrix
        TransT[0][0]=Trans[0][0]; TransT[0][1]=Trans[1][0]; TransT[0][2]=Trans[2][0];

        TransT[1][0]=Trans[0][1]; TransT[1][1]=Trans[1][1]; TransT[1][2]=Trans[2][1];

        TransT[2][0]=Trans[0][2]; TransT[2][1]=Trans[1][2]; TransT[2][2]=Trans[2][2];

        // Multiplying Q' (Transpose) by W (StressStrain)
        for (i_disp = 0; i_disp< rowFirst; ++i_disp) {
            for(j_disp = 0; j_disp< columnSecond; ++j_disp) {
                for(k_disp=0; k_disp<columnFirst; ++k_disp) {
                    WorkingMatrix[i_disp][j_disp] += TransT[i_disp][k_disp]*
                            Strain[k_disp][j_disp];
                }
            }
        }
        // Multiplying Q'W by Q
        for(i_disp = 0; i_disp< rowFirst; ++i_disp) {
            for(j_disp = 0; j_disp< columnSecond; ++j_disp) {
                for(k_disp=0; k_disp<columnFirst; ++k_disp) {
                    StrainTrans[i_disp][j_disp] += WorkingMatrix[i_disp][k_disp] *
                            Trans[k_disp][j_disp];
                }
            }
        }
        //faceCoefficients[iKappa3]
        //Populate visc stress tensor in cartesian normal frame
        ViscStressTrans[0][0]=-1/2*faceCoefficients[iEta0]*
                (StrainTrans[0][0] + StrainTrans[1][1])
                -1/2*faceCoefficients[iEta1]*
                (StrainTrans[0][0] - StrainTrans[1][1])
                -faceCoefficients[iEta3]*(StrainTrans[0][1]);
        ViscStressTrans[0][1]=-faceCoefficients[iEta1]*StrainTrans[0][1]
                +1/2*faceCoefficients[iEta3]*
                (StrainTrans[0][0] - StrainTrans[1][1]);

        ViscStressTrans[1][0]= ViscStressTrans[0][1];

        ViscStressTrans[0][2]=-faceCoefficients[iEta2]*StrainTrans[0][2]
                - faceCoefficients[iEta4]*StrainTrans[1][2];

        ViscStressTrans[2][0]= ViscStressTrans[0][2];

        ViscStressTrans[1][1]=-1/2*faceCoefficients[iEta0]*
                (StrainTrans[0][0] + StrainTrans[1][1])
                -1/2*faceCoefficients[iEta1]*
                (StrainTrans[1][1] - StrainTrans[0][0])
                +faceCoefficients[iEta3]*StrainTrans[0][1];

        ViscStressTrans[1][2]=-faceCoefficients[iEta2]*StrainTrans[1][2] +
                faceCoefficients[iEta4]*StrainTrans[0][2];

        ViscStressTrans[2][1]=ViscStressTrans[1][2];

        ViscStressTrans[2][2]=-faceCoefficients[iEta0]*StrainTrans[2][2];

        for (i_disp=0; i_disp<3; ++i_disp) {//Set to zero
            for (j_disp=0; j_disp<3; ++j_disp) {
                WorkingMatrix[i_disp][j_disp] = 0.;
                ViscStress[i_disp][j_disp] = 0.;
            }
        }
        // Multiplying Q (Trans) by PI' (ViscStressTrans)
        for(i_disp = 0; i_disp< rowFirst; ++i_disp) {
            for(j_disp = 0; j_disp< columnSecond; ++j_disp) {
                for(k_disp=0; k_disp<columnFirst; ++k_disp) {
                    WorkingMatrix[i_disp][j_disp] += Trans[i_disp][k_disp] *
                            ViscStressTrans[k_disp][j_disp];
                }
            }
        }
        // Multiplying Q*PI' by Q^T
        for(i_disp = 0; i_disp< rowFirst; ++i_disp) {
            for(j_disp = 0; j_disp< columnSecond; ++j_disp) {
                for(k_disp=0; k_disp<columnFirst; ++k_disp) {
                    ViscStress[i_disp][j_disp] += WorkingMatrix[i_disp][k_disp] *
                            TransT[k_disp][j_disp];
                }
            }
        }
    }

    //Storing
    //NOTE STORAGE ACCORDING TO THE TAUXX, TAUYY, TAUZZ, TAUXY OR TAUYX,
    // TAUYZ OR TAUZY, TAUXZ OR TAUZX
    ViscTens[0] = ViscStress[0][0];
    ViscTens[1] = ViscStress[1][1];
    ViscTens[2] = ViscStress[2][2];
    ViscTens[3] = ViscStress[0][1];
    ViscTens[4] = ViscStress[1][2];
    ViscTens[5] = ViscStress[0][2];

    return ;
}


// ===================================================================================
// function for calculating the viscous stress tensor and heat flux vector
// for the braginskii transport
void BraginskiiCTU::IsotropicBraginskiiViscousTensorHeatFlux(FluxSpecies flux_type,
                                                             Vector<Real> u_rel,
                                                             Real dTdx, Real dTdy, Real dTdz,
                                                             Real dudx, Real dudy, Real dudz,
                                                             Real dvdx, Real dvdy, Real dvdz,
                                                             Real dwdx, Real dwdy, Real dwdz,
                                                             Array<Real,+ElectronDiffusionCoeffs::NUM_ELE_DIFF_COEFFS> faceCoefficients,
                                                             Array<Real, 6> &ViscTens,
                                                             Array<Real, 3> &q_flux) {
    BL_PROFILE("BraginskiiCTU::IsotropicBraginskiiViscousTensorHeatFlux");
    //Note all the properties used in here need to be for the interface, not just the cell i!!!


    Real divu = dudx + dvdy + dwdz ;
    int i_disp, j_disp;
    //Using the formulation of Li 2018 "High-order two-fluid plasma
    // solver for direct numerical simulations of plasma flowswith full
    // transport phenomena"  --- this is outdated, i think i was eefering
    // to the coulomb loagrithm which i just ended up taking fro braginskii

    int iTemp=-1, iEta0=-1, iEta1=-1, iEta2=-1, iEta3=-1, iEta4=-1, iKappa1=-1,
            iKappa2=-1, iKappa3=-1, iBeta1=-1, iBeta2=-1, iBeta3=-1;

    if (flux_type == FluxSpecies::IonFlux) {
        iTemp = +IonDiffusionCoeffs::IonTemp;
        iEta0 = +IonDiffusionCoeffs::IonEta0;
        iEta1 = +IonDiffusionCoeffs::IonEta1;
        iEta2 = +IonDiffusionCoeffs::IonEta2;
        iEta3 = +IonDiffusionCoeffs::IonEta3;
        iEta4 = +IonDiffusionCoeffs::IonEta4;
        iKappa1=+IonDiffusionCoeffs::IonKappa1;
        iKappa2=+IonDiffusionCoeffs::IonKappa2;
        iKappa3 = +IonDiffusionCoeffs::IonKappa3;

    } else {
        iTemp = +ElectronDiffusionCoeffs::EleTemp;
        iEta0 = +ElectronDiffusionCoeffs::EleEta0;
        iEta1 = +ElectronDiffusionCoeffs::EleEta1;
        iEta2 = +ElectronDiffusionCoeffs::EleEta2;
        iEta3 = +ElectronDiffusionCoeffs::EleEta3;
        iEta4 = +ElectronDiffusionCoeffs::EleEta4;
        iKappa1=+ElectronDiffusionCoeffs::EleKappa1;
        iKappa2=+ElectronDiffusionCoeffs::EleKappa2;
        iKappa3= +ElectronDiffusionCoeffs::EleKappa3;
        iBeta1 =+ElectronDiffusionCoeffs::EleBeta1;
        iBeta2= +ElectronDiffusionCoeffs::EleBeta2;
        iBeta3 = +ElectronDiffusionCoeffs::EleBeta3;

    }


    /* Matrix representations of viscous stree tensor and associated
    quantities are defined as:
    Strain      - Strain rate matrix mate W
    ViscStress  - Viscous stress tensor in lab frame, PI
    */
    Array<Array<Real,3>,3> Strain;
    Array<Array<Real,3>,3> ViscStress;

    //if (global_idx == electron_idx)
    if (flux_type == FluxSpecies::ElectronFlux) {
        q_flux[0] = faceCoefficients[iBeta1]*u_rel[0] - faceCoefficients[iKappa1]*dTdx;
        q_flux[1] = faceCoefficients[iBeta1]*u_rel[1] - faceCoefficients[iKappa1]*dTdy;
        q_flux[2] = faceCoefficients[iBeta1]*u_rel[2] - faceCoefficients[iKappa1]*dTdz;
    } else {
        q_flux[0] = -faceCoefficients[iKappa1]*dTdx;
        q_flux[1] = -faceCoefficients[iKappa1]*dTdy;
        q_flux[2] = -faceCoefficients[iKappa1]*dTdz;
    }

    //Calculate the viscous stress tensor
    //Populate strain rate tensor in B unit aligned cartesian frame
    //This is braginskii's  - the strain rate tensor multplied by negative one later to Li
    // Livescue formulation
    Strain[0][0] = 2*dudx - 2./3.*divu;
    Strain[0][1] = dudy + dvdx;
    Strain[0][2] = dwdx + dudz;
    Strain[1][0] = Strain[0][1];
    Strain[1][1] = 2*dvdy - 2./3.*divu;
    Strain[1][2] = dvdz + dwdy;
    Strain[2][0] = Strain[0][2];
    Strain[2][1] = Strain[1][2];
    Strain[2][2] = 2*dwdz - 2./3.*divu;

    for (i_disp=0; i_disp<3; ++i_disp) { // set to zero
        for (j_disp=0; j_disp<3; ++j_disp) {
            Strain[i_disp][j_disp] = - Strain[i_disp][j_disp] ; //make negtive for Li Livescue
        }
    }
    ViscStress[0][0] = faceCoefficients[iEta0]*Strain[0][0];

    ViscStress[0][1] = faceCoefficients[iEta0]*Strain[0][1];
    ViscStress[1][0] = ViscStress[0][1];

    ViscStress[0][2] = faceCoefficients[iEta0]*Strain[0][2];
    ViscStress[2][0] = ViscStress[0][2];

    ViscStress[1][1] = faceCoefficients[iEta0]*Strain[1][1];
    ViscStress[1][2] = faceCoefficients[iEta0]*Strain[1][2];
    ViscStress[2][1] = ViscStress[1][2];

    ViscStress[2][2] = faceCoefficients[iEta0]*Strain[2][2];

    //Storing
    //NOTE STORAGE ACCORDING TO THE TAUXX, TAUYY, TAUZZ, TAUXY OR TAUYX,
    // TAUYZ OR TAUZY, TAUXZ OR TAUZX
    ViscTens[0] = ViscStress[0][0];
    ViscTens[1] = ViscStress[1][1];
    ViscTens[2] = ViscStress[2][2];
    ViscTens[3] = ViscStress[0][1];
    ViscTens[4] = ViscStress[1][2];
    ViscTens[5] = ViscStress[0][2];

    return ;
}



void BraginskiiCTU::calc_spatial_derivative(MFP* mfp, Vector<UpdateData>& update, const Real time, const Real dt, const Real flux_register_scale)
{
    BL_PROFILE("BraginskiiCTU::calc_spatial_derivative");

    const Geometry& geom = mfp->Geom();
    const int level = mfp->get_level();
    const int finest_level = mfp->get_parent()->finestLevel();

    const size_t n_states = states.size();
    Vector<int> data_indexes;
    Vector<EulerianState*> data_states;

    for (const int& global_idx : state_indexes) {
        EulerianState &istate = EulerianState::get_state_global(global_idx);
        data_indexes.push_back(istate.data_idx);
        data_states.push_back(&istate);
    }


    // ==========================================================================
    // all of level storage

    Vector<MultiFab*> local_old(n_states); // local data that covers the whole level
    for (size_t i=0; i<n_states; ++i) {
        local_old[i] = &update[data_indexes[i]].U;
        update[data_indexes[i]].dU_status = UpdateData::Status::Changed;
    }
    Vector<FabType> active(n_states); // flag for calculation

    Vector<YAFluxRegister*> fr_as_crse(n_states, nullptr);
    Vector<YAFluxRegister*> fr_as_fine(n_states, nullptr);

    // get the maximum stencil size for ghost cells across all components
    int num_grow = 0;
    for (int idx=0; idx<n_states; ++idx) {
        EulerianState &istate = *data_states[idx];
        num_grow = std::max(num_grow, istate.num_grow);
    }

    for (int idx=0; idx<n_states; ++idx) {
        EulerianState &istate = *data_states[idx];
        const int data_idx = data_indexes[idx];


        if (istate.reflux && level < finest_level) {
            MFP& fine_level = mfp->getLevel(level + 1);
            fr_as_crse[idx] = &fine_level.flux_reg[data_idx];
        }

        if (istate.reflux && level > 0) {
            fr_as_fine[idx] = &mfp->flux_reg[data_idx];
        }
    }

    // ==========================================================================
    // per state storage

    Vector<FArrayBox*> conserved(n_states);
    Vector<FArrayBox> primitives(n_states);
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> R_lo(n_states), R_hi(n_states);
    Vector<Array<FArrayBox, AMREX_SPACEDIM>> fluxes(n_states);

    int nu; // per state size

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    // ==========================================================================
    // iterate over all of the FABs within the level performing reconstruction, flux calculation,
    // and updating the cell-centred data according to the resulting divergence


    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();

        // region over which to perform reconstruction for face values
        const Box rbox = amrex::grow(box, 1);

        // ==========================================================================
        // 1. iterate over all states to set-up the data required for flux calculation
        for (int idx=0; idx<n_states; ++idx) {
            EulerianState &istate = *data_states[idx];

            // get a pointer to the conserved quantities
            conserved[idx] = &(*local_old[idx])[mfi];

            active[idx] = FabType::regular;


            // region over which to get cell centered primitives for reconstruction
            const Box pbox = amrex::grow(box, num_grow);

            // ===============================================
            // 1.1 Calculate primitive values within each cell

            FArrayBox& cons = (*local_old[idx])[mfi];
            FArrayBox& prim = primitives[idx];

            istate.calc_primitives(pbox,
                                   cons,
                                   prim,
                                   dx,
                                   time,
                                   prob_lo);

            //            plot_FAB_2d(prim, 0, "prim[0] - "+istate.name, false, true);


            // fill in any cells that need special boundary values
            istate.update_boundary_cells(pbox,
                                         geom,
                                         prim,
                                         time);



            // =======================================
            // 1.2 Calculate reconstructed face values

            // each cell has a hi and lo side in each direction

            // calculate the reconstructed face values
            istate.calc_reconstruction(rbox,
                                       prim,
                                       R_lo[idx],
                                       R_hi[idx]
                                       );

            // ===================================================================
            // update the face values to time t+1/2 based on the local wave speeds

            // TODO: this step currently assumes a single speed for all components and
            // should be updated to calculate the correct characteristic speeds
            // for each component individually
            if (do_CTU) {
                istate.calc_time_averaged_faces(rbox,
                                                prim,
                                                R_lo[idx],
                                                R_hi[idx],
                                                dx,
                                                dt);
            }

        }

        // 3.1 Setup for flux calculation

        // resize the flux arrays before any get used
        for (int idx=0; idx<n_states; ++idx) {

            if (active[idx] == FabType::covered)
                continue;

            nu = local_old[idx]->nComp();
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                // get a node centered box in direction d that encloses the original box
                Box fbox = surroundingNodes(box, d);

                if (do_CTU) {
                    // enlarge the fluxes for face corrections
                    // grow in all directions but that of the flux
                    IntVect expand(1);
                    expand.setVal(d, 0);
                    fbox = amrex::grow(fbox,expand);
                }

                fluxes[idx][d].resize(fbox, nu);
            }
        }

        // ==========================================================================
        // 3.2 Calculate fluxes

        //---
        for (int idx=0; idx<n_states; ++idx) {
            EulerianState &istate = *data_states[idx];

            if (active[idx] == FabType::covered) {
                continue;
            }

            istate.update_face_prim(box,
                                    geom,
                                    R_lo[idx],
                                    R_hi[idx],
                                    time);
        }

        //---
        for (int idx=0; idx<n_states; ++idx) {
            EulerianState &istate = *data_states[idx];


            if (active[idx] == FabType::covered) continue;

            istate.calc_fluxes(box,
                               *conserved[idx],
                               R_lo[idx],
                               R_hi[idx],
                               fluxes[idx],
                               dx, dt);

#if AMREX_SPACEDIM > 1
            if (do_CTU) {
                // now that we have done the fluxes, do we need to correct them???
                // correction is done according to the following steps:
                // 1. calculate the fluxes (above)
                // 2. modify the reconstructed face states using the flux information
                // 3. recalculate the fluxes with the updated face values

                istate.correct_face_prim(box,
                                         R_lo[idx],
                                         R_hi[idx],
                                         fluxes[idx],
                                         dx, dt);

                // following the update of the face values we need to update any boundary conditions

                istate.update_face_prim(box,
                                        geom,
                                        R_lo[idx],
                                        R_hi[idx],
                                        time,
                                        true);
            }
#endif
        }

        //---
#if AMREX_SPACEDIM > 1
        if (do_CTU) {
            for (int idx=0; idx<n_states; ++idx) {
                EulerianState &istate = *data_states[idx];
                if (active[idx] == FabType::covered) continue;

                // recalculate the fluxes
                istate.calc_fluxes(box,
                                   *conserved[idx],
                                   R_lo[idx],
                                   R_hi[idx],
                                   fluxes[idx],
                                   dx, dt);
            }
        }
#endif

        // now calculate any viscous fluxes

        const Box pbox = grow(box, num_grow);

        //        plot_FAB_2d(primitives[+BraginskiiStateIdx::Ion], 0, "ions[0]", false, false);
        //        plot_FAB_2d(primitives[+BraginskiiStateIdx::Electron], 0, "electrons[0]", false, true);

        //        calc_ion_viscous_fluxes(box,
        //                                fluxes[+BraginskiiStateIdx::Ion],
        //                pbox,
        //                primitives[+BraginskiiStateIdx::Ion].const_array(),
        //                primitives[+BraginskiiStateIdx::Electron].const_array(),
        //                primitives[+BraginskiiStateIdx::Field].const_array(),
        //                dx);

        //        calc_electron_viscous_fluxes(box,
        //                                     fluxes[+BraginskiiStateIdx::Electron],
        //                pbox,
        //                primitives[+BraginskiiStateIdx::Ion].const_array(),
        //                primitives[+BraginskiiStateIdx::Electron].const_array(),
        //                primitives[+BraginskiiStateIdx::Field].const_array(),
        //                dx);

        //---
        for (int idx=0; idx<n_states; ++idx) {
            EulerianState &istate = *data_states[idx];
            const int data_idx = data_indexes[idx];

            nu = local_old[idx]->nComp();

            FArrayBox& du = update[data_idx].dU[mfi];


            // calculate divergence

            istate.calc_divergence(box,
                                   fluxes[idx],
                                   du,
                                   dx,
                                   dt);

            if (fr_as_crse[idx]) {
                fr_as_crse[idx]->CrseAdd(mfi, {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])}, dx, flux_register_scale, RunOn::Cpu);
            }

            if (fr_as_fine[idx]) {
                fr_as_fine[idx]->FineAdd(mfi, {AMREX_D_DECL(&fluxes[idx][0], &fluxes[idx][1], &fluxes[idx][2])}, dx, flux_register_scale, RunOn::Cpu);
            }


        }

        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }

}


// calculate slopes and pack them serially into a vector
void BraginskiiCTU::calc_slopes(const Box& box,
                                const FArrayBox& src,
                                FArrayBox& slopes,
                                Reconstruction& reco,
                                const Real *dx) const
{
    BL_PROFILE("BraginskiiCTU::calc_slopes");

    const Box slope_box = grow(box, 1);

    slopes.resize(slope_box, AMREX_SPACEDIM);

    const Dim3 lo = amrex::lbound(slope_box);
    const Dim3 hi = amrex::ubound(slope_box);


    FArrayBox buffer(slope_box, 1); //num_slopes() box with spatial info i,j,k and the slopes
    Array4<Real> const& b4 = buffer.array();


    Array4<const Real> const& h4 = src.array();

    const int nc = electron_state->n_cons();
    Vector<Real> U(nc);

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

                // grab the conserved quantities
                for (int n=0; n<nc; ++n) {
                    U[n] = h4(i,j,k,n);
                }

                // get temperature
                b4(i,j,k) = electron_state->gas->get_temperature_from_cons(U);
            }
        }
    }

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        EulerianState::calc_slope(box, buffer, slopes, dx, 0, d, d, reco);
    }

    return;
}

void BraginskiiCTU::get_alpha_beta_coefficients(Real mass_e, Real T_e, Real charge_e,
                                                Real charge_i, Real nd_e, Real nd_i, Real& alpha_0, Real& alpha_1,
                                                Real& alpha_2,
                                                Real& beta_0, Real& beta_1, Real& beta_2, Real& t_c_e, Real p_lambda,
                                                Real Bx, Real By, Real Bz)
{

    //collision time nondimensional
    Real Debye = MFP::Debye, Larmor = MFP::Larmor, n0_ref=MFP::n0;

    Real pi_num = 3.14159265358979323846;

    t_c_e = std::pow(Debye,4)*n0_ref
            *(6*std::sqrt(2*mass_e)*std::pow(pi_num*T_e, 3./2.)) /
            (p_lambda*std::pow((charge_i/-charge_e),2)*nd_i);


    if (not std::isfinite(t_c_e)) Abort();


    Real omega_ce = -charge_e * std::sqrt( Bx*Bx + By*By + Bz*Bz ) / mass_e / Larmor;

    //TODO create a table lookup for the coefficients based off the plasma constituents
    Real delta_coef, x_coef ;// coefficients used exclusively in the braginskii
    // transport terms
    x_coef = omega_ce*t_c_e;
    Real delta_0 =3.7703, delta_1 = 14.79;
    delta_coef = x_coef*x_coef*x_coef*x_coef+delta_1*x_coef*x_coef + delta_0;// TODO delta0 tables


    Real a_0 = 0.5129, a_0_p =1.837, a_1_p =6.416, a_0_pp =0.7796, a_1_pp =1.704;

    alpha_0 = mass_e*nd_e/t_c_e*a_0;

    if (braginskii_anisotropic) {
        alpha_1 = mass_e*nd_e/t_c_e*( 1 - (a_1_p*x_coef*x_coef + a_0_p)/delta_coef);
        alpha_2 = mass_e*nd_e/t_c_e*x_coef*(a_1_pp*x_coef*x_coef+a_0_pp)/delta_coef;
    } else {
        alpha_1 = 0; alpha_2 = 0;
    }


    Real b_0 = 0.711, b_0_pp = 3.053, b_0_p=2.681, b_1_p=5.101, b_1_pp=3./2.;

    beta_0 = nd_e*b_0;
    if (braginskii_anisotropic) {
        beta_1 = nd_e*(b_1_p*x_coef*x_coef+b_0_p)/delta_coef;
        beta_2 = nd_e*x_coef*(b_1_pp*x_coef*x_coef+b_0_pp)/delta_coef;
    } else {
        beta_1 = 0;
        beta_2 = 0;
    }

    return;
}

int BraginskiiCTU::rhs(Real t, Array<Real,+VectorIdx::NUM> y, Array<Real,+VectorIdx::NUM> &dydt, Array<Real,+DataIdx::NUM>& data)
{
    BL_PROFILE("BraginskiiSource::rhs");

    //TODO find away around the variable declaration in the case of isotropic - may just need to bite the bullet and hav a spearate function
    Array<Real,3> B_unit; //Magnetic field unit vector
    Array<Real,3> u_para; //Velocity parallel to B_unit
    Array<Real,3> u_perp; //Velocity perpendicular to B_unit
    Array<Real,3> u_chev; //unit vector perp to u and B_unit
    Array<Real,3> TG_para;//Temp grad parallel to B_unit
    Array<Real,3> TG_perp;//Temp grad perpendicular to B_unit
    Array<Real,3> TG_chev;//unit vector perp to gradT and B_unit

    Real B_p=0.,B_pp=0.,bx_pp=0.,by_pp=0.,bz_pp=0.,bx_p=0.,by_p=0., xB=0, yB=0, zB=0; //initialised and set to zero to allow get_alpha_beta_coeffs to be run without any changes for both iso and aniso cae


    // magnetic field
    xB = data[+DataIdx::Bx];
    yB = data[+DataIdx::By];
    zB = data[+DataIdx::Bz];

    Real B = xB*xB + yB*yB + zB*zB;

    if (B < effective_zero) {
        B_pp = 0.;
        B_p  = 0.;
    } else if ( (std::abs(xB) < effective_zero) && (std::abs(yB) < effective_zero) && (std::abs(zB) > effective_zero) ) {
        B_pp = 1/sqrt(B); // B prime prime
        B_p  = 0.;
    } else {
        B_pp = 1/sqrt(B); // B prime prime
        B_p  = 1/sqrt(xB*xB + yB*yB); // B prime
    }

    bx_pp = xB*B_pp; bx_p = xB*B_p;
    by_pp = yB*B_pp; by_p = yB*B_p;
    bz_pp = zB*B_pp;

    B_unit[0] = bx_pp; B_unit[1] = by_pp; B_unit[2] = bz_pp;


    const Real m_e = data[+DataIdx::ElectronMass];
    const Real q_e = data[+DataIdx::ElectronCharge];
    const Real gam_e =  data[+DataIdx::ElectronGamma];

    const Real rho_e = data[+DataIdx::ElectronDensity];

    const Real n_e = rho_e/m_e;


    if (rho_e <= effective_zero)
        amrex::Abort("density less than effective zero");

    const Real inv_rho_e = 1/rho_e;
    const Real u_e = y[+VectorIdx::ElectronXmom]*inv_rho_e;
    const Real v_e = y[+VectorIdx::ElectronYmom]*inv_rho_e;
    const Real w_e = y[+VectorIdx::ElectronZmom]*inv_rho_e;
    const Real p_e = (gam_e - 1.0)*(y[+VectorIdx::ElectronEden] - 0.5*rho_e*(u_e*u_e + v_e*v_e + w_e*w_e));
    const Real T_e = p_e*m_e*inv_rho_e;

    const Real dT_dx = y[+DataIdx::dTdx];

#if AMREX_SPACEDIM >= 2
    const Real dT_dy = y[+DataIdx::dTdy];
#else
    const Real dT_dy=0;
#endif

#if AMREX_SPACEDIM == 3
    const Real dT_dz = y[+DataIdx::dTdz];
#else
    const Real dT_dz=0;
#endif


    const Real m_i = data[+DataIdx::IonMass];
    const Real q_i = data[+DataIdx::IonCharge];
    const Real gam_i =  data[+DataIdx::IonGamma];

    const Real rho_i = data[+DataIdx::IonDensity];

    const Real n_i = rho_i/m_i;


    if (rho_i <= effective_zero)
        amrex::Abort("density less than effective zero");

    const Real inv_rho_i = 1/rho_i;
    const Real u_i = y[+VectorIdx::IonXmom]*inv_rho_i;
    const Real v_i = y[+VectorIdx::IonYmom]*inv_rho_i;
    const Real w_i = y[+VectorIdx::IonZmom]*inv_rho_i;
    const Real p_i = (gam_i - 1.0)*(y[+VectorIdx::IonEden] - 0.5*rho_i*(u_i*u_i + v_i*v_i + w_i*w_i));
    const Real T_i = p_i*m_i*inv_rho_i;

    const Real du = u_e - u_i;
    const Real dv = v_e - v_i;
    const Real dw = w_e - w_i;

    //Braginskii directionality stuff.
    if (braginskii_anisotropic) {
        Real dot_B_unit_TG, dot_B_unit_U ;//temp variables

        dot_B_unit_U = bx_pp*du + by_pp*dv + bz_pp*dw;
        dot_B_unit_TG = bx_pp*dT_dx + by_pp*dT_dy + bz_pp*dT_dz;

        for (int i_disp = 0; i_disp < 3; ++i_disp) {//i_disp - i disposable
            u_para[i_disp] = B_unit[i_disp]*dot_B_unit_U ;
            TG_para[i_disp]= B_unit[i_disp]*dot_B_unit_TG ;
            if (i_disp==0){
                u_perp[i_disp] = du - u_para[0]; // couuld be automated with
                // index = prim_vel_id[0] + i_disp
                u_chev[i_disp] = B_unit[1]*dw-B_unit[2]*dv;
                TG_perp[i_disp]= dT_dx - TG_para[0];  //...automated with dT_i vector...
                TG_chev[i_disp]= B_unit[1]*dT_dz-B_unit[2]*dT_dy;
            }
            else if (i_disp==1) {
                u_perp[1] = dv - u_para[1];
                u_chev[1] = -(B_unit[0]*dw-B_unit[2]*du);
                TG_perp[1] = dT_dy - TG_para[1];
                TG_chev[1]= -(B_unit[0]*dT_dz-B_unit[2]*dT_dx);
            }
            else {
                u_perp[i_disp] = dw - u_para[2];
                u_chev[i_disp] = B_unit[0]*dv-B_unit[1]*du;
                TG_perp[i_disp] = dT_dz - TG_para[2];
                TG_chev[i_disp]= B_unit[0]*dT_dy-B_unit[1]*dT_dx;
            }
        }
    }

    //---------------Braginskii Momentum source
    Real alpha_0, alpha_1, alpha_2, beta_0, beta_1, beta_2, t_c_a;
    Real p_lambda = get_coulomb_logarithm(T_i,T_e,n_e);

    get_alpha_beta_coefficients(m_e, T_e, q_e, q_i, n_e, n_i, alpha_0, alpha_1, alpha_2,
                                beta_0, beta_1, beta_2, t_c_a, p_lambda, xB, yB, zB);


    Array<Real,3> R_u, R_T;
    if (braginskii_anisotropic) {
        //frictional force
        R_u[0] = -alpha_0*u_para[0] - alpha_1*u_perp[0] + alpha_2*u_chev[0];
        R_u[1] = -alpha_0*u_para[1] - alpha_1*u_perp[1] + alpha_2*u_chev[1];
        R_u[2] = -alpha_0*u_para[2] - alpha_1*u_perp[2] + alpha_2*u_chev[2];
        //Thermal force

        R_T[0] = -beta_0*TG_para[0] - beta_1*TG_perp[0] - beta_2*TG_chev[0];
        R_T[1] = -beta_0*TG_para[1] - beta_1*TG_perp[1] - beta_2*TG_chev[1];
        R_T[2] = -beta_0*TG_para[2] - beta_1*TG_perp[2] - beta_2*TG_chev[2];
    } else {
        R_u[0] = -alpha_0*du;
        R_u[1] = -alpha_0*dv;
        R_u[2] = -alpha_0*dw;
        //Thermal force

        R_T[0] = -beta_0*dT_dx;
        R_T[1] = -beta_0*dT_dy;
        R_T[2] = -beta_0*dT_dz;
    }
    //Thermal equilibration
    Real Q_delta = 3*m_e/m_i*n_e/t_c_a*(T_e-T_i);
    Real Q_fric  = (R_u[0]+R_T[0])*du + (R_u[1]+R_T[1])*dv + (R_u[2]+R_T[2])*dw;

    const Real d_mx = R_u[0]+R_T[0];
    const Real d_my = R_u[1]+R_T[1];
    const Real d_mz = R_u[2]+R_T[2];
    const Real d_ed = -Q_delta + Q_fric;

    dydt[+VectorIdx::ElectronXmom] = d_mx;
    dydt[+VectorIdx::ElectronYmom] = d_my;
    dydt[+VectorIdx::ElectronZmom] = d_mz;
    dydt[+VectorIdx::ElectronEden] = d_ed;

    dydt[+VectorIdx::IonXmom] = -d_mx;
    dydt[+VectorIdx::IonYmom] = -d_my;
    dydt[+VectorIdx::IonZmom] = -d_mz;
    dydt[+VectorIdx::IonEden] = -d_ed;



    int idx=0;
    for (const auto& val : dydt) {
        if (not std::isfinite(val)) {
            return 1;
        }
        idx++;
    }


    return 0;

}

bool BraginskiiCTU::check_invalid(Array<Real,+VectorIdx::NUM> &y, Array<Real,+DataIdx::NUM>& data)
{
    if (y[+VectorIdx::IonEden] < effective_zero) {
        return true;
    }

    if (y[+VectorIdx::ElectronEden] < effective_zero) {
        return true;
    }

    int i = 0;
    for (const auto& val : y) {
        if (not std::isfinite(val)) {
            return true;
        }
        i++;
    }



    const Real gam_e =  data[+DataIdx::ElectronGamma];
    const Real rho_e = data[+DataIdx::ElectronDensity];
    const Real inv_rho_e = 1/rho_e;
    const Real u_e = y[+VectorIdx::ElectronXmom]*inv_rho_e;
    const Real v_e = y[+VectorIdx::ElectronYmom]*inv_rho_e;
    const Real w_e = y[+VectorIdx::ElectronZmom]*inv_rho_e;
    const Real p_e = (gam_e - 1.0)*(y[+VectorIdx::ElectronEden] - 0.5*rho_e*(u_e*u_e + v_e*v_e + w_e*w_e));

    if (p_e < effective_zero) {
        return true;
    }


    const Real gam_i =  data[+DataIdx::IonGamma];
    const Real rho_i = data[+DataIdx::IonDensity];
    const Real inv_rho_i = 1/rho_i;
    const Real u_i = y[+VectorIdx::IonXmom]*inv_rho_i;
    const Real v_i = y[+VectorIdx::IonYmom]*inv_rho_i;
    const Real w_i = y[+VectorIdx::IonZmom]*inv_rho_i;
    const Real p_i = (gam_i - 1.0)*(y[+VectorIdx::IonEden] - 0.5*rho_i*(u_i*u_i + v_i*v_i + w_i*w_i));

    if (p_i < effective_zero) {
        return true;
    }

    return false;
}

void BraginskiiCTU::calc_time_derivative(MFP* mfp, Vector<UpdateData> &update, const Real time, const Real dt)
{
    BL_PROFILE("Plasma5::explicit_solve");

    const int nc_i = ion_state->n_cons();
    const int nc_e = electron_state->n_cons();

    // collect all of the MultiFabs that we need
    MultiFab& cost = mfp->get_new_data(MFP::Cost_Idx);

    MultiFab& ion_data = update[ion_state->data_idx].U;
    MultiFab& field_data = update[field_state->data_idx].U;
    MultiFab& electron_data = update[electron_state->data_idx].U;

    update[ion_state->data_idx].dU_status = UpdateData::Status::Changed;
    update[electron_state->data_idx].dU_status = UpdateData::Status::Changed;

    Vector<Real> U_i(nc_i), U_e(nc_e);

    const Real* dx = mfp->Geom().CellSize();

    Array<Real,+VectorIdx::NUM> y;
    Array<Real,+DataIdx::NUM> data;


    for (MFIter mfi(cost); mfi.isValid(); ++mfi) {

        Real wt = ParallelDescriptor::second();

        const Box& box = mfi.tilebox();
        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);



        Array4<const Real> const& field4 = field_data.array(mfi);

        Array4<const Real> const& ion4 = ion_data.array(mfi);
        Array4<Real> const& ion_dU4 = update[ion_state->data_idx].dU.array(mfi);

        Array4<const Real> const& electron4 = electron_data.array(mfi);
        Array4<Real> const& electron_dU4 = update[electron_state->data_idx].dU.array(mfi);

        // get the electron temperature slopes
        FArrayBox slopes;

        calc_slopes(box, electron_data[mfi], slopes, *(electron_state->reconstructor.get()), dx);

        Array4<const Real> const& slopes4 = slopes.array();



        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {


                    // load all of the states into vectors
                    for (size_t l=0; l<nc_i; ++l) {
                        U_i[l] = ion4(i,j,k,l);
                    }

                    for (size_t l=0; l<nc_e; ++l) {
                        U_e[l] = electron4(i,j,k,l);
                    }

                    data[+DataIdx::IonMass] = ion_state->gas->get_mass_from_cons(U_i);;
                    data[+DataIdx::IonCharge] = ion_state->gas->get_charge_from_cons(U_i);
                    data[+DataIdx::IonGamma] = ion_state->gas->get_gamma_from_cons(U_i);
                    data[+DataIdx::IonDensity] = U_i[+HydroDef::ConsIdx::Density];
                    data[+DataIdx::ElectronMass] = electron_state->gas->get_mass_from_cons(U_e);
                    data[+DataIdx::ElectronCharge] = electron_state->gas->get_charge_from_cons(U_e);
                    data[+DataIdx::ElectronGamma] = electron_state->gas->get_gamma_from_cons(U_e);
                    data[+DataIdx::ElectronDensity] = U_e[+HydroDef::ConsIdx::Density];
                    data[+DataIdx::Bx] = field4(i,j,k,+FieldDef::ConsIdx::Bx);
                    data[+DataIdx::By] = field4(i,j,k,+FieldDef::ConsIdx::By);
                    data[+DataIdx::Bz] = field4(i,j,k,+FieldDef::ConsIdx::Bz);

                    for (size_t l=0; l<AMREX_SPACEDIM; ++l) {
                        data[+DataIdx::dTdx + l] = slopes4(i,j,k,l);
                    }

                    y[+VectorIdx::IonXmom] = U_i[+HydroDef::ConsIdx::Xmom];
                    y[+VectorIdx::IonYmom] = U_i[+HydroDef::ConsIdx::Ymom];
                    y[+VectorIdx::IonZmom] = U_i[+HydroDef::ConsIdx::Zmom];
                    y[+VectorIdx::IonEden] = U_i[+HydroDef::ConsIdx::Eden];
                    y[+VectorIdx::ElectronXmom] = U_e[+HydroDef::ConsIdx::Xmom];
                    y[+VectorIdx::ElectronYmom] = U_e[+HydroDef::ConsIdx::Ymom];
                    y[+VectorIdx::ElectronZmom] = U_e[+HydroDef::ConsIdx::Zmom];
                    y[+VectorIdx::ElectronEden] = U_e[+HydroDef::ConsIdx::Eden];

                    Real t = time;
                    Real h = dt;
                    int depth = 0;
                    rk4_adaptive(t, y, data, BraginskiiCTU::rhs, BraginskiiCTU::check_invalid, h, time_refinement_factor, depth, max_time_refinement);

                    ion_dU4(i,j,k,+HydroDef::ConsIdx::Xmom) = y[+VectorIdx::IonXmom] - U_i[+HydroDef::ConsIdx::Xmom];
                    ion_dU4(i,j,k,+HydroDef::ConsIdx::Ymom) = y[+VectorIdx::IonYmom] - U_i[+HydroDef::ConsIdx::Ymom];
                    ion_dU4(i,j,k,+HydroDef::ConsIdx::Zmom) = y[+VectorIdx::IonZmom] - U_i[+HydroDef::ConsIdx::Zmom];
                    ion_dU4(i,j,k,+HydroDef::ConsIdx::Eden) = y[+VectorIdx::IonEden] - U_i[+HydroDef::ConsIdx::Eden];
                    electron_dU4(i,j,k,+HydroDef::ConsIdx::Xmom) = y[+VectorIdx::ElectronXmom] - U_e[+HydroDef::ConsIdx::Xmom];
                    electron_dU4(i,j,k,+HydroDef::ConsIdx::Ymom) = y[+VectorIdx::ElectronYmom] - U_e[+HydroDef::ConsIdx::Ymom];
                    electron_dU4(i,j,k,+HydroDef::ConsIdx::Zmom) = y[+VectorIdx::ElectronZmom] - U_e[+HydroDef::ConsIdx::Zmom];
                    electron_dU4(i,j,k,+HydroDef::ConsIdx::Eden) = y[+VectorIdx::ElectronEden] - U_e[+HydroDef::ConsIdx::Eden];
                }
            }
        }

        // update the cost function
        wt = (ParallelDescriptor::second() - wt) / box.d_numPts();
        cost[mfi].plus(wt, box);
    }
}

#endif
