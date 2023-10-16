#ifndef CTU_H
#define CTU_H
#ifndef AMREX_USE_EB

#include "MFP_action.H"
#include "MFP_hydro.H"
#include "MFP_field.H"

class BraginskiiCTU : public Action
{
public:
    BraginskiiCTU();
    BraginskiiCTU(const int idx, const sol::table &def);
    ~BraginskiiCTU();

    virtual std::string get_tag() const override {return tag;}
    virtual ActionType get_type() const override {return ActionType::CTU;}

    virtual void get_data(MFP* mfp, Vector<UpdateData>& update, const Real time) const override;

    static Real get_coulomb_logarithm(const Real& T_i, const Real& T_e, const Real& nd_e);

    enum class IonDiffusiveVariable {
        Temp = 0,
        Eta0,
        Eta1,
        Eta2,
        Eta3,
        Eta4,
        Kappa1,
        Kappa2,
        Kappa3,
        NUM_ION_DIFF_COEFFS
    };

    enum class IonDiffusionCoeffs {
        IonTemp = 0,
        IonEta0,
        IonEta1,
        IonEta2,
        IonEta3,
        IonEta4,
        IonKappa1,
        IonKappa2,
        IonKappa3,
        NUM_ION_DIFF_COEFFS
    };



    void get_ion_coeffs(const Vector<Real>& Q_i,
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
                        int& truncatedTau);

    Real get_max_speed_ions(const Vector<amrex::Real>& U_i,
                            const Vector<amrex::Real>& U_e,
                            const Vector<amrex::Real>& U_f);

    void calc_ion_diffusion_terms(const Box& box,
                                  Array4<const Real> const& prim_i4,
                                  Array4<const Real> const& prim_e4,
                                  Array4<const Real> const& prim_f4,
                                  FArrayBox& diff
                                  );

    void calc_ion_viscous_fluxes(const Box& box,
                                 Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                 const Box& pbox,
                                 Array4<const Real> const& prim_i4,
                                 Array4<const Real> const& prim_e4,
                                 Array4<const Real> const& prim_f4,
                                 const Real* dx);

    enum class ElectronDiffusionCoeffs {
        EleTemp = 0,
        EleEta0,
        EleEta1,
        EleEta2,
        EleEta3,
        EleEta4,
        EleKappa1,
        EleKappa2,
        EleKappa3,
        EleBeta1,
        EleBeta2,
        EleBeta3,
        NUM_ELE_DIFF_COEFFS
    };

    enum class EleDiffusiveVariable {
        Temp = 0,
        Eta0,
        Eta1,
        Eta2,
        Eta3,
        Eta4,
        Kappa1,
        Kappa2,
        Kappa3,
        Beta1,
        Beta2,
        Beta3,
        NUM_ELE_DIFF_COEFFS
    };

    void get_electron_coeffs(
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
            int& truncatedTau);

    Real get_max_speed_electrons(const Vector<Real>& U_e,
                                 const Vector<Real>& U_i,
                                 const Vector<Real>& U_f);

    void calc_electron_diffusion_terms(const Box& box,
                                       Array4<const Real> const& prim_i4,
                                       Array4<const Real> const& prim_e4,
                                       Array4<const Real> const& prim_f4,
                                       FArrayBox& diff
                                       );

    void calc_electron_viscous_fluxes(const Box& box,
                                      Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                      const Box& pbox,
                                      Array4<const Real> const& prim_i4,
                                      Array4<const Real> const& prim_e4,
                                      Array4<const Real> const& prim_f4,
                                      const Real* dx);

    enum class FluxSpecies {
        IonFlux=0,
        ElectronFlux
    };

    void calc_charged_viscous_fluxes(FluxSpecies flux_type,
                                     const Box& box,
                                     Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                     Array4<const Real> const& p4,
                                     Array4<const Real> const& prim_i4,
                                     Array4<const Real> const& prim_e4,
                                     Array4<const Real> const& prim_f4,
                                     const Real* dx, FArrayBox& diff);

    void BraginskiiViscousTensorHeatFlux(FluxSpecies flux_type,
                                         Real xB, Real yB, Real zB,
                                         Vector<Real> u_rel,
                                         Real dTdx, Real dTdy, Real dTdz,
                                         Real dudx, Real dudy, Real dudz,
                                         Real dvdx, Real dvdy, Real dvdz,
                                         Real dwdx, Real dwdy, Real dwdz,
                                         Array<Real,+ElectronDiffusionCoeffs::NUM_ELE_DIFF_COEFFS> faceCoefficients,
                                         Array<Real,6>& ViscTens,
                                         Array<Real,3>& q_flux);

    void IsotropicBraginskiiViscousTensorHeatFlux(FluxSpecies flux_type,
                                                  Vector<Real> u_rel,
                                                  Real dTdx, Real dTdy, Real dTdz,
                                                  Real dudx, Real dudy, Real dudz,
                                                  Real dvdx, Real dvdy, Real dvdz,
                                                  Real dwdx, Real dwdy, Real dwdz,
                                                  Array<Real,+ElectronDiffusionCoeffs::NUM_ELE_DIFF_COEFFS> faceCoefficients,
                                                  Array<Real, 6> &ViscTens,
                                                  Array<Real, 3> &q_flux);



    virtual void calc_spatial_derivative(MFP* mfp, Vector<UpdateData>& dU, const Real time, const Real dt, const Real flux_register_scale) override;

    void calc_slopes(const Box& box,
                     const FArrayBox& src,
                     FArrayBox& slopes,
                     Reconstruction& reco,
                     const Real *dx) const;

    static void get_alpha_beta_coefficients(const Real& Z_i, Real mass_e, Real T_e, Real charge_e,
                                     Real charge_i, Real nd_e, Real nd_i, Real& alpha_0, Real& alpha_1,
                                     Real& alpha_2,
                                     Real& beta_0, Real& beta_1, Real& beta_2, Real& t_c_e, Real p_lambda,
                                     Real Bx, Real By, Real Bz);


    enum class DataIdx {
        IonMass=0,
        IonCharge,
        IonGamma,
        IonDensity,
        ElectronMass,
        ElectronCharge,
        ElectronGamma,
        ElectronDensity,
        Bx,
        By,
        Bz,
        AMREX_D_DECL(dTdx,dTdy,dTdz),
        NUM
    };

    enum class VectorIdx {
        IonXmom=0,
        IonYmom,
        IonZmom,
        IonEden,
        ElectronXmom,
        ElectronYmom,
        ElectronZmom,
        ElectronEden,
        NUM
    };

    static int rhs(Real t, Array<Real,+VectorIdx::NUM> y, Array<Real,+VectorIdx::NUM> &dydt, Array<Real,+DataIdx::NUM>& data);
    static bool check_invalid(Array<Real,+VectorIdx::NUM> &y, Array<Real,+DataIdx::NUM>& data);
    virtual void calc_time_derivative(MFP* mfp, Vector<UpdateData>& dU, const Real time, const Real dt) override;

    enum class BraginskiiStateIdx {
        Ion=0,
        Electron,
        Field
    };

    Array<EulerianState*,3> states;

    HydroState *ion_state, *electron_state;
    FieldState* field_state;

    static constexpr Real effective_zero = 1e-14;
    static bool srin_switch;
    static bool braginskii_anisotropic;
    Real cfl = 1.0;

    struct IonCoeffs {
        Real mu_0, T0, BT, forceViscosityValue;
        bool forceViscosity;
        Real T_i,eta0,eta1,eta2,eta3,eta4,kappa1,kappa2,kappa3;
    };

    IonCoeffs ion_coeffs;

    struct ElectronCoeffs {
        Real mu_0, T0, BT, forceViscosityValue;
        bool forceViscosity;
        Real T_e,eta0,eta1,eta2,eta3,eta4,kappa1,kappa2,kappa3,beta1,beta2,beta3;
    };

    ElectronCoeffs electron_coeffs;

    bool do_CTU;

    int time_refinement_factor;
    int max_time_refinement;

    static bool registered;
    static std::string tag;
};

#endif
#endif // CTU_H
