#ifdef EILMER_GAS
    #ifndef EILMERGASMODEL_H
        #define EILMERGASMODEL_H

        #include "MFP_hydro_gas.H"

class EilmerGasModel : public HydroGas
{
  public:
    EilmerGasModel();
    EilmerGasModel(const int global_idx, const sol::table& def);
    ~EilmerGasModel();

    void set_eilmer_gas_state_from_cons(const Vector<Real>& U) const;
    void set_eilmer_gas_state_from_prim(const Vector<Real>& Q) const;

    virtual Real get_gamma_from_cons(const Vector<Real>& U,
                                     const int density_idx = +HydroDef::ConsIdx::Density,
                                     const int tracer_idx = +HydroDef::ConsIdx::NUM) const override;
    virtual Real get_cp_from_cons(const Vector<Real>& U,
                                  const int density_idx = +HydroDef::ConsIdx::Density,
                                  const int tracer_idx = +HydroDef::ConsIdx::NUM) const override;

    virtual Real get_gamma_from_prim(const Vector<Real>& Q,
                                     const int idx = +HydroDef::PrimIdx::NUM) const override;
    virtual Real get_cp_from_prim(const Vector<Real>& Q,
                                  const int tracer_idx = +HydroDef::PrimIdx::NUM) const override;

    virtual bool cons2prim(Vector<Real>& U, Vector<Real>& Q) const override;
    virtual void prim2cons(Vector<Real>& Q, Vector<Real>& U) const override;

    virtual void define_rho_p_T(Vector<Real>& Q) const override;

    virtual Real get_temperature_from_cons(const Vector<Real>& U) const override;
    virtual RealArray get_speed_from_cons(const Vector<Real>& U) const override;
    virtual RealArray get_speed_from_prim(const Vector<Real>& Q) const override;

    virtual std::string get_tag() const override { return tag; }
    virtual void write_info(nlohmann::json& js) const override;

    static bool registered;
    static std::string tag;

    struct SpeciesInfo {
        std::string name;
        int alpha_idx;
    };

    Vector<SpeciesInfo> species_info;
    int n_species;

    int gas_model_id, gas_state_id;
};

    #endif  // EILMERGASMODEL_H
#endif
