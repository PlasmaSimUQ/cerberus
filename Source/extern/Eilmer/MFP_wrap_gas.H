#ifdef EILMER_GAS
extern "C" int cwrap_gas_init();
extern "C" int gas_model_new(const char* file_name);
extern "C" int gas_model_n_species(int gm_i);
extern "C" int gas_model_species_name_and_length(int gm_i, int isp, char* dest_name, int* n);
extern "C" int gas_state_new(int gm_i);
extern "C" int thermochemical_reactor_new(int gm_i, const char* filename1, const char* filename2);

extern "C" int gas_state_set_scalar_field(int gs_i, const char* field_name, double value);
extern "C" int gas_state_get_scalar_field(int gs_i, const char* field_name, double* value);
extern "C" int gas_state_set_array_field(int gs_i, const char* field_name, double* values, int n);
extern "C" int gas_state_get_array_field(int gs_i, const char* field_name, double* values, int n);
extern "C" int gas_model_gas_state_update_thermo_from_rhou(int gm_i, int gs_i);
extern "C" int gas_model_gas_state_update_thermo_from_rhop(int gm_i, int gs_i);
extern "C" int gas_model_gas_state_update_thermo_from_rhoT(int gm_i, int gs_i);
extern "C" int gas_model_gas_state_update_thermo_from_pT(int gm_i, int gs_i);
extern "C" int thermochemical_reactor_gas_state_update(int cr_i, int gs_i, double t_interval, double* dt_suggest);

extern "C" int gas_model_mol_masses(int gm_i, double* mm);
extern "C" int gas_model_gas_state_gamma(int gm_i, int gs_i, double* result);
extern "C" int gas_model_gas_state_Cp(int gm_i, int gs_i, double* result);

extern "C" int gas_model_gas_state_update_sound_speed(int gm_i, int gs_i);

class EilmerGas {
public:
    static bool initialised;
    static void initialise();
};

#endif
