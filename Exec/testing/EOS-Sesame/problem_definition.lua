-- ======== ONE-ZONE SELF-TEST: SESAME-EXTRACTED TABLES ==========
--
-- No time steps (max_step = 0): the eos_table_self_test(...) calls run at
-- config time and print "EOSTAB-SELFTEST[...] PASS/FAIL" lines that
-- check.py parses. DEBUG builds only (see run).
--
-- Every table here is real conditioned SESAME data (tier-2 policy, as for
-- D_fpeos in EOS-Table): the robustness checks (reader, round trips,
-- identities, hull, corner) are gated; 'fd-vs-blocks' is reported but not
-- gated — the cavitated-response floor (doc/eos_sesame_plan.md §3.0)
-- makes the derivative blocks deliberately stiffer than the value surface
-- inside the crossover band.

-- === REFERENCE QUANTITIES (unused by the self-test; needed for a valid config) ===

ref_length = 1.0
ref_density = 1.0
ref_mass = 1.6726219e-27
ref_temp = 1e5

-- === SETTINGS ===

verbosity = 1
cfl = 0.5
time_integration_scheme = 'RK2'

-- === TABLE SELF-TESTS (doc/eos_sesame_plan.md WS5) ===

eos_table_self_test('data/copper_3337_s311.eostab', 48)
eos_table_self_test('data/deuterium_5267_s301.eostab', 48)
eos_table_self_test('data/diamond_7834_s301.eostab', 48)
eos_table_self_test('data/hydrogen_5251_s301.eostab', 48)
eos_table_self_test('data/ti-beta-21s_2963_s311.eostab', 48)
eos_table_self_test('data/ti-beta-21s_2963_trackP.eostab', 48)
eos_table_self_test('data/ti-beta-21s_2963_coldext.eostab', 48)

-- common-energy-reference set (eref295): the three mixture members on one
-- shared energy gauge (e re-referenced to each material's 295 K fill
-- state + ONE shared positivity shift; HANDOFF_common_energy_reference.md)
eos_table_self_test('data/ti-beta-21s_2963_coldext_eref295.eostab', 48)
eos_table_self_test('data/deuterium_5267_s301_eref295.eostab', 48)
eos_table_self_test('data/dry-air_5031_s301_eref295.eostab', 48)
eos_table_self_test('data/aluminum_3720_coldext_eref295.eostab', 48)
eos_table_self_test('data/diamond_7834_s301_eref295.eostab', 48)
-- sub-floor T extension of the highest-floor member (air_lowT_extension/): same gauge, T floor 17.7 K
eos_table_self_test('air_lowT_extension/data/dry-air_5031_s301_eref295_Tf1p25.eostab', 48)

-- === MINIMAL STATE (so the config is valid; never advanced) ===

states = {
  fluid = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = 1.0,
      charge = 0.0,
      gamma = 1.4,
    },
    reconstruction = 'minmod',
    flux = 'HLLE',
    value = {
      rho = 1.0,
      x_vel = 0,
      y_vel = 0,
      z_vel = 0,
      p = 1.0,
    },
  },
}
