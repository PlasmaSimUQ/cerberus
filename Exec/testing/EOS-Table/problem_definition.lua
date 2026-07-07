-- ======== ONE-ZONE TABULATED-EOS SELF-TEST (Stage 2, W5) ==========
--
-- This case runs NO time steps (max_step = 0 in onezone.inputs). Its whole
-- purpose is the two eos_table_self_test(...) calls below, which execute at
-- config time (this script runs while the code reads its configuration) and
-- print one "EOSTAB-SELFTEST[...] PASS/FAIL" line per check; check.py
-- parses those lines from run_log.txt.
--
-- The self-test function only exists in DEBUG builds (the run script builds
-- with DEBUG=TRUE) and exercises the table in DIMENSIONAL mode — reference
-- quantities play no role here (see STAGE2.md, review amendment #3).
--
-- Tier 1 (the gate): synthetic ideal-gas table — closed-form truth, so the
--   thresholds inside the self-test are meaningful absolutely.
-- Tier 2 (informative): the conditioned FPEOS deuterium table from Stage 1 —
--   consistency checks on real warm-dense-matter data (ragged hull, 76.9%
--   coverage). check.py gates tier 2 only on the robustness checks; see
--   there for what is allowed to differ.

-- === REFERENCE QUANTITIES (unused by the self-test; needed for a valid config) ===

ref_length = 1.0
ref_density = 1.0
ref_mass = 1.6726219e-27
ref_temp = 1e5

-- === SETTINGS ===

verbosity = 1
cfl = 0.5
time_integration_scheme = 'RK2'

-- === TABLE SELF-TESTS ===

eos_table_self_test('data/ideal_synthetic.eostab', 48)
eos_table_self_test('data/D_fpeos.eostab', 48)

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
