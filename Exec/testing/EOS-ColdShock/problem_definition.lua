-- ======== COLD-START SHOCKS ON THE SPLICED WIDE-RANGE TABLE (SS4) ========
--
-- The target is cryogenic liquid deuterium — the physical initial state of
-- shock-driven ICF work (doc/shock-initialization-notes.md: cold start +
-- EOS table is mandatory; ionization is EMERGENT from the jump). A hot
-- driver slab launches a single shock; three driver strengths land the
-- shocked state in the gas-gun, multi-Mbar and classical-plasma regimes.
-- check.py verifies each measured plateau against the Rankine-Hugoniot
-- locus computed from the SAME table and the run's own pre-shock state,
-- and — the new SS4 gate — that the hull-clamp counter stays at ZERO for
-- the physical runs: the entire shock path from 24 K liquid to plasma
-- must live inside the table hull.
--
-- Initial state: rho = 0.171 g/cc (code 1.0), T = 24 K (code 2.4e-4).
-- 24 K rather than ~20 K because the (0.171 g/cc, T) cell is fully
-- in-hull only for T >= ~23.5 K: the two-phase dome edge (satL(20 K) =
-- 0.1718 > 0.171) plus the hull-flag shoulder of the splice's rho
-- hand-off. A slightly-elevated UN-IONIZED start on the cold curve is the
-- legitimate numerical-relief variant (shock-notes section 5); satL(24 K)
-- = 0.164 < 0.171, so this is a single-phase compressed liquid at ~90 bar.
--
-- MODE = 'shock'       (default) driver/target shock tube, needs T_DRIVER
-- MODE = 'rarefaction' abusive expansion of the cold liquid into the
--                      two-phase dome (hull-0): clamps must fire, run
--                      must complete — the robustness gate.
--
-- Table hull in code units: rho in [5.8e-4, 5848], T in [1.78e-4, 1e4].

ref_length = 1.0            -- m
ref_density = 171.0         -- kg/m^3 -> code rho=1 is 0.171 g/cc
ref_mass = 3.34358e-27      -- kg (deuteron, matching the table material)
ref_temp = 1e5              -- K; u_ref = sqrt(kB T_ref / m) ~ 20.3 km/s

verbosity = 2               -- >=2 turns on the hull-clamp per-step report
cfl = 0.5
time_integration_scheme = 'RK2'

T_COLD = 2.4e-4             -- 24 K in code units

MODE = MODE or 'shock'

if MODE == 'shock' then
  if not T_DRIVER then error("MODE='shock' requires T_DRIVER (code units)") end

  -- diaphragm at x = -0.2 (as EOS-Hugoniot): stop times in `run` keep the
  -- shock inside the domain AND the driver rarefaction head off the left
  -- boundary, so the naive-sum conservation gates stay exact
  function T0(dat)
    if dat['x'] < -0.2 then
      return T_DRIVER
    else
      return T_COLD
    end
  end

  init = { rho = 1.0, x_vel = 0, y_vel = 0, z_vel = 0, T = T0 }

elseif MODE == 'rarefaction' then
  -- VACUUM-FORMING expansion (|u| ~ 40x the cold liquid's cs ~ 0.052):
  -- the wide-range table's 17.8 K floor is hydrodynamically UNREACHABLE
  -- (pdV cooling cannot under-run it — measured: a 5x-cs dome expansion
  -- completes with zero clamps), so the only abusive regime left is
  -- density below the table's 1e-4 g/cc edge. The centre cavitates, the
  -- hull rho-clamps and hull-edge floors (W8.1) must carry the run to
  -- completion. This is deliberately the regime the FPEOS-table case
  -- deferred (its absolute floors collapsed dt, W16); the tabulated gas's
  -- hull-physical floors are supposed to close it — this run is the gate.
  function u0(dat)
    if dat['x'] < 0.0 then
      return -2.0
    else
      return 2.0
    end
  end

  init = { rho = 1.0, x_vel = u0, y_vel = 0, z_vel = 0, T = T_COLD }

else
  error("MODE must be 'shock' or 'rarefaction'")
end

states = {
  fluid = {
    type = 'hydro',
    gas = {
      type = 'tabulated',
      table = '../EOS-Table/data/D_spliced.eostab',
      mass = 1.0,
      charge = 0.0,
    },
    reconstruction = 'minmod',
    -- no explicit solver -> the tabulated default HLLC_general_eos
    flux = FLUX,
    value = init,
  },
}

actions = {
  fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = { 'fluid' },
  },
}
