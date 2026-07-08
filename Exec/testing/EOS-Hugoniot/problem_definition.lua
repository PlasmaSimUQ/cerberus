-- ======== STRONG-SHOCK HUGONIOT ON THE REAL FPEOS TABLE (Stage 4, W9) ====
--
-- A hot deuterium driver slab launches a strong shock into colder deuterium
-- near the table's temperature floor. Each driver temperature T_DRIVER
-- produces one post-shock state = one point on the Hugoniot (the locus of
-- states reachable by a single shock from a given initial state); check.py
-- verifies the measured points lie on the locus predicted by the SAME
-- table via the Rankine-Hugoniot jump conditions (conservation of mass,
-- momentum and energy across the front).
--
-- The published PIMC/FPEOS Hugoniot starts from cryogenic liquid D
-- (0.171 g/cc, 19.6 K) — below the table floor (15625 K), so the
-- *simulated* locus is centred at our pre-shock state (0.171 g/cc,
-- 17000 K) instead; the offline Stage-1 QA artifact anchors the table to
-- the published locus and check.py re-asserts that agreement separately.
--
-- MODE = 'shock'       (default) driver/target shock tube, needs T_DRIVER
-- MODE = 'rarefaction' abusive double rarefaction (diverging velocities
--                      drive the centre's internal energy below the table
--                      hull; the W8 clamps must fire and the run must
--                      complete — the Stage-4 robustness gate)
--
-- The table hull (code units below): rho in [0.0115, 9333], T in
-- [0.15625, 640]. All initial states sit inside it.

ref_length = 1.0            -- m
ref_density = 171.0         -- kg/m^3 -> code rho=1 is 0.171 g/cc
ref_mass = 3.34358e-27      -- kg (deuteron, matching the table material)
ref_temp = 1e5              -- K; u_ref = sqrt(kB T_ref / m) ~ 20.3 km/s

verbosity = 2               -- >=2 turns on the hull-clamp per-step report
cfl = 0.5
time_integration_scheme = 'RK2'

MODE = MODE or 'shock'

if MODE == 'shock' then
  if not T_DRIVER then error("MODE='shock' requires T_DRIVER (code units)") end

  -- diaphragm at x = -0.2: enough target run for a wide plateau while the
  -- driver rarefaction stays inside the domain until stop_time (so the
  -- naive-sum conservation gate remains exact — nothing leaves)
  function T0(dat)
    if dat['x'] < -0.2 then
      return T_DRIVER
    else
      return 0.17          -- 17000 K, just above the 15625 K table floor
    end
  end

  init = { rho = 1.0, x_vel = 0, y_vel = 0, z_vel = 0, T = T0 }

elseif MODE == 'rarefaction' then
  -- diverging flow at ~3x the local sound speed (cs(1, 0.3) ~ 0.66): the
  -- centre expands and cools through the table's 15625 K floor -> the W8
  -- hull clamps must fire while the run completes. Deliberately NOT a
  -- vacuum-forming |u| (that regime needs the deferred W16 per-variable
  -- floors / vacuum-velocity treatment; the FPEOS hull bottom is reached
  -- at ~3x rarefaction, long before vacuum, so this is the physical gate).
  function u0(dat)
    if dat['x'] < 0.0 then
      return -2.0
    else
      return 2.0
    end
  end

  -- T0 only ~15% above the 15625 K table floor: the very first expansion
  -- cooling pushes the centre below the hull -> clamps fire immediately
  init = { rho = 1.0, x_vel = u0, y_vel = 0, z_vel = 0, T = 0.18 }

else
  error("MODE must be 'shock' or 'rarefaction'")
end

states = {
  fluid = {
    type = 'hydro',
    gas = {
      type = 'tabulated',
      table = '../EOS-Table/data/D_fpeos.eostab',
      mass = 1.0,
      charge = 0.0,
    },
    reconstruction = 'minmod',
    flux = 'HLLC',
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
