-- ======== SOD SHOCK TUBE: TABULATED-EOS TWIN RUN (Stage 3, W9) ==========
--
-- The classic Sod problem (a diaphragm between a high-pressure and a
-- low-pressure gas bursts at t=0, launching a shock, a contact surface and
-- a rarefaction), run TWICE by the run script with identical numerics:
--
--   GAS_TYPE = 'ideal'      -> thermally_perfect gas, gamma = 1.4
--   GAS_TYPE = 'tabulated'  -> tabulated gas on a synthetic gamma=1.4 table
--
-- Stage 5 (W10) adds a FLUX global: each .inputs file may set the Riemann
-- solver explicitly (FLUX = 'HLLC' or 'HLLC_general_eos'); leaving it unset
-- on the tabulated run exercises the config default (tabulated states with
-- no 'flux' key select 'HLLC_general_eos' — STAGE5.md D-e/G5). The four
-- runs the script drives give the solver-isolation and A/B gates:
--   ideal              TPG + HLLC              (baseline)
--   ideal_geos         TPG + HLLC_general_eos  (G1: round-off vs baseline)
--   tabulated          table + default solver  (G2/G5)
--   tabulated_effgamma table + HLLC            (Stage-3 effective_gamma A/B)
--
-- Because the table IS the ideal gas (generated in closed form by
-- eos_table_prep.py), the two runs must agree to table-interpolation
-- accuracy; any larger difference is an EOS-path bug with everything else
-- identical. check.py also compares both against the exact Riemann solution
-- and checks conservation. GAS_TYPE is set in each .inputs file just before
-- this file is loaded.
--
-- Unit choices (so the code-unit states land inside the table's hull):
-- the table stores CGS with rho in [1e-4, 10] g/cc and T in [1e3, 1e7] K.
-- With ref_density = 1000 kg/m^3 (= 1 g/cc), ref_mass = m_D and
-- ref_temp = 1e5 K, the Sod states rho = {1, 0.125}, p = {1, 0.1} in code
-- units sit at 1 and 0.125 g/cc and T = 1e5 / 0.8e5 K — comfortably inside.

-- === REFERENCE QUANTITIES ===

ref_length = 1.0
ref_density = 1000.0        -- kg/m^3 -> code rho=1 is 1 g/cc
ref_mass = 3.34358e-27      -- kg (deuteron mass, matching the table)
ref_temp = 1e5              -- K

-- === SETTINGS ===

verbosity = 1
cfl = 0.5
time_integration_scheme = 'RK2'

-- === PROBLEM ===

if GAS_TYPE == 'tabulated' then
  gas_def = {
    type = 'tabulated',
    table = 'data/sod_synthetic.eostab',
    mass = 1.0,
    charge = 0.0,
  }
elseif GAS_TYPE == 'ideal' then
  gas_def = {
    type = 'thermally_perfect',
    mass = 1.0,
    charge = 0.0,
    gamma = 1.4,
  }
else
  error("GAS_TYPE must be 'ideal' or 'tabulated' (set in the .inputs file)")
end

function rho0(dat)
  if dat['x'] < 0.0 then
    return 1.0
  else
    return 0.125
  end
end

function p0(dat)
  if dat['x'] < 0.0 then
    return 1.0
  else
    return 0.1
  end
end

states = {
  fluid = {
    type = 'hydro',
    gas = gas_def,
    reconstruction = 'minmod',
    flux = FLUX,  -- nil (unset global) -> key absent -> config default
    value = {
      rho = rho0,
      x_vel = 0,
      y_vel = 0,
      z_vel = 0,
      p = p0,
    },
  },
}

actions = {
  fluxes = {
    type = 'CTU',
    corner_transport = true,
    states = { 'fluid' },
  },
}
