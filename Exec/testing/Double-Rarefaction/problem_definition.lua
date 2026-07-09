-- ======== 3D DOUBLE-RAREFACTION (EINFELDT) + TRACER STEP ==========
--
-- Two symmetric rarefaction (expansion) waves pull the gas apart from x = 0,
-- driving density and pressure toward vacuum at the centre of the domain.
-- This is the standard positivity stress test for Godunov-type schemes:
--
--   [1] B. Einfeldt, C.D. Munz, P.L. Roe, B. Sjogreen,
--       "On Godunov-type methods near low densities",
--       J. Comput. Phys. 92 (1991) 273-295.
--   [2] E.F. Toro, "Riemann Solvers and Numerical Methods for Fluid
--       Dynamics", 3rd ed., Springer (2009), Test 2 of section 4.3.3.
--   [3] X. Zhang, C.-W. Shu, "On positivity-preserving high order
--       discontinuous Galerkin schemes for compressible Euler equations on
--       rectangular meshes", J. Comput. Phys. 229 (2010) 8918-8934.
--       (context: positive face values do not by themselves guarantee
--       positive cell averages after the update)
--
-- Initial conditions: rho = 1, p = 0.4, gamma = 1.4 everywhere,
-- u_x = -u0 for x < 0 and +u0 for x > 0. A passive tracer mass fraction
-- alpha is stepped 0 -> 1 at x = 0.25 so that an unlimited reconstruction
-- scheme is guaranteed to undershoot below zero on the low side of the jump.
-- The domain is fully periodic; the wrap-around at |x| = 0.5 forms a
-- colliding-stream compression which is benign for positivity.
--
-- PURPOSE: regression test for the per-state 'reconstruction_fallback'
-- option. The primary scheme here is the UNLIMITED sixth-order 'O6', which
-- undershoots at the tracer step every step; the 'minmod' fallback redoes
-- exactly those face values. Guarded components (see
-- EulerianState::get_positive_prim) are rho and p (floor = effective_zero)
-- and each tracer alpha (floor = 0, since alpha == 0 is a legitimate value).
--
-- Expected result at u0 = 1, t = 0.2, 160x16x16 (undershoot magnitude is
-- resolution dependent): with the fallback the final min(alpha) is ~ -7e-5
-- (small negative residual because face positivity does not strictly bound
-- the updated cell average, see [3]); without 'reconstruction_fallback' the
-- same run gives min(alpha) ~ -1.2e-3, and coarser grids give larger
-- undershoot (~ -7e-2 at 64x8x8). check.py asserts min(alpha) > -3e-4,
-- rho, p > 0, and that the fallback actually fired.
--
-- Escalation notes (documented behaviour, not exercised by this test):
--   u0 = 2: O6 face values reach negative p -> HLLC aborts (SIGFPE with
--           amrex.fpe_trap_* on) unless the fallback is enabled; with the
--           fallback the flux path survives but cell centres can still be
--           evacuated by the update.
--   u0 = 5: true vacuum forms at cell centres; the run dies in the
--           timestep calculation (get_speed_from_cons) regardless of
--           reconstruction, which requires the planned density/pressure
--           floor layers (apply_prim_floor) rather than this feature.

-- === REFERENCE QUANTITIES ===

ref_length = 1.0
ref_density = 1.0
ref_mass = 1.6726219000e-27
ref_temp = 1e5

-- === SETTINGS ===

verbosity = 2 -- >= 2 so the fallback reports its per-box hit counts
cfl = 0.5
time_integration_scheme = 'RK2'

-- === PROBLEM ===

u0 = 1.0 -- rarefaction strength (see escalation notes above)

function ux(dat)
  if dat['x'] < 0.0 then
    return -u0
  else
    return u0
  end
end

-- tracer step: unlimited schemes undershoot below zero on the low side
function ax(dat)
  if dat['x'] > 0.25 then
    return 1.0
  else
    return 0.0
  end
end

states = {
  fluid = {
    type = 'hydro',
    gas = {
      type = 'thermally_perfect',
      mass = 1.0,
      charge = 0.0,
      gamma = 1.4,
    },
    reconstruction = 'O6', -- unlimited: guaranteed to overshoot at the step
    reconstruction_fallback = 'minmod', -- feature under test
    -- positivity floor threshold: used by the fallback guard for rho and p
    -- (tracers are guarded at exactly 0), by the cons2prim cell-centre
    -- floors, the wave-speed floors and the pre-Riemann face clamp (the
    -- floors compile in only with USE_PRIM_FLOOR=TRUE, the default - see
    -- the run script). 1e-14 is the default; stated explicitly here.
    effective_zero = 1e-14,
    flux = 'HLLC',
    value = {
      rho = 1.0,
      x_vel = ux,
      y_vel = 0,
      z_vel = 0,
      p = 0.4,
      alpha = ax,
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
