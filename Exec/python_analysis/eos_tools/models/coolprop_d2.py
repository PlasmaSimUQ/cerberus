"""Fluid deuterium via CoolProp (Richardson et al. 2014 fundamental EOS).

I. A. Richardson, J. W. Leachman, E. W. Lemmon, "Fundamental Equation of
State for Deuterium", J. Phys. Chem. Ref. Data 43, 013103 (2014) — the
reference Helmholtz-explicit EOS for fluid D2, valid from the melting line
to 600 K at pressures to 2 GPa. CoolProp's "Deuterium" fluid implements it.

CoolProp is a *generation-time* dependency only (guarded import); committed
.eostab artifacts keep the test gates CoolProp-free (plan §1.4).

Outputs are CGS-specific: p [erg/cc], e [erg/g], s [erg/g/K], cv [erg/g/K],
cs [cm/s]. Energy/entropy zero points are CoolProp's own reference state —
arbitrary, handled by the splice energy alignment; only within-model
consistency matters here (it is exact: single Helmholtz functional).
Out-of-range queries return NaN + valid=False rather than extrapolating.
"""

import numpy as np

try:
    import CoolProp.CoolProp as CP
    HAVE_COOLPROP = True
except ImportError:  # pragma: no cover
    CP = None
    HAVE_COOLPROP = False

FLUID = "Deuterium"

# SI -> CGS conversions
PA = 10.0          # Pa -> barye
JKG = 1.0e4        # J/kg -> erg/g
JKGK = 1.0e4       # J/(kg K) -> erg/(g K)
MS = 1.0e2         # m/s -> cm/s
KGM3 = 1.0e-3      # kg/m3 -> g/cc


class CoolPropD2:
    """rho-T evaluations of fluid D2 (not a HelmholtzModel: CoolProp already
    guarantees consistency internally; we only wrap and unit-convert)."""

    name = "coolprop_d2"

    def __init__(self):
        if not HAVE_COOLPROP:
            raise ImportError(
                "CoolProp is required to evaluate the fluid-D2 model "
                "(pip install coolprop; see environment.yml)")
        self.T_min = CP.PropsSI("Tmin", FLUID)
        self.T_max = CP.PropsSI("Tmax", FLUID)
        self.p_max = CP.PropsSI("pmax", FLUID)

    def eval(self, rho, T):
        """Dict of CGS arrays + 'valid' mask; NaN where CoolProp refuses."""
        rho, T = np.broadcast_arrays(np.asarray(rho, float),
                                     np.asarray(T, float))
        shape = rho.shape
        out = {k: np.full(shape, np.nan) for k in
               ("p", "e", "s", "cv", "cs", "a")}
        valid = np.zeros(shape, bool)
        it = np.nditer([rho, T], flags=["multi_index"])
        for r, t in it:
            idx = it.multi_index
            if not (self.T_min <= t <= self.T_max):
                continue
            try:
                rho_si = float(r) / KGM3
                args = ("D", rho_si, "T", float(t), FLUID)
                p = CP.PropsSI("P", *args)
                if not (0.0 < p <= self.p_max):
                    continue
                # reject two-phase evaluations; Q raises or returns -1 for
                # single-phase states depending on backend/region
                try:
                    q = CP.PropsSI("Q", *args)
                except ValueError:
                    q = -1.0
                if 0.0 <= q <= 1.0:
                    continue
                e = CP.PropsSI("UMASS", *args)
                s = CP.PropsSI("SMASS", *args)
                out["p"][idx] = p * PA
                out["e"][idx] = e * JKG
                out["s"][idx] = s * JKGK
                out["cv"][idx] = CP.PropsSI("CVMASS", *args) * JKGK
                out["cs"][idx] = CP.PropsSI("A", *args) * MS
                out["a"][idx] = (e - float(t) * s) * JKG
                valid[idx] = True
            except ValueError:
                continue
        out["valid"] = valid
        out["rho"], out["T"] = rho, T
        return out
