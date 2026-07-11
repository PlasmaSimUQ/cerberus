"""Deuterium cold/low-T composite model (plan SS2, §1.3 ingredient 1).

Three aligned Helmholtz-level sources, blended with exact derivative
formulas (models/composite.py):

  fluid side = LogRamp in T (CoolProp/Richardson D2 below ~550 K,
               ideal D2 rotor gas above)                       [T-only ramp]
  total      = LogRamp in rho (fluid side at low rho, QEOS solid —
               LLNL Vinet cold curve + Slater-Debye ion thermal —
               at high rho)                                    [rho-only ramp]

The Slater-Debye temperature comes from the fitted cold curve itself,
theta(rho) = (hbar/kB) (6 pi^2 n)^(1/3) sqrt(B_cold(rho)/rho) — no
hand-entered Debye constants (evaluates to ~119 K at rho0, consistent with
solid-D2 literature).

Alignment: CoolProp -> rotor gas in the dilute-gas window (both molecular
D2; fixes CoolProp's arbitrary reference state), then QEOS -> aligned fluid
in a compressed-liquid window near melt. The second shift absorbs the melt
entropy — the v1 composite deliberately smears melting (no strength, no
latent heat; documented). Residual stds of both alignments are QA gates.

Two-phase dome and any CoolProp-refused cell inside its nominal weight
region fall back to the rotor-gas values and are marked hull = 0.
"""

import os

import numpy as np

from ..constants import GPA_CGS, KB, M_D, NA
from ..models.coldcurve import VinetColdCurve, fit_vinet, vinet_p
from ..models.composite import LogRamp, align_constants, blend, shifted
from ..models.qeos import SlaterDebye
from ..models.rotor_gas import RotorGasD2
from .. import sources as _sources

# reference state (cryogenic liquid D2; matches EOS-Hugoniot ref_density)
RHO0 = 0.171   # g/cc
T0 = 20.0      # K

# blend placement (log10): fluid->gas ramp ends at CoolProp's 600 K ceiling;
# fluid->solid ramp centred where the cold curve reaches ~1.2 GPa
GAS_RAMP = LogRamp(center=2.74, halfwidth=0.035)    # 505..600 K band
SOLID_RAMP_P_GPA = 1.2                              # sets the rho center
SOLID_RAMP_HALFWIDTH = 0.06                         # dex in rho

# alignment windows
W_GAS = dict(rho=np.geomspace(1e-3, 1e-2, 6), T=np.linspace(420.0, 590.0, 6))
W_MELT = dict(rho=np.geomspace(0.24, 0.31, 6), T=np.linspace(25.0, 50.0, 6))


def coldcurve_csv_path():
    return os.path.join(_sources.repo_root(), "Exec", "testing", "EOS-Table",
                        "data", "raw", "llnl_coldcurve", "h_coldcurve_0K.csv")


class DeuteriumColdModel:
    """eval(rho, T) -> dict(a, p, e, s, hull, src) on the aligned scale."""

    name = "cold_model_d2"

    def __init__(self, csv_path=None, use_coolprop=True):
        csv_path = csv_path or coldcurve_csv_path()
        P_GPa, V_mol = np.loadtxt(csv_path, delimiter=",", skiprows=4,
                                  unpack=True)
        v = V_mol / (M_D * NA)  # cc/g, deuterium mass basis
        (v0, B0, B0p), rms = fit_vinet(P_GPa * GPA_CGS, v)
        self.vinet = (v0, B0, B0p)
        self.vinet_rms = rms
        cold = VinetColdCurve(v0, B0, B0p)
        self.solid_parts = (cold, SlaterDebye(M_D, cold))
        self.gas = RotorGasD2()

        # rho-ramp center from the cold curve: rho at SOLID_RAMP_P_GPA
        from scipy.optimize import brentq
        rho_sw = brentq(lambda r: vinet_p(1.0 / r, v0, B0, B0p)
                        - SOLID_RAMP_P_GPA * GPA_CGS, 0.19, 2.0)
        self.solid_ramp = LogRamp(np.log10(rho_sw), SOLID_RAMP_HALFWIDTH)
        self.rho_switch = rho_sw

        self.cp = None
        if use_coolprop:
            from ..models.coolprop_d2 import CoolPropD2
            self.cp = CoolPropD2()

        self._align()

    # --- source eval helpers -------------------------------------------
    def _eval_solid(self, rho, T):
        cold, ion = self.solid_parts
        return {"a": cold.a(rho, T) + ion.a(rho, T),
                "p": cold.p(rho, T) + ion.p(rho, T),
                "e": cold.a(rho, T) * np.ones_like(np.asarray(T, float))
                     + ion.e(rho, T),
                "s": ion.s(rho, T)}

    def _eval_gas(self, rho, T):
        g = self.gas
        return {"a": g.a(rho, T), "p": g.p(rho, T),
                "e": g.e(rho, T), "s": g.s(rho, T)}

    def _satl_fill(self, T):
        """Saturated-liquid props at T (CGS eval-dict), NaN above Tc.

        Fill for CoolProp-refused cells below the critical temperature:
        the physical liquid-branch continuation for the two-phase dome
        (rarefactions into the dome, and the Hugoniot reference state at
        rho0 = 0.171 g/cc which sits marginally inside it)."""
        import CoolProp.CoolProp as CP
        T = np.asarray(T, float)
        Tc = CP.PropsSI("Tcrit", "Deuterium")
        Tmin = self.cp.T_min

        def one(tt):
            if tt >= 0.999 * Tc:
                return (np.nan,) * 4
            # below the EOS floor (18.72 K) evaluate the saturated liquid
            # AT the floor: a constant extension over the table's bottom
            # rows (17.8-18.7 K), far better than the rotor-gas fallback
            # (which puts ~60 bar of ideal-gas pressure at liquid density
            # and then leaks through the dome via the monotonise cummax)
            tt = max(float(tt), Tmin)
            try:
                args = ("T", tt, "Q", 0.0, "Deuterium")
                p = CP.PropsSI("P", *args) * 10.0
                e = CP.PropsSI("UMASS", *args) * 1e4
                s = CP.PropsSI("SMASS", *args) * 1e4
                return p, e, s, e - float(tt) * s
            except ValueError:
                return (np.nan,) * 4

        p, e, s, a = np.vectorize(one, otypes=[float] * 4)(T)
        return {"p": p, "e": e, "s": s, "a": a}

    def _eval_fluid(self, rho, T):
        """CoolProp with saturated-liquid (below Tc) then rotor-gas
        fallback; returns (dict, valid_mask)."""
        gasd = self._eval_gas(rho, T)
        if self.cp is None:
            return gasd, np.zeros(np.shape(gasd["p"]), bool)
        out = self.cp.eval(rho, T)
        valid = out["valid"]
        # fallback ladder: saturated liquid where defined, else rotor gas
        satl = self._satl_fill(T)
        fall = {k: np.where(np.isfinite(satl[k]), satl[k], gasd[k])
                for k in ("a", "p", "e", "s")}
        fall_is_cp = np.isfinite(satl["p"])  # satL is CoolProp-scale
        cpd = {}
        for k in ("a", "p", "e", "s"):
            cpd[k] = np.where(valid, out[k], fall[k])
        # CoolProp-scale values (real or satL fill) get the alignment shift;
        # rotor-gas fill is already on the reference scale
        cp_scale = valid | fall_is_cp
        sh = shifted(cpd, self.cE_cp, self.cS_cp, T)
        for k in ("a", "p", "e", "s"):
            cpd[k] = np.where(cp_scale, sh[k], cpd[k])
        return cpd, valid

    # --- alignment ------------------------------------------------------
    def _align(self):
        self.cE_cp = self.cS_cp = 0.0
        self.align_stats = {}
        rg, Tg = np.meshgrid(W_GAS["rho"], W_GAS["T"], indexing="ij")
        if self.cp is not None:
            cpo = self.cp.eval(rg, Tg)
            m = cpo["valid"]
            gaso = self._eval_gas(rg, Tg)
            cE, cS, se, ss = align_constants(gaso["e"][m], gaso["s"][m],
                                             cpo["e"][m], cpo["s"][m])
            self.cE_cp, self.cS_cp = cE, cS
            # constancy stds normalised by local kT-scale energy/entropy
            kT_e = KB / (2.0 * M_D) * Tg[m]
            self.align_stats["cp_vs_gas"] = dict(
                cE=cE, cS=cS, std_e_over_kT=float(np.mean(se / kT_e)),
                std_s_over_R=float(ss / (KB / (2.0 * M_D))))
        rm, Tm = np.meshgrid(W_MELT["rho"], W_MELT["T"], indexing="ij")
        fl, valid = self._eval_fluid(rm, Tm)
        sol = self._eval_solid(rm, Tm)
        sel = valid if valid.any() else np.ones_like(Tm, bool)
        cE, cS, se, ss = align_constants(fl["e"][sel], fl["s"][sel],
                                         sol["e"][sel], sol["s"][sel])
        self.cE_sol, self.cS_sol = cE, cS
        kT_e = KB / (2.0 * M_D) * Tm[sel]
        self.align_stats["solid_vs_fluid"] = dict(
            cE=cE, cS=cS, std_e_over_kT=float(np.mean(se / kT_e)),
            std_s_over_R=float(ss / (KB / (2.0 * M_D))))

    # --- public eval -----------------------------------------------------
    def eval(self, rho, T):
        rho, T = np.broadcast_arrays(np.asarray(rho, float),
                                     np.asarray(T, float))
        gasd = self._eval_gas(rho, T)  # already the alignment reference
        fluid, cp_valid = self._eval_fluid(rho, T)

        # fluid side: CoolProp (cold) <-> rotor gas (hot), T-only ramp
        w_gas = GAS_RAMP.w(T)
        dwg_dT = GAS_RAMP.dw_dval(T)
        fl = blend(gasd, fluid, w_gas, dw_dT=dwg_dT, T=T)

        # total: solid (high rho) <-> fluid side, rho-only ramp
        sol = shifted(self._eval_solid(rho, T), self.cE_sol, self.cS_sol, T)
        w_sol = self.solid_ramp.w(rho)
        dws_drho = self.solid_ramp.dw_dval(rho)
        tot = blend(sol, fl, w_sol, dw_drho=dws_drho, rho=rho)

        # hull: 0 where the fluid side leaned on a CoolProp fallback
        # (two-phase dome / out-of-range) while carrying real weight
        fluid_weight = (1.0 - w_sol) * (1.0 - w_gas)
        hull = ~((fluid_weight > 0.05) & ~cp_valid)
        src = np.where(w_sol > 0.5, 2, np.where(w_gas > 0.5, 1, 0))
        tot["hull"] = hull.astype(float)
        tot["src"] = src  # 0 = coolprop, 1 = rotor gas, 2 = qeos solid
        tot["rho"], tot["T"] = rho, T
        return tot
