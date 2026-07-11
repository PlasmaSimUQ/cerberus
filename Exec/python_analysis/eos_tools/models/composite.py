"""Helmholtz-exact blending of aligned analytic sources (plan §1.3 step 1).

For weight w(rho, T) on ALIGNED specific Helmholtz energies a1 (w -> 1) and
a2, the blend a = w a1 + (1-w) a2 has exact derivatives

    p = w p1 + (1-w) p2 + rho^2 (dw/drho) (a1 - a2)
    s = w s1 + (1-w) s2 -        (dw/dT)   (a1 - a2)
    e = a + T s = w e1 + (1-w) e2 - T (dw/dT) (a1 - a2)

so consistency is preserved WITHOUT numerical differentiation of the blend.
Alignment (a -> a + cE - T cS) shifts e by cE and s by cS, leaving p and cv
untouched — the constants are fixed in overlap windows and their residual
spreads are QA outputs.

Weights are tanh ramps in one log10 axis with analytic derivatives; a ramp
in log rho has dw/dT = 0 and vice versa, so each blend adds exactly one
correction term.
"""

import numpy as np

LN10 = np.log(10.0)
K_TANH = 1.47  # w: 0.05 -> 0.95 across center +- halfwidth


class LogRamp:
    """w = 0.5 (1 + tanh(K (x - c)/D)) with x = log10(axis value)."""

    def __init__(self, center, halfwidth):
        self.c = center
        self.d = halfwidth

    def w(self, val):
        x = np.log10(np.asarray(val, float))
        return 0.5 * (1.0 + np.tanh(K_TANH * (x - self.c) / self.d))

    def dw_dval(self, val):
        """dw/d(val) — chain rule through x = log10(val)."""
        val = np.asarray(val, float)
        x = np.log10(val)
        sech2 = 1.0 / np.cosh(K_TANH * (x - self.c) / self.d) ** 2
        dw_dx = 0.5 * K_TANH / self.d * sech2
        return dw_dx / (val * LN10)


def blend(out1, out2, w, dw_drho=None, dw_dT=None, rho=None, T=None):
    """Blend two aligned eval-dicts (keys a, p, e, s) with weight w on out1."""
    da = out1["a"] - out2["a"]
    res = {"a": w * out1["a"] + (1.0 - w) * out2["a"],
           "p": w * out1["p"] + (1.0 - w) * out2["p"],
           "e": w * out1["e"] + (1.0 - w) * out2["e"],
           "s": w * out1["s"] + (1.0 - w) * out2["s"]}
    if dw_drho is not None:
        res["p"] = res["p"] + np.asarray(rho, float) ** 2 * dw_drho * da
    if dw_dT is not None:
        res["s"] = res["s"] - dw_dT * da
        res["e"] = res["e"] - np.asarray(T, float) * dw_dT * da
    return res


def align_constants(e_ref, s_ref, e_src, s_src):
    """(cE, cS) making src match ref in an overlap sample; residual stds.

    src' = (e_src + cE, s_src + cS); returns cE, cS, std_e, std_s over the
    sample (the constancy stds are the §1.2-step-2 alignment gate).
    """
    de = np.asarray(e_ref, float) - np.asarray(e_src, float)
    ds = np.asarray(s_ref, float) - np.asarray(s_src, float)
    return float(np.mean(de)), float(np.mean(ds)), \
        float(np.std(de)), float(np.std(ds))


def shifted(out, cE, cS, T):
    """Apply the alignment shift a -> a + cE - T cS to an eval-dict."""
    T = np.asarray(T, float)
    return {"a": out["a"] + cE - T * cS,
            "p": out["p"],
            "e": out["e"] + cE,
            "s": out["s"] + cS}
