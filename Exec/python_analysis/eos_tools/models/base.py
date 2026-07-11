"""Base for the Helmholtz-level analytic EOS models (plan SS2).

Every model defines specific Helmholtz free energy a(rho, T) [erg/g] and
inherits thermodynamically CONSISTENT outputs by differentiating it:

    s  = -da/dT|rho          [erg/g/K]
    e  = a + T s             [erg/g]
    p  = rho^2 da/drho|T     [erg/cc]
    cv = -T d2a/dT2|rho      [erg/g/K]

Derivatives are central differences with relative steps — consistency by
construction is the point of the composite cold model (plan §1.3), and the
per-model closed forms in the unit tests keep the differencing honest.
Models may override any output with an analytic form (same contract).
"""

import numpy as np

# relative finite-difference steps (central, second order). 1e-5 keeps the
# truncation error ~1e-10 relative while staying far above double-precision
# cancellation noise on O(1e10) erg/g free energies.
H_REL = 1.0e-5


class HelmholtzModel:
    """Derive p, e, s, cv from a(rho, T) by central differences."""

    name = "base"

    def a(self, rho, T):  # specific Helmholtz, erg/g
        raise NotImplementedError

    def s(self, rho, T):
        rho, T = np.asarray(rho, float), np.asarray(T, float)
        hT = H_REL * T
        return -(self.a(rho, T + hT) - self.a(rho, T - hT)) / (2.0 * hT)

    def e(self, rho, T):
        return self.a(rho, T) + np.asarray(T, float) * self.s(rho, T)

    def p(self, rho, T):
        rho, T = np.asarray(rho, float), np.asarray(T, float)
        hr = H_REL * rho
        return rho ** 2 * (self.a(rho + hr, T) - self.a(rho - hr, T)) / (2.0 * hr)

    def cv(self, rho, T):
        rho, T = np.asarray(rho, float), np.asarray(T, float)
        hT = H_REL * T
        return -T * (self.a(rho, T + hT) - 2.0 * self.a(rho, T)
                     + self.a(rho, T - hT)) / hT ** 2

    def eval(self, rho, T):
        """Bundle everything (dict of arrays broadcast over rho, T)."""
        rho, T = np.broadcast_arrays(np.asarray(rho, float),
                                     np.asarray(T, float))
        return {
            "rho": rho, "T": T,
            "a": self.a(rho, T), "p": self.p(rho, T),
            "e": self.e(rho, T), "s": self.s(rho, T),
            "cv": self.cv(rho, T),
        }


class SumModel(HelmholtzModel):
    """Helmholtz-additive combination (e.g. cold curve + ion thermal)."""

    def __init__(self, parts, name="sum"):
        self.parts = list(parts)
        self.name = name

    def a(self, rho, T):
        out = self.parts[0].a(rho, T)
        for m in self.parts[1:]:
            out = out + m.a(rho, T)
        return out
