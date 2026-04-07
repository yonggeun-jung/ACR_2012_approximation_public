"""
distributions.py — Productivity/technology distribution classes.

DGPs:
  1. FrechetPareto  — Exact benchmark (ε constant)
  2. LogNormal      — Slow tail convergence (2RV β=0)
  3. ParetoMixture  — Contaminated Pareto (density bunching)
  4. TruncatedPareto — Upper-truncated Pareto (2RV β<0)

Each class provides:
  - pdf, cdf, sf, rvs
  - Phi(phi_star, sigma): revenue-weighted aggregate
  - eps_melitz(phi_star, sigma): Melitz trade elasticity
  - hazard_rate(z): scaled hazard z f(z)/(1-F(z))
  - eps_ek(z): EK trade elasticity (= hazard rate)

Reference:
  Jung (2026), "When Are Sufficient Statistics Sufficient?"
"""

import numpy as np
from scipy import integrate, special
from abc import ABC, abstractmethod


class BaseDistribution(ABC):

    @abstractmethod
    def pdf(self, x): pass

    @abstractmethod
    def cdf(self, x): pass

    def sf(self, x):
        return 1.0 - self.cdf(x)

    @abstractmethod
    def rvs(self, size, rng=None): pass

    # --- EK class ---
    def hazard_rate(self, z):
        """z f(z) / (1 - F(z))."""
        z = np.asarray(z, dtype=float)
        f = self.pdf(z)
        s = self.sf(z)
        return np.where(s > 0, z * f / s, np.inf)

    def eps_ek(self, z):
        """EK trade elasticity = hazard rate."""
        return self.hazard_rate(z)

    # --- Melitz class ---
    def Phi(self, phi_star, sigma):
        """Φ(φ*) = ∫_{φ*}^∞ (φ/φ*)^{σ-1} g(φ) dφ.  Numerical."""
        s = sigma - 1.0
        def integrand(phi):
            return (phi / phi_star) ** s * self.pdf(phi)
        result, _ = integrate.quad(integrand, phi_star, np.inf,
                                   limit=200, epsrel=1e-10)
        return result

    def eps_melitz(self, phi_star, sigma):
        """ε(φ*) = (σ-1) + φ* g(φ*) / Φ(φ*)."""
        s = sigma - 1.0
        Phi_val = self.Phi(phi_star, sigma)
        if Phi_val <= 0:
            return np.inf
        return s + phi_star * self.pdf(phi_star) / Phi_val


# -----------------------------------------------------------------------
# DGP 1: Fréchet / Pareto
# -----------------------------------------------------------------------

class FrechetPareto(BaseDistribution):
    """
    Pareto: G(φ) = 1 - (φ_min/φ)^κ,  φ ≥ φ_min.
    Fréchet interpretation for EK class.
    """
    def __init__(self, kappa, phi_min=1.0):
        self.kappa = kappa
        self.phi_min = phi_min

    def pdf(self, x):
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        m = x >= self.phi_min
        out[m] = self.kappa * self.phi_min**self.kappa * x[m]**(-(self.kappa + 1))
        return out

    def cdf(self, x):
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        m = x >= self.phi_min
        out[m] = 1.0 - (self.phi_min / x[m])**self.kappa
        return out

    def sf(self, x):
        x = np.asarray(x, dtype=float)
        out = np.ones_like(x)
        m = x >= self.phi_min
        out[m] = (self.phi_min / x[m])**self.kappa
        return out

    def rvs(self, size, rng=None):
        rng = rng or np.random.default_rng()
        u = rng.uniform(size=size)
        return self.phi_min * (1.0 - u)**(-1.0 / self.kappa)

    def Phi(self, phi_star, sigma):
        s = sigma - 1.0
        phi_star = max(phi_star, self.phi_min)
        return self.kappa / (self.kappa - s) * (self.phi_min / phi_star)**self.kappa

    def eps_melitz(self, phi_star, sigma):
        return float(self.kappa)

    def eps_ek(self, z):
        return np.full_like(np.asarray(z, dtype=float), self.kappa)


# -----------------------------------------------------------------------
# DGP 2: Log-Normal
# -----------------------------------------------------------------------

class LogNormal(BaseDistribution):
    """ln(φ) ~ N(μ, σ_ln²)."""

    def __init__(self, mu=0.0, sigma_ln=1.0):
        self.mu = mu
        self.sigma_ln = sigma_ln

    def pdf(self, x):
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        m = x > 0
        lnx = np.log(x[m])
        out[m] = np.exp(-0.5 * ((lnx - self.mu) / self.sigma_ln)**2) \
                 / (x[m] * self.sigma_ln * np.sqrt(2 * np.pi))
        return out

    def cdf(self, x):
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        m = x > 0
        out[m] = 0.5 * (1 + special.erf(
            (np.log(x[m]) - self.mu) / (self.sigma_ln * np.sqrt(2))))
        return out

    def rvs(self, size, rng=None):
        rng = rng or np.random.default_rng()
        return np.exp(rng.normal(self.mu, self.sigma_ln, size=size))


# -----------------------------------------------------------------------
# DGP 3: Pareto Mixture
# -----------------------------------------------------------------------

class ParetoMixture(BaseDistribution):
    """
    π Pareto(κ₁) + (1-π) Pareto(κ₂), common φ_min.
    Requires κ₁, κ₂ > σ-1.
    """
    def __init__(self, kappa1, kappa2, pi=0.9, phi_min=1.0):
        self.kappa1 = kappa1
        self.kappa2 = kappa2
        self.pi = pi
        self.phi_min = phi_min
        self._c1 = FrechetPareto(kappa1, phi_min)
        self._c2 = FrechetPareto(kappa2, phi_min)

    def pdf(self, x):
        return self.pi * self._c1.pdf(x) + (1 - self.pi) * self._c2.pdf(x)

    def cdf(self, x):
        return self.pi * self._c1.cdf(x) + (1 - self.pi) * self._c2.cdf(x)

    def sf(self, x):
        return self.pi * self._c1.sf(x) + (1 - self.pi) * self._c2.sf(x)

    def rvs(self, size, rng=None):
        rng = rng or np.random.default_rng()
        sel = rng.binomial(1, 1 - self.pi, size=size)
        return np.where(sel == 0, self._c1.rvs(size, rng), self._c2.rvs(size, rng))

    def Phi(self, phi_star, sigma):
        return self.pi * self._c1.Phi(phi_star, sigma) \
             + (1 - self.pi) * self._c2.Phi(phi_star, sigma)

    def eps_melitz(self, phi_star, sigma):
        s = sigma - 1.0
        Phi_val = self.Phi(phi_star, sigma)
        if Phi_val <= 0:
            return np.inf
        return s + phi_star * self.pdf(phi_star) / Phi_val


# -----------------------------------------------------------------------
# DGP 4: Truncated Pareto
# -----------------------------------------------------------------------

class TruncatedPareto(BaseDistribution):
    """
    Pareto(κ) truncated at φ_max.
    G(φ) = [1-(φ_min/φ)^κ] / [1-(φ_min/φ_max)^κ],  φ ∈ [φ_min, φ_max].
    """
    def __init__(self, kappa, phi_min=1.0, phi_max=1e4):
        self.kappa = kappa
        self.phi_min = phi_min
        self.phi_max = phi_max
        self._Z = 1.0 - (phi_min / phi_max)**kappa

    def pdf(self, x):
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        m = (x >= self.phi_min) & (x <= self.phi_max)
        out[m] = self.kappa * self.phi_min**self.kappa \
                 * x[m]**(-(self.kappa + 1)) / self._Z
        return out

    def cdf(self, x):
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x)
        m = (x >= self.phi_min) & (x <= self.phi_max)
        out[m] = (1.0 - (self.phi_min / x[m])**self.kappa) / self._Z
        out[x > self.phi_max] = 1.0
        return out

    def rvs(self, size, rng=None):
        rng = rng or np.random.default_rng()
        u = rng.uniform(size=size)
        return self.phi_min * (1.0 - u * self._Z)**(-1.0 / self.kappa)

    def Phi(self, phi_star, sigma):
        s = sigma - 1.0
        k = self.kappa
        pm, px = self.phi_min, self.phi_max
        phi_star = max(phi_star, pm)
        if phi_star >= px:
            return 0.0
        return k / (self._Z * (k - s)) \
               * ((pm / phi_star)**k - (pm / px)**k * (px / phi_star)**s)

    def eps_melitz(self, phi_star, sigma):
        s = sigma - 1.0
        Phi_val = self.Phi(phi_star, sigma)
        if Phi_val <= 0:
            return np.inf
        return s + phi_star * self.pdf(phi_star) / Phi_val
