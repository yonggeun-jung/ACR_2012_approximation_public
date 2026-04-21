"""Reusable ACR diagnostic for gains-from-trade bounds.

This module packages the core diagnostic into a single function so users can
apply it to their own bilateral trade data without running the full pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import statsmodels.api as sm
from statsmodels.genmod.families import Poisson


@dataclass
class AcrBoundResult:
    epsilon_pool: float
    beta_2: float
    p_value: float
    delta_epsilon_hat: float
    bound_ek: float
    bound_melitz: float
    ln_distance_min: float
    ln_distance_max: float
    n_obs: int

    def to_dict(self) -> dict[str, float]:
        return {
            "epsilon_pool": self.epsilon_pool,
            "beta_2": self.beta_2,
            "p_value": self.p_value,
            "delta_epsilon_hat": self.delta_epsilon_hat,
            "bound_ek": self.bound_ek,
            "bound_melitz": self.bound_melitz,
            "ln_distance_min": self.ln_distance_min,
            "ln_distance_max": self.ln_distance_max,
            "n_obs": self.n_obs,
        }


def _as_1d(name: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 0:
        raise ValueError(f"{name} must be non-empty.")
    return arr


def _build_fe_dummies(ids: Any) -> np.ndarray:
    arr = np.asarray(ids).reshape(-1)
    levels, inv = np.unique(arr, return_inverse=True)
    if levels.size <= 1:
        return np.empty((arr.size, 0))
    dummies = np.zeros((arr.size, levels.size - 1), dtype=float)
    mask = inv > 0
    dummies[mask, inv[mask] - 1] = 1.0
    return dummies


def _stack_features(*parts: np.ndarray) -> np.ndarray:
    kept = [p for p in parts if p is not None and p.size > 0]
    if not kept:
        raise ValueError("No regressors were provided.")
    return np.column_stack(kept)


def compute_acr_bound(
    trade_flow: Any,
    ln_distance: Any,
    sigma: float,
    ln_lambda_hat_abs: float,
    bilateral_controls: Any | None = None,
    exporter_ids: Any | None = None,
    importer_ids: Any | None = None,
) -> dict[str, float]:
    """Compute the ACR approximation-error diagnostic and implied bounds.

    Parameters
    ----------
    trade_flow
        Bilateral trade flow vector (non-negative).
    ln_distance
        Log distance vector.
    sigma
        Elasticity of substitution used in the welfare formula.
    ln_lambda_hat_abs
        Absolute value of the domestic-share change term, |ln lambda_hat|.
    bilateral_controls
        Optional 2D array of bilateral controls (contiguity, language, etc.).
    exporter_ids, importer_ids
        Optional exporter/importer identifiers to include origin and
        destination fixed effects as dummies.

    Returns
    -------
    dict
        Keys: epsilon_pool, beta_2, p_value, delta_epsilon_hat, bound_ek,
        bound_melitz, ln_distance_min, ln_distance_max, n_obs.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> n = 100
    >>> ln_d = rng.uniform(5.0, 10.0, n)
    >>> trade = np.exp(8.0 - 1.5 * ln_d + rng.normal(0.0, 0.3, n))
    >>> out = compute_acr_bound(
    ...     trade_flow=trade,
    ...     ln_distance=ln_d,
    ...     sigma=4.0,
    ...     ln_lambda_hat_abs=0.05,
    ... )
    >>> out["bound_ek"] > 0
    True
    """
    y = _as_1d("trade_flow", trade_flow)
    ln_d = _as_1d("ln_distance", ln_distance)
    if y.size != ln_d.size:
        raise ValueError("trade_flow and ln_distance must have the same length.")
    if sigma <= 1:
        raise ValueError("sigma must be greater than 1.")
    if ln_lambda_hat_abs < 0:
        raise ValueError("ln_lambda_hat_abs must be non-negative.")

    ctrl = None
    if bilateral_controls is not None:
        ctrl = np.asarray(bilateral_controls, dtype=float)
        if ctrl.ndim == 1:
            ctrl = ctrl.reshape(-1, 1)
        if ctrl.shape[0] != y.size:
            raise ValueError("bilateral_controls must have the same rows as trade_flow.")

    exp_fe = np.empty((y.size, 0))
    imp_fe = np.empty((y.size, 0))
    if exporter_ids is not None:
        exp_fe = _build_fe_dummies(exporter_ids)
        if exp_fe.shape[0] != y.size:
            raise ValueError("exporter_ids must have the same length as trade_flow.")
    if importer_ids is not None:
        imp_fe = _build_fe_dummies(importer_ids)
        if imp_fe.shape[0] != y.size:
            raise ValueError("importer_ids must have the same length as trade_flow.")

    x_linear = sm.add_constant(_stack_features(ln_d, ctrl, exp_fe, imp_fe))
    x_quad = sm.add_constant(_stack_features(ln_d, ln_d**2, ctrl, exp_fe, imp_fe))

    res_linear = sm.GLM(y, x_linear, family=Poisson()).fit(disp=False, cov_type="HC1")
    res_quad = sm.GLM(y, x_quad, family=Poisson()).fit(disp=False, cov_type="HC1")

    epsilon_pool = float(-res_linear.params[1])
    beta_2 = float(res_quad.params[2])
    p_value = float(res_quad.pvalues[2])

    ln_d_min = float(np.min(ln_d))
    ln_d_max = float(np.max(ln_d))
    delta_eps = float(2.0 * abs(beta_2) * (ln_d_max - ln_d_min))

    bound_ek = float((delta_eps / (sigma - 1.0)) * ln_lambda_hat_abs * 100.0)
    if epsilon_pool > 1e-6:
        bound_melitz = float(bound_ek / epsilon_pool)
    else:
        bound_melitz = float("nan")

    return AcrBoundResult(
        epsilon_pool=epsilon_pool,
        beta_2=beta_2,
        p_value=p_value,
        delta_epsilon_hat=delta_eps,
        bound_ek=bound_ek,
        bound_melitz=bound_melitz,
        ln_distance_min=ln_d_min,
        ln_distance_max=ln_d_max,
        n_obs=int(y.size),
    ).to_dict()
