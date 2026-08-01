"""ODR fitting of the u(d) = rate * d^2 models M0-M4, model comparison stats,
and reconstruction of the fully-corrected source constant K_est(d).
"""
from functools import partial

import numpy as np
import pandas as pd
from scipy.odr import ODR, Model, RealData
from scipy.stats import chi2 as chi2_dist
from scipy.stats import t as t_dist

from .deadtime import dtrue_drobs_paralyzable, true_rate_paralyzable
from .geometry import solid_angle_disk
from .models import (
    DetectorConfig,
    m0_constant,
    m1_offset,
    m2_solid_angle,
    m3_dead_time,
    m4_air_transmission,
)
from .stats import weighted_linear_fit, weighted_mean

MODEL_ORDER = [
    "M0 constant",
    "M1 + offset",
    "M2 + solid angle",
    "M3 + dead time",
    "M4 + air transmission",
]


def _far_field_k_guess(d, u, su, min_d):
    mask = d >= min_d
    w = 1.0 / su[mask] ** 2
    return float(np.sum(w * u[mask]) / np.sum(w))


def _pack(out, f, x, k):
    # out.sum_square is the ODR objective: the weighted sum of squares in the
    # *orthogonal* sense, i.e. it accounts for uncertainty in both x and y.
    # A naive sum(((y - yhat) / sy) ** 2) ignores the x-uncertainty ODR is
    # explicitly fitting against, and drastically overstates lack of fit for
    # any model with x-dependence, making model comparisons meaningless.
    chi2 = float(out.sum_square)
    dof = len(x) - k
    return {
        "beta": out.beta,
        "sd_beta": out.sd_beta,
        "n": len(x),
        "k": k,
        "chi2": chi2,
        "dof": dof,
        "chi2_red": chi2 / dof,
        "aic": chi2 + 2 * k,
        "bic": chi2 + k * np.log(len(x)),
        "model_fn": f,
    }


def fit_models(d, sd, u, su, cfg: DetectorConfig, far_field_min_d=20.0):
    """Fit models M0 (far-field only) through M4 (all corrections) by ODR
    (uncertainties in both x and y). Returns a dict keyed by model name.
    """
    k_guess = _far_field_k_guess(d, u, su, far_field_min_d)
    results = {}

    far = d >= far_field_min_d
    data_far = RealData(d[far], u[far], sx=sd[far], sy=su[far])
    out0 = ODR(data_far, Model(m0_constant), beta0=[k_guess]).run()
    results["M0 constant"] = _pack(out0, m0_constant, d[far], 1)

    data = RealData(d, u, sx=sd, sy=su)

    out1 = ODR(data, Model(m1_offset), beta0=[k_guess, 0.5]).run()
    results["M1 + offset"] = _pack(out1, m1_offset, d, 2)

    f2 = partial(m2_solid_angle, cfg=cfg)
    out2 = ODR(data, Model(f2), beta0=[k_guess, 0.5]).run()
    results["M2 + solid angle"] = _pack(out2, f2, d, 2)

    f3 = partial(m3_dead_time, cfg=cfg)
    out3 = ODR(data, Model(f3), beta0=[k_guess, 0.5]).run()
    results["M3 + dead time"] = _pack(out3, f3, d, 2)

    f4 = partial(m4_air_transmission, cfg=cfg)
    out4 = ODR(data, Model(f4), beta0=[k_guess, 0.5, -4.0]).run()
    results["M4 + air transmission"] = _pack(out4, f4, d, 3)

    return results


def comparison_table(results, order=MODEL_ORDER):
    rows = [
        {"Model": name, **{k: results[name][k] for k in ("k", "chi2", "dof", "chi2_red", "aic", "bic")}}
        for name in order
    ]
    return pd.DataFrame(rows)


def nested_chi2_tests(results, order=MODEL_ORDER):
    """Delta-chi2 tests between successive nested models (meaningful only when delta_k>0)."""
    rows = []
    for a, b in zip(order[:-1], order[1:]):
        ra, rb = results[a], results[b]
        dchi = ra["chi2"] - rb["chi2"]
        dk = rb["k"] - ra["k"]
        p = chi2_dist.sf(dchi, dk) if (dk > 0 and dchi >= 0) else np.nan
        rows.append({"comparison": f"{a} -> {b}", "delta_chi2": dchi, "delta_k": dk, "p_value": p})
    return pd.DataFrame(rows)


def corrected_constant(d, sd, r, sr, cfg: DetectorConfig, K, d0, f_low, tau):
    """Back out the fully-corrected per-point source constant K_est(d): invert
    dead time, then divide out solid angle and air transmission, propagating
    uncertainty from sigma_r and sigma_d.
    """

    def k_of(d_val, r_obs):
        rsep = d_val + d0
        omega = solid_angle_disk(rsep, cfg.a_nom_cm)
        t, _, _ = cfg.air_model.transmission(rsep, f_low)
        r_det = true_rate_paralyzable(r_obs, tau)
        return r_det * cfg.a_det_cm2 / (omega * t)

    rsep = d + d0
    omega = solid_angle_disk(rsep, cfg.a_nom_cm)
    t_comb, _, _ = cfg.air_model.transmission(rsep, f_low)
    r_det = true_rate_paralyzable(r, tau)
    k_est = r_det * cfg.a_det_cm2 / (omega * t_comb)

    drdrobs = dtrue_drobs_paralyzable(r, r_det, tau)
    s_rdet = np.abs(drdrobs) * sr

    k_plus = k_of(d + sd, r)
    k_minus = k_of(d - sd, r)
    dk_dd = (k_plus - k_minus) / (2 * sd)

    s_k = np.sqrt((cfg.a_det_cm2 / (omega * t_comb) * s_rdet) ** 2 + (dk_dd * sd) ** 2)
    return k_est, s_k


def constancy_and_trend_tests(d, k_est, s_k):
    """Test whether K_est(d) is consistent with a constant (chi2 about the
    weighted mean), and whether it has a residual linear trend with distance.
    """
    k_bar, s_k_bar = weighted_mean(k_est, s_k)
    chi2_k = float(np.sum(((k_est - k_bar) / s_k) ** 2))
    dof_k = len(d) - 1
    p_constancy = chi2_dist.sf(chi2_k, dof_k)

    m, sm, c, sc = weighted_linear_fit(d, k_est, s_k)
    t_stat = m / sm
    p_slope = 2 * t_dist.sf(np.abs(t_stat), len(d) - 2)

    return {
        "k_bar": k_bar,
        "s_k_bar": s_k_bar,
        "chi2_constancy": chi2_k,
        "dof_constancy": dof_k,
        "p_constancy": p_constancy,
        "slope": m,
        "s_slope": sm,
        "intercept": c,
        "s_intercept": sc,
        "p_slope": p_slope,
    }
