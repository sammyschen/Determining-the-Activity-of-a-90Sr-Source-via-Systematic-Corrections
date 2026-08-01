"""Paralyzable detector dead-time model."""
import numpy as np
from scipy.special import lambertw


def observed_rate_paralyzable(r_true, tau):
    """Forward paralyzable dead-time model: r_obs = r_true * exp(-r_true*tau)."""
    r_true = np.asarray(r_true, dtype=float)
    return r_true * np.exp(-r_true * tau)


def true_rate_paralyzable(r_obs, tau):
    """Invert r_obs = r*exp(-r*tau) -> r = -W(-tau*r_obs)/tau."""
    r_obs = np.asarray(r_obs, dtype=float)
    if tau == 0:
        return r_obs
    z = -tau * r_obs
    return -lambertw(z).real / tau


def dtrue_drobs_paralyzable(r_obs, r_true, tau):
    """Derivative dr_true/dr_obs, for uncertainty propagation through the inversion."""
    if tau == 0:
        return np.ones_like(np.asarray(r_obs, dtype=float))
    denom = 1.0 - r_true * tau
    denom = np.where(np.abs(denom) < 1e-12, np.nan, denom)
    return np.exp(r_true * tau) / denom
