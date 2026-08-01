"""Nested models for u(d) = rate(d) * d^2, from a bare inverse-square law (M0)
up to a full model with geometric offset, solid angle, dead time, and air
transmission (M4). Each fn(beta, x) matches the signature scipy.odr.Model expects.
"""
from dataclasses import dataclass, field

import numpy as np

from .beta_spectrum import AirTransmissionModel, logistic
from .deadtime import observed_rate_paralyzable
from .geometry import solid_angle_disk


@dataclass
class DetectorConfig:
    tau_s: float = 3.12e-6
    a_det_cm2: float = 1.0
    rho_air_g_cm3: float = 0.001204
    beta_e_low_mev: float = 0.54
    beta_e_high_mev: float = 2.27

    a_nom_cm: float = field(init=False)
    air_model: AirTransmissionModel = field(init=False)

    def __post_init__(self):
        self.a_nom_cm = np.sqrt(self.a_det_cm2 / np.pi)
        self.air_model = AirTransmissionModel(
            self.beta_e_low_mev, self.beta_e_high_mev, rho_air=self.rho_air_g_cm3
        )


def m0_constant(beta, x):
    (K,) = beta
    return K * np.ones_like(x)


def m1_offset(beta, x):
    K, d0 = beta
    return K * (x / (x + d0)) ** 2


def m2_solid_angle(beta, x, cfg: DetectorConfig):
    K, d0 = beta
    omega = solid_angle_disk(x + d0, cfg.a_nom_cm)
    return K * x**2 * omega / cfg.a_det_cm2


def m3_dead_time(beta, x, cfg: DetectorConfig):
    K, d0 = beta
    omega = solid_angle_disk(x + d0, cfg.a_nom_cm)
    r_det = K * omega / cfg.a_det_cm2
    r_obs = observed_rate_paralyzable(r_det, cfg.tau_s)
    return x**2 * r_obs


def m4_air_transmission(beta, x, cfg: DetectorConfig):
    K, d0, p = beta
    f_low = logistic(p)
    rsep = x + d0
    omega = solid_angle_disk(rsep, cfg.a_nom_cm)
    T, _, _ = cfg.air_model.transmission(rsep, f_low)
    r_det = K * omega / cfg.a_det_cm2 * T
    r_obs = observed_rate_paralyzable(r_det, cfg.tau_s)
    return x**2 * r_obs
