"""Beta-particle transmission through air, via an allowed-spectrum model and
the empirical range-energy relation for electrons.
"""
import numpy as np

ELECTRON_MASS_MEV = 0.511


def range_to_energy(R_g_cm2):
    """Invert the empirical range-energy relation R = 0.11*(sqrt(1+22.4E^2)-1).

    R is in g/cm^2, returned energy is the kinetic energy in MeV.
    """
    R = np.asarray(R_g_cm2, dtype=float)
    val = ((1 + R / 0.11) ** 2 - 1) / 22.4
    val = np.maximum(val, 0.0)
    return np.sqrt(val)


def allowed_spectrum(E, E0, me=ELECTRON_MASS_MEV):
    """Unnormalised allowed beta spectrum shape (Fermi function ignored)."""
    E = np.asarray(E, dtype=float)
    W = E + me
    p = np.sqrt(np.maximum(W**2 - me**2, 0.0))
    return p * W * (E0 - E) ** 2 * (E >= 0) * (E <= E0)


class BetaEndpoint:
    """Precomputed CDF for an allowed beta spectrum with endpoint energy E0,
    for fast lookups of the survival fraction above an energy cut.
    """

    def __init__(self, E0, n_grid=6000):
        self.E0 = E0
        E = np.linspace(0.0, E0, n_grid)
        s = allowed_spectrum(E, E0)
        dE = np.diff(E)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (s[:-1] + s[1:]) * dE)])
        self.E_grid = E
        self.cdf = cum / cum[-1] if cum[-1] > 0 else np.zeros_like(cum)

    def survival_fraction(self, E_cut):
        """Fraction of the spectrum above E_cut."""
        E_cut = np.asarray(E_cut, dtype=float)
        E_clip = np.clip(E_cut, 0.0, self.E0)
        below = np.interp(E_clip, self.E_grid, self.cdf)
        above = 1.0 - below
        above = np.where(E_cut <= 0.0, 1.0, above)
        above = np.where(E_cut >= self.E0, 0.0, above)
        return above


class AirTransmissionModel:
    """Beta transmission through air for a two-endpoint mixture (e.g. the
    Sr-90/Y-90 decay chain), using the range-energy relation to convert an
    air path length into an energy cutoff.
    """

    def __init__(self, e_low, e_high, rho_air=0.001204, n_grid_low=5000, n_grid_high=8000):
        self.rho_air = rho_air
        self.low = BetaEndpoint(e_low, n_grid_low)
        self.high = BetaEndpoint(e_high, n_grid_high)

    def transmission(self, path_length_cm, f_low):
        """Combined transmission T(d) = f_low*T_low(d) + (1-f_low)*T_high(d).

        Returns (T_combined, T_low, T_high).
        """
        sigma_air = self.rho_air * np.asarray(path_length_cm, dtype=float)
        e_cut = range_to_energy(sigma_air)
        t_low = self.low.survival_fraction(e_cut)
        t_high = self.high.survival_fraction(e_cut)
        return f_low * t_low + (1 - f_low) * t_high, t_low, t_high


def logistic(p):
    """Map real -> (0,1), used to keep a fitted mixing weight f_low in range."""
    return 1.0 / (1.0 + np.exp(-p))
