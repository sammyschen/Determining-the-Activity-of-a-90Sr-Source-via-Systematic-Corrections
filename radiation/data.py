"""Loading and collapsing rate-vs-distance measurements."""
import numpy as np
import pandas as pd

from .stats import weighted_mean

# Ruler resolution = 1 cm -> standard uncertainty (uniform rounding +/- 0.5 cm)
SIGMA_D_CM_1CM_RULER = 0.5 / np.sqrt(3)


def load_rate_vs_distance(path, sigma_d_cm=SIGMA_D_CM_1CM_RULER):
    """Load a rate-vs-distance CSV and collapse repeat measurements at the same
    distance into a single weighted-mean point.

    Expects a `distance_cm` column, and either
      (rate_cps, sigma_rate_cps) or (u_cps_cm2, sigma_u_cps_cm2).

    Returns a dict of numpy arrays: d, sd, r, sr, u, su.
    """
    df = pd.read_csv(path)
    df = df[df["distance_cm"] > 0].copy()

    use_u_direct = "u_cps_cm2" in df.columns and "sigma_u_cps_cm2" in df.columns
    use_rate = "rate_cps" in df.columns and "sigma_rate_cps" in df.columns
    if not (use_u_direct or use_rate):
        raise ValueError(
            "CSV must contain either (rate_cps, sigma_rate_cps) or "
            "(u_cps_cm2, sigma_u_cps_cm2)."
        )

    rows = []
    for dist, s in df.groupby("distance_cm"):
        if use_u_direct:
            u_m, su_m = weighted_mean(s["u_cps_cm2"].values, s["sigma_u_cps_cm2"].values)
            rows.append({"distance_cm": dist, "u": u_m, "su": su_m})
        else:
            r_m, sr_m = weighted_mean(s["rate_cps"].values, s["sigma_rate_cps"].values)
            rows.append({"distance_cm": dist, "rate_cps": r_m, "sigma_rate_cps": sr_m})

    by = pd.DataFrame(rows).sort_values("distance_cm").reset_index(drop=True)
    d = by["distance_cm"].to_numpy()
    sd = np.full_like(d, sigma_d_cm, dtype=float)

    if use_u_direct:
        u = by["u"].to_numpy()
        su = by["su"].to_numpy()
        r = u / d**2
        sr = su / d**2
    else:
        r = by["rate_cps"].to_numpy()
        sr = by["sigma_rate_cps"].to_numpy()
        u = r * d**2
        su = np.sqrt((d**2 * sr) ** 2 + (2 * d * r * sd) ** 2)

    return {"d": d, "sd": sd, "r": r, "sr": sr, "u": u, "su": su}
