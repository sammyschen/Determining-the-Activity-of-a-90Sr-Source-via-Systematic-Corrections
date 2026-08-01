#!/usr/bin/env python3
"""Standalone cross-check of beta transmission through air: grid-search the
low/high-energy mixing weight f_low against far-field u(d) data, independent
of the full ODR model fit in fit_inverse_square.py.

Usage:
    python scripts/air_transmission_explore.py [--data data/filtered_data.csv]
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from radiation.beta_spectrum import AirTransmissionModel
from radiation.data import load_rate_vs_distance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/filtered_data.csv")
    parser.add_argument("--out-dir", default="results/figures")
    parser.add_argument("--far-field-min-d", type=float, default=15.0)
    parser.add_argument("--e-low", type=float, default=0.54)
    parser.add_argument("--e-high", type=float, default=2.27)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ds = load_rate_vs_distance(args.data)
    d, u, su = ds["d"], ds["u"], ds["su"]

    air = AirTransmissionModel(args.e_low, args.e_high)

    f_grid = np.linspace(0, 1, 501)
    far = d >= args.far_field_min_d
    w = 1.0 / su[far] ** 2

    best_chi2, best_f, best_k = np.inf, None, None
    for f in f_grid:
        T, _, _ = air.transmission(d[far], f)
        k = np.sum(w * u[far] * T) / np.sum(w * T**2)
        chi2 = np.sum(w * (u[far] - k * T) ** 2)
        if chi2 < best_chi2:
            best_chi2, best_f, best_k = chi2, f, k

    print(f"Best-fit f_low = {best_f:.3f}")
    print(f"Best-fit constant K = {best_k:.1f} cps.cm^2")

    T_comb, T_low, T_high = air.transmission(d, best_f)
    u_corr = u / T_comb
    su_corr = su / T_comb

    plt.figure()
    plt.plot(d, T_high, "o-", label=f"T_high (Emax={args.e_high} MeV)")
    plt.plot(d, T_low, "o--", label=f"T_low (Emax={args.e_low} MeV)")
    plt.plot(d, T_comb, "o-.", label=f"T_combined (f_low={best_f:.3f})")
    plt.xlabel("Distance d (cm)")
    plt.ylabel("Transmission T(d)")
    plt.title("Beta transmission through air (range-energy model)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "air_transmission_curves.png", dpi=300)
    plt.close()

    plt.figure()
    plt.errorbar(d, u, yerr=su, fmt="o", capsize=3, label="measured u = rate . d^2")
    plt.errorbar(d, u_corr, yerr=su_corr, fmt="o", capsize=3, label="u corrected by /T(d)")
    plt.axhline(best_k, linestyle="--", label=f"best-fit constant K = {best_k:.0f}")
    plt.xlabel("Distance d (cm)")
    plt.ylabel("u (cps.cm^2)")
    plt.title("Effect of air-transmission correction on u(d)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_dir / "air_transmission_correction_effect.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
