#!/usr/bin/env python3
"""Fit systematic-correction models M0-M4 to u(d) = rate * d^2 for the Sr-90
inverse-square-law measurement, and report the fully-corrected source
constant K_est(d).

Usage:
    python scripts/fit_inverse_square.py [--data data/filtered_data.csv]
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from radiation.beta_spectrum import logistic
from radiation.data import load_rate_vs_distance
from radiation.fitting import (
    comparison_table,
    constancy_and_trend_tests,
    corrected_constant,
    fit_models,
    nested_chi2_tests,
)
from radiation.models import DetectorConfig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="data/filtered_data.csv")
    parser.add_argument("--out-dir", default="results")
    parser.add_argument("--far-field-min-d", type=float, default=20.0)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)

    ds = load_rate_vs_distance(args.data)
    d, sd, r, sr, u, su = ds["d"], ds["sd"], ds["r"], ds["sr"], ds["u"], ds["su"]

    cfg = DetectorConfig()
    results = fit_models(d, sd, u, su, cfg, far_field_min_d=args.far_field_min_d)

    table = comparison_table(results).sort_values("aic")
    print("\n=== Model fit comparison (sorted by AIC) ===")
    print(table.to_string(index=False))
    table.to_csv(out_dir / "tables" / "model_fit_comparison_stats.csv", index=False)

    tests = nested_chi2_tests(results)
    print("\n=== Nested delta-chi2 tests (meaningful only when delta_k>0) ===")
    print(tests.to_string(index=False))
    tests.to_csv(out_dir / "tables" / "nested_model_delta_chi2_tests.csv", index=False)

    best = results["M4 + air transmission"]
    K, d0, p = best["beta"]
    f_low = float(logistic(p))
    print("\n=== M4 best-fit parameters ===")
    print(f"K = {K:.1f} cps.cm^2, d0 = {d0:.3f} cm, f_low = {f_low:.4f}")

    k_est, s_k = corrected_constant(d, sd, r, sr, cfg, K, d0, f_low, cfg.tau_s)
    stats = constancy_and_trend_tests(d, k_est, s_k)
    print("\n=== Post-correction constancy tests ===")
    print(f"Weighted mean K_bar = {stats['k_bar']:.0f} +/- {stats['s_k_bar']:.0f} cps.cm^2")
    print(
        f"Constancy chi2 = {stats['chi2_constancy']:.2f} for dof={stats['dof_constancy']} "
        f"-> p = {stats['p_constancy']:.3f}"
    )
    print(
        f"Slope = {stats['slope']:.2f} +/- {stats['s_slope']:.2f} (cps.cm^2)/cm "
        f"-> p = {stats['p_slope']:.3f}"
    )

    pd.DataFrame({"distance_cm": d, "K_est": k_est, "sigma_K": s_k}).to_csv(
        out_dir / "tables" / "K_est_all_corrections.csv", index=False
    )

    _plot_model_comparison(d, sd, u, su, results, out_dir / "figures" / "u_vs_d_model_comparison.png")
    _plot_k_est(d, sd, k_est, s_k, stats, out_dir / "figures" / "K_est_constancy_test_plot.png")


def _plot_model_comparison(d, sd, u, su, results, path):
    d_grid = np.linspace(d.min(), d.max(), 800)
    plt.figure()
    plt.errorbar(d, u, xerr=sd, yerr=su, fmt="o", capsize=3, label="data: u(d) = rate . d^2")
    m0, m3, m4 = results["M0 constant"], results["M3 + dead time"], results["M4 + air transmission"]
    plt.plot(d_grid, m0["model_fn"](m0["beta"], d_grid), label="M0: constant", color="red")
    plt.plot(d_grid, m3["model_fn"](m3["beta"], d_grid), label="M3: geometry + dead time", color="cyan")
    plt.plot(d_grid, m4["model_fn"](m4["beta"], d_grid), label="M4: + air transmission", color="magenta")
    plt.xlabel("distance d (cm)")
    plt.ylabel(r"$u(d)=(n/\Delta t)\,d^2$ (cps.cm$^2$)")
    plt.title("Systematic-correction model comparison")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def _plot_k_est(d, sd, k_est, s_k, stats, path):
    d_grid = np.linspace(d.min(), d.max(), 300)
    k_line = stats["slope"] * d_grid + stats["intercept"]
    plt.figure()
    plt.errorbar(d, k_est, xerr=sd, yerr=s_k, fmt="o", capsize=3, label="K_est(d) after all corrections")
    plt.plot(d_grid, k_line, color="red", label="WLS trend")
    plt.axhline(stats["k_bar"], linestyle="--", label="weighted mean")
    plt.xlabel("distance d (cm)")
    plt.ylabel("K_est (cps.cm^2)")
    plt.title("Post-correction constancy test")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
