#!/usr/bin/env python3
"""Measure detector dead time from oscilloscope pulse traces (D*.CSV), by
detecting each pulse's start/end time and reporting the mean duration.

Usage:
    python scripts/measure_deadtime.py [--data-dir data/deadtime]
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from radiation.pulses import guess_time_voltage_columns, pulse_start_end


def _load_trace(path):
    df = pd.read_csv(path)
    t_col, v_col, df = guess_time_voltage_columns(df)
    t = pd.to_numeric(df[t_col], errors="coerce").to_numpy()
    v = pd.to_numeric(df[v_col], errors="coerce").to_numpy()
    mask = np.isfinite(t) & np.isfinite(v)
    return t[mask], v[mask]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="data/deadtime")
    parser.add_argument("--out-dir", default="results")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    fig_dir = Path(args.out_dir) / "figures" / "deadtime"
    table_dir = Path(args.out_dir) / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(data_dir.glob("D*.CSV")) + sorted(data_dir.glob("D*.csv"))
    if not files:
        raise FileNotFoundError(f"No D*.CSV files found in {data_dir}")

    rows = []
    for f in files:
        t, v = _load_trace(f)
        result = pulse_start_end(t, v)

        plt.figure()
        plt.plot(t * 1e6, v, linewidth=1)
        plt.xlabel("Time (us)")
        plt.ylabel("Voltage (V)")
        plt.title(f.name)
        plt.grid(True)
        if result["ok"]:
            plt.axvline(result["t_start"] * 1e6, linestyle="--")
            plt.axvline(result["t_end"] * 1e6, linestyle="--")
            plt.axhline(result["base_pre"], linestyle=":")
            plt.axhline(result["base_post"], linestyle=":")
        plt.tight_layout()
        plt.savefig(fig_dir / f"{f.stem}.png", dpi=200)
        plt.close()

        duration_us = (result["t_end"] - result["t_start"]) * 1e6 if result["ok"] else np.nan
        rows.append({"file": f.name, "ok": result["ok"], "duration_us": duration_us})

    plt.figure()
    for f in files:
        t, v = _load_trace(f)
        order = np.argsort(t)
        plt.plot(t[order] * 1e6, v[order], linewidth=1, label=f.stem)
    plt.xlabel("Time (us)")
    plt.ylabel("Voltage (V)")
    plt.title("Overlay: pulse traces")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_dir / "overlay.png", dpi=200)
    plt.close()

    summary = pd.DataFrame(rows).sort_values("file")
    summary.to_csv(table_dir / "peak_durations_summary.csv", index=False)
    print(summary.to_string(index=False))

    durations = summary.loc[summary["ok"], "duration_us"].to_numpy()
    tau_mean = durations.mean()
    tau_sem = durations.std(ddof=1) / np.sqrt(len(durations))
    print(f"\nMean pulse duration (dead time) tau = {tau_mean:.3f} +/- {tau_sem:.3f} us")


if __name__ == "__main__":
    main()
