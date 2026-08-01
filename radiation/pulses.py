"""Detecting pulse start/end times in oscilloscope traces, for dead-time
measurement from Geiger-counter output pulses.
"""
import numpy as np
import pandas as pd


def guess_time_voltage_columns(df):
    """Return (time_col, voltage_col, df) for a scope-trace CSV, preferring the
    known LeCroy-style header, otherwise falling back to the first two numeric
    columns.
    """
    cols = [c.strip() for c in df.columns]
    df = df.copy()
    df.columns = cols
    if "in s" in cols and "C1 in V" in cols:
        return "in s", "C1 in V", df
    numeric_cols = []
    for c in cols:
        try:
            pd.to_numeric(df[c].iloc[:50], errors="raise")
            numeric_cols.append(c)
        except Exception:
            pass
    if len(numeric_cols) >= 2:
        return numeric_cols[0], numeric_cols[1], df
    raise ValueError("Could not infer time/voltage columns")


def _first_sustained(mask, n, frac):
    n = max(2, int(n))
    counts = np.convolve(mask.astype(int), np.ones(n, dtype=int), mode="valid")
    hits = np.where(counts >= frac * n)[0]
    return None if len(hits) == 0 else int(hits[0])


def pulse_start_end(
    t,
    v,
    pre_plateau_us=1.0,
    post_plateau_us=1.0,
    k=5,
    start_sustain_ns=80,
    end_sustain_us=0.8,
    frac_start=0.95,
    frac_end=0.98,
):
    """Detect a single pulse's start/end times in an oscilloscope trace, by
    finding a sustained crossing above the pre-pulse baseline and a sustained
    return to the post-pulse baseline.

    Returns a dict with keys: ok, t_start, t_end, base_pre, base_post.
    """
    order = np.argsort(t)
    t, v = t[order], v[order]
    dt = float(np.median(np.diff(t)))
    peak_idx = int(np.argmax(v))
    tmin, tmax = float(t.min()), float(t.max())

    pre_mask = t < (tmin + pre_plateau_us * 1e-6)
    post_mask = t > (tmax - post_plateau_us * 1e-6)

    base_pre = float(np.median(v[pre_mask])) if pre_mask.sum() > 10 else float(np.median(v))
    sig_pre = float(np.std(v[pre_mask], ddof=1)) if pre_mask.sum() > 10 else float(np.std(v, ddof=1))
    base_post = float(np.median(v[post_mask])) if post_mask.sum() > 10 else float(np.median(v))
    sig_post = float(np.std(v[post_mask], ddof=1)) if post_mask.sum() > 10 else float(np.std(v, ddof=1))

    n_start = int(round((start_sustain_ns * 1e-9) / dt))
    n_end = int(round((end_sustain_us * 1e-6) / dt))

    thr_start = base_pre + k * sig_pre
    start_idx = _first_sustained(v > thr_start, n_start, frac_start)

    low, high = base_post - k * sig_post, base_post + k * sig_post
    within_post = (v >= low) & (v <= high)
    within_post[:peak_idx] = False
    end_idx = _first_sustained(within_post, n_end, frac_end)

    ok = start_idx is not None and end_idx is not None
    return {
        "ok": ok,
        "t_start": t[start_idx] if ok else np.nan,
        "t_end": t[end_idx] if ok else np.nan,
        "base_pre": base_pre,
        "base_post": base_post,
    }
