# Imperial Physics Year Two Radiation Experiment

Measures how a ⁹⁰Sr beta source's count rate falls off with distance, and fits
a chain of physical corrections (geometry, dead time, air attenuation) needed
to recover a clean inverse-square law from real detector data.

![Model comparison](results/figures/u_vs_d_model_comparison.png)

*Raw inverse-square-scaled rate `u(d) = rate(d)·d²` (points) against three
nested correction models — a bare constant (red) badly fails both near the
source and far from it; adding geometry, dead time, and air attenuation
(magenta) tracks the data across the full 0–90 cm range.*

## Why it matters

For an ideal point source and point detector, count rate should scale
exactly as `1/d²`, so `u(d) = rate(d)·d²` should be flat. Real measurements
never are: near the source the detector's finite size and its dead time
distort the count rate, and far away beta particles lose energy in air
before reaching the detector. Fitting a bare `1/d²` law without these
corrections silently biases the inferred source activity — the same failure
mode that shows up in any inverse-square measurement (radioactive sources,
light sources, RF power) taken with a real, non-ideal detector. This repo
works through that correction chain end-to-end and quantifies how much each
term actually matters, rather than assuming the ideal law and moving on.

## Key technical features

- **Nested nonlinear model family** (`radiation/models.py`) — five models
  (M0–M4) building up from a bare constant to full geometry + dead-time +
  air-attenuation, fit by orthogonal distance regression (`radiation/fitting.py`)
  so both distance and rate uncertainties are honored, not just rate.
- **Physical beta-spectrum transmission model** (`radiation/beta_spectrum.py`)
  — converts an air path length to an energy cutoff via the empirical
  range-energy relation, then integrates an allowed-spectrum shape (two mixed
  endpoints for the Sr-90/Y-90 decay chain) to get a transmission fraction —
  a physical model, not a fitted empirical curve.
- **Paralyzable dead-time inversion** (`radiation/deadtime.py`) — inverts
  `r_obs = r·e^(-rτ)` in closed form via the Lambert W function, with
  analytic uncertainty propagation through the inversion.
- **Model comparison, not just model fitting** — AIC/BIC and nested Δχ² tests
  (`radiation/fitting.py`) to justify each added correction term rather than
  just reporting the fanciest model's numbers.
- **Independent cross-check** — `scripts/air_transmission_explore.py`
  re-derives the beta-mixing weight `f_low` by a totally separate method
  (grid search on far-field data only) as a sanity check on the main ODR fit.

## Installation and quick start

```bash
pip install -r requirements.txt
# or: pip install -e .
```

```bash
# Fit models M0-M4 to rate-vs-distance data and report the fully
# corrected source constant K_est(d).
python scripts/fit_inverse_square.py --data data/filtered_data.csv

# Independent cross-check: grid-search the beta-spectrum mixing weight
# against far-field data only.
python scripts/air_transmission_explore.py --data data/filtered_data.csv

# Measure detector dead time from oscilloscope pulse traces.
python scripts/measure_deadtime.py --data-dir data/deadtime
```

Each script prints its results and (re)writes the corresponding
figures/tables under `results/`.

## Main results

The five nested models for `u(d) = rate(d)·d²`, which should be flat if the
inverse-square law and no systematics held exactly:

- **M0** — constant (far-field only)
- **M1** — + geometric distance offset `d0`
- **M2** — + finite-detector solid angle
- **M3** — + paralyzable dead-time correction
- **M4** — + beta transmission through air

Dead time, from five oscilloscope pulse traces: **τ = 3.12 ± 0.06 μs**.

Fitting all five to 30 distance points (0–90 cm), the full model (M4) is
overwhelmingly preferred over any model missing the air term:

| Model | k | χ² | dof | χ²_red |
|---|---|---|---|---|
| M1 + offset | 2 | 1065.7 | 28 | 38.1 |
| M2 + solid angle | 2 | 1035.8 | 28 | 37.0 |
| M3 + dead time | 2 | 990.1 | 28 | 35.4 |
| **M4 + air transmission** | **3** | **67.9** | **27** | **2.52** |

Adding the air-transmission term alone (M3 → M4) drops χ² by 922 for one
extra parameter (p ≈ 1×10⁻²⁰²). The fit gives a geometric offset
`d0 = 1.34 ± 0.08 cm` and an effective low-energy beta-spectrum weight
`f_low = 0.169 ± 0.016`.

Backing out the fully corrected source constant point-by-point and taking
its weighted mean gives **K̄ = 151,580 ± 540 cps·cm²**, consistent with a
constant across the full distance range (χ² = 26.4, dof = 29, p = 0.61) with
no significant residual trend with distance (slope p = 0.25):

![Post-correction constancy](results/figures/K_est_constancy_test_plot.png)

With a 1 cm² detector area, this implies `Aη = 4πK̄ ≈ (1.905 ± 0.007)×10⁶ Bq`,
and — relative to the source's labeled 5 MBq activity — an overall detection
efficiency **η ≈ 0.381 ± 0.001**.

## Repository structure

```
radiation/        reusable analysis package (physics models, fitting, I/O)
scripts/          thin CLI entry points that call into radiation/
data/             raw measurements (rate vs. distance, oscilloscope traces)
results/          generated figures and tables, committed for visibility
docs/             lab script and final report PDFs
archive/          original notebooks and first-pass outputs, kept for reference
```

## Tests and reproducibility

- The ODR fit is verified robust, not just lucky: refitting M4 from five
  different initial guesses for the mixing-weight parameter converges to the
  same optimum every time (`K≈151,280–151,300`, `d0≈1.336–1.337`,
  `f_low≈0.169`).
- `scripts/air_transmission_explore.py` cross-checks `f_low` against
  `scripts/fit_inverse_square.py`'s ODR fit using an unrelated method (grid
  search on far-field-only data) rather than trusting a single fitting
  routine.
- `results/` is committed with the exact output of the current code, so a
  stale figure or number shows up as a diff the next time the scripts are
  rerun.
- There is no automated (`pytest`-style) unit test suite yet — see
  Limitations below.

## Limitations and future work

- Detector solid angle is modeled as an equivalent-area **disk**
  (`radiation/geometry.py`), not the exact rectangular-aperture formula for
  the 1cm×1cm detector derived in the lab script's appendix — a reasonable
  approximation at these distances, but not an exact match.
- Even the best model (M4) does not formally pass a χ² goodness-of-fit test
  (χ²_red = 2.52 gives p ≈ 2×10⁻⁵): some smaller, unmodeled systematic
  remains beyond geometry, dead time, and air attenuation.
- Air density and the Sr-90/Y-90 beta endpoint energies are fixed defaults in
  `DetectorConfig`, not fit to ambient conditions or re-derived from the
  source's decay scheme.
- No automated test suite: correctness currently rests on the cross-checks
  described above rather than pinned regression tests. Adding a small
  `pytest` suite (e.g. asserting the ODR fit reproduces known parameters
  within tolerance) would catch silent regressions from future refactors.

## License

MIT — see [LICENSE](LICENSE).
