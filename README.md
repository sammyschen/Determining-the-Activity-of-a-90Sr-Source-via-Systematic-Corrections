# Radiation

Analysis of a Sr-90 beta source's count rate vs. distance, testing the
inverse-square law and correcting for detector dead time, solid-angle
geometry, and beta attenuation in air.

## Layout

```
radiation/        reusable analysis package (physics models, fitting, I/O)
scripts/          thin CLI entry points that call into radiation/
data/             raw measurements (rate vs. distance, oscilloscope traces)
results/          generated figures and tables (not committed, see below)
docs/             lab script and final report PDFs
archive/          original notebooks and first-pass outputs, kept for reference
```

## Usage

```bash
# Fit models M0 (bare inverse-square) through M4 (+ offset, solid angle,
# dead time, air transmission) to rate-vs-distance data, and report the
# fully corrected source constant K_est(d).
python scripts/fit_inverse_square.py --data data/filtered_data.csv

# Standalone cross-check: grid-search the beta-spectrum mixing weight
# against far-field data only.
python scripts/air_transmission_explore.py --data data/filtered_data.csv

# Measure detector dead time from oscilloscope pulse traces.
python scripts/measure_deadtime.py --data-dir data/deadtime
```

Each script prints its results and writes figures/tables under `results/`.

## Physics models

`radiation/models.py` defines five nested models for
`u(d) = rate(d) * d^2`, which should be flat if the inverse-square law and no
systematics held exactly:

- **M0** — constant (far-field only)
- **M1** — + geometric distance offset `d0`
- **M2** — + finite-detector solid angle
- **M3** — + paralyzable dead-time correction
- **M4** — + beta transmission through air (two-endpoint allowed-spectrum
  model for the Sr-90/Y-90 decay chain)

`radiation/fitting.py` fits all five by orthogonal distance regression
(uncertainties in both distance and rate), compares them by AIC/BIC and
nested delta-chi2 tests, and reconstructs the corrected constant `K_est(d)`
with propagated uncertainty, testing it for constancy and residual trend.

## Notes

- `results/` is committed with the current run's figures/tables so the
  analysis is visible without executing anything; rerunning the scripts
  regenerates them in place.
- `archive/` holds the original exploratory notebooks and first-pass dead-time
  outputs the scripts above were distilled from.

## License

MIT — see [LICENSE](LICENSE).
