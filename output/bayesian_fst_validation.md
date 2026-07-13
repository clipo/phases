# Bayesian F_ST validation battery

Config: FULL (f_grid=[0.02, 0.05, 0.1, 0.2, 0.35, 0.5], reps=40, draws=1000 x tune=1000 x chains=4; large sizes=[200, 200, 200], small sizes=[15, 15, 15]).

## 1. Coverage of the 95% credible interval (large sizes)

Coverage is assessed against the REALIZED Gini-Simpson F_ST of each simulated dataset (the quantity the model estimates), not the F-model concentration parameter used to generate it.

| F-model F_ST | mean realized F_ST | reps | coverage |
|---|---|---|---|
| 0.020 | 0.0137 | 40 | 0.975 |
| 0.050 | 0.0377 | 40 | 0.925 |
| 0.100 | 0.0606 | 40 | 0.950 |
| 0.200 | 0.1387 | 40 | 1.000 |
| 0.350 | 0.2856 | 40 | 0.950 |
| 0.500 | 0.4329 | 40 | 0.975 |

Nominal target is 0.95 coverage at each grid point. With few reps (fast mode) coverage estimates are noisy; the full battery should be run for a publication-grade calibration check.

## 2. Small-sample bias vs the realized F_ST (sizes = [15, 15, 15])

| F-model F_ST | mean realized F_ST | reps | coverage | plug-in bias | posterior bias |
|---|---|---|---|---|---|
| 0.020 | 0.0126 | 40 | 1.000 | +0.0403 | +0.0084 |
| 0.050 | 0.0393 | 40 | 0.875 | +0.0397 | -0.0032 |
| 0.100 | 0.0608 | 40 | 0.875 | +0.0507 | -0.0039 |
| 0.200 | 0.1358 | 40 | 0.700 | +0.0312 | -0.0300 |
| 0.350 | 0.2592 | 40 | 0.800 | +0.0298 | -0.0435 |
| 0.500 | 0.3680 | 40 | 0.925 | +0.0356 | -0.0345 |

Mean absolute bias at small sizes (vs the realized F_ST): plug-in = 0.0590, posterior mean = 0.0434 (posterior closer to truth).

## 3. MCMC health on the real St. Francis basin fit

- 5 spatial clusters, 10 decorated types, draws=1000 x tune=1000 x chains=4.
- max R-hat = 1.0086 (want < 1.01).
- min ESS = 835.
- divergences = 0.

Divergences, if any, are not addressed here: reparameterization or prior tuning is a separate follow-on task; this script only reports the diagnostic.

## 4. Frequentist cross-read on the real basin data

- Posterior mean F_ST = 0.0327, 95% credible interval [0.0302, 0.0353].
- Frequentist plug-in (variance.cultural_fst, direct) = 0.0330.
- Plug-in falls inside the 95% credible interval.
