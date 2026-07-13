# Bayesian F_ST validation battery (Balding-Nichols model)

Config: FULL (f_grid=[0.02, 0.05, 0.1, 0.2, 0.35, 0.5], reps=40, draws=1000 x tune=1000 x chains=2; large sizes=[200, 200, 200], small sizes=[15, 15, 15]; SBC n=200).

## 1-2. Coverage (large sizes)

Coverage of the BN F parameter is checked against the generating F; coverage of the Gini-Simpson readout is checked against the REALIZED Gini-Simpson F_ST of each simulated dataset.

| true F_ST | mean realized G-S F_ST | reps | coverage (BN F) | coverage (G-S readout) |
|---|---|---|---|---|
| 0.020 | 0.0112 | 40 | 1.000 | 0.975 |
| 0.050 | 0.0279 | 40 | 1.000 | 0.950 |
| 0.100 | 0.0664 | 40 | 0.950 | 0.900 |
| 0.200 | 0.1221 | 40 | 0.950 | 0.925 |
| 0.350 | 0.2986 | 40 | 0.975 | 0.950 |
| 0.500 | 0.3486 | 40 | 0.650 | 0.950 |

## 3. Small-sample bias vs the realized Gini-Simpson F_ST (sizes = [15, 15, 15])

| true F_ST | mean realized | reps | plug-in bias | posterior bias |
|---|---|---|---|---|
| 0.020 | 0.0169 | 40 | +0.0428 | +0.0192 |
| 0.050 | 0.0300 | 40 | +0.0478 | +0.0207 |
| 0.100 | 0.0728 | 40 | +0.0506 | +0.0163 |
| 0.200 | 0.1740 | 40 | +0.0406 | -0.0016 |
| 0.350 | 0.3141 | 40 | +0.0171 | -0.0297 |
| 0.500 | 0.3214 | 40 | +0.0578 | +0.0216 |

Mean absolute bias at small sizes: plug-in = 0.0621, posterior median = 0.0481 (posterior closer to truth).

## 4. Simulation-based calibration (BN F parameter)

- 200 datasets drawn from the prior; rank of true F within its posterior binned into 20 bins.
- Chi-square uniformity: chi2 = 19.60, p = 0.419 (want p > 0.05, i.e. not distinguishable from uniform).
- KS uniformity: D = 0.050, p = 0.680.

## 5. MCMC health on the real St. Francis basin fit

- 5 spatial clusters, 10 decorated types; default target_accept (no funnel-taming needed).
- max R-hat = 1.0007 (want < 1.01); min ESS = 4672; divergences = 0.
- BN F_ST median = 0.0665, 95% CI [0.0341, 0.1140].
- Gini-Simpson readout median = 0.0327, 95% CI [0.0301, 0.0352]; plug-in = 0.0330 (inside the interval).

## 6. Bayes-factor sanity

- Panmixia-generated data: 2 ln BF10 = -8.61 (strong for panmixia).
- Strongly-structured data (F=0.30): 2 ln BF10 = 392.29 (very strong for structure).
