# ABC-SMC inference of the transmission bias

Basin curated set (n = 29), 10 decorated types, 6 ordinal bins. ABC-SMC with local-linear regression adjustment, 800 particles, 6 rounds (full config). Same model and summary statistics as 19_abc_transmission.py.

## Posterior of the transmission-bias parameter b

- Posterior mean b = +0.006, 95% interval [-0.020, +0.035].
- The 95% interval INCLUDES the neutral value b = 0.
- P(b > 0) = 0.70 (0.5 = no directional information).
- Posterior SD 0.014 vs prior SD 0.289 (ratio 0.05).

## Early vs late halves

- Early-half b = -0.008 [-0.071, +0.035].
- Late-half b = +0.012 [-0.015, +0.044].
- Shift early->late: +0.020.

## Round diagnostics (whole-sequence target)

- Tolerance schedule: [1.5062, 1.0175, 0.9916, 0.915, 0.7286, 0.5365].
- Acceptance rates: [1.0, 0.3901, 0.2277, 0.1905, 0.122, 0.1283].
- Simulations per round: [800, 2051, 3514, 4200, 6556, 6237].
