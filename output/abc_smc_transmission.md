# ABC-SMC inference of the transmission bias

Basin curated set (n = 29), 10 decorated types, 6 ordinal bins. ABC-SMC with local-linear regression adjustment, 800 particles, 6 rounds (full config). Same model and summary statistics as 19_abc_transmission.py.

## Posterior of the transmission-bias parameter b

- Posterior mean b = +0.010, 95% interval [-0.020, +0.052].
- The 95% interval INCLUDES the neutral value b = 0.
- P(b > 0) = 0.77 (0.5 = no directional information).
- Posterior SD 0.016 vs prior SD 0.289 (ratio 0.06).

## Early vs late halves

- Early-half b = -0.016 [-0.065, +0.023].
- Late-half b = +0.013 [-0.015, +0.049].
- Shift early->late: +0.029.

## Round diagnostics (whole-sequence target)

- Tolerance schedule: [1.5062, 1.0175, 0.9916, 0.9103, 0.7178, 0.5257].
- Acceptance rates: [1.0, 0.3901, 0.2306, 0.1945, 0.1437, 0.1339].
- Simulations per round: [800, 2051, 3469, 4113, 5567, 5974].
