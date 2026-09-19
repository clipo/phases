# ABC-SMC inference of the transmission bias

Basin curated set (n = 29), 10 decorated types, 6 ordinal bins. ABC-SMC with local-linear regression adjustment, 800 particles, 6 rounds (full config). Same model and summary statistics as 19_abc_transmission.py.

## Posterior of the transmission-bias parameter b

- Posterior mean b = -0.003, 95% interval [-0.096, +0.031].
- The 95% interval INCLUDES the neutral value b = 0.
- P(b > 0) = 0.51 (0.5 = no directional information).
- Posterior SD 0.025 vs prior SD 0.289 (ratio 0.09).

## Early vs late halves

- Early-half b = -0.023 [-0.083, +0.020].
- Late-half b = +0.018 [-0.011, +0.054].
- Shift early->late: +0.040.

## Round diagnostics (whole-sequence target)

- Tolerance schedule: [1.5081, 1.0166, 0.9869, 0.9129, 0.7173, 0.5187].
- Acceptance rates: [1.0, 0.4016, 0.229, 0.1898, 0.1512, 0.1372].
- Simulations per round: [800, 1992, 3493, 4214, 5292, 5833].
