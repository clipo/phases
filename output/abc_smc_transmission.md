# ABC-SMC inference of the transmission bias

Basin curated set (n = 28), 10 decorated types, 6 ordinal bins. ABC-SMC with local-linear regression adjustment, 800 particles, 6 rounds (full config). Same model and summary statistics as 19_abc_transmission.py.

## Posterior of the transmission-bias parameter b

- Posterior mean b = +0.003, 95% interval [-0.035, +0.039].
- The 95% interval INCLUDES the neutral value b = 0.
- P(b > 0) = 0.62 (0.5 = no directional information).
- Posterior SD 0.017 vs prior SD 0.289 (ratio 0.06).

## Early vs late halves

- Early-half b = -0.000 [-0.071, +0.045].
- Late-half b = +0.012 [-0.017, +0.046].
- Shift early->late: +0.012.

## Round diagnostics (whole-sequence target)

- Tolerance schedule: [1.5998, 0.9138, 0.8902, 0.8298, 0.6798, 0.4995].
- Acceptance rates: [1.0, 0.4111, 0.2308, 0.1828, 0.1354, 0.1053].
- Simulations per round: [800, 1946, 3466, 4376, 5909, 7600].
