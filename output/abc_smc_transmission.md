# ABC-SMC inference of the transmission bias

Basin curated set (n = 43), 10 decorated types, 6 ordinal bins. ABC-SMC with local-linear regression adjustment, 800 particles, 6 rounds (full config). Same model and summary statistics as 19_abc_transmission.py.

## Posterior of the transmission-bias parameter b

- Posterior mean b = +0.001, 95% interval [-0.026, +0.029].
- The 95% interval INCLUDES the neutral value b = 0.
- P(b > 0) = 0.52 (0.5 = no directional information).
- Posterior SD 0.014 vs prior SD 0.289 (ratio 0.05).

## Early vs late halves

- Early-half b = -0.008 [-0.066, +0.035].
- Late-half b = +0.017 [-0.010, +0.051].
- Shift early->late: +0.025.

## Round diagnostics (whole-sequence target)

- Tolerance schedule: [1.451, 1.0813, 1.0447, 0.9359, 0.7112, 0.5317].
- Acceptance rates: [1.0, 0.3974, 0.2237, 0.1999, 0.1587, 0.1612].
- Simulations per round: [800, 2013, 3576, 4002, 5041, 4962].
