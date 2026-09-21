# ABC inference of the transmission bias (dynamically-sufficient test)

Basin curated set (n = 43), 10 decorated types, 6 ordinal bins. 40000 prior draws, closest 400 accepted. Time-averaging window inferred.

## Posterior of the transmission-bias parameter b

- Posterior mean b = +0.029, 95% interval [-0.028, +0.065].
- The 95% interval INCLUDES the neutral value b = 0.
- P(b > 0) = 0.91 (0.5 = no directional information).
- Posterior SD 0.025 vs prior SD 0.287 (ratio 0.09; much less than 1 means the data strongly constrain b rather than leaving it at the prior).

## Early vs late halves (transition / mode-shift test)

- Early-half posterior b = +0.005 [-0.091, +0.051].
- Late-half posterior b = +0.043 [+0.002, +0.082].
- Shift early->late: +0.037; no resolved shift (intervals overlap).

## Verdict

The posterior of the transmission bias is tightly constrained near neutral (95% interval [-0.028, +0.065], including b = 0; posterior SD 0.025 vs prior 0.287), which rules out the strong conformist bias a marker or assortment dynamic would require. A weak conformist lean (P(b>0) = 0.91) is not resolved as a departure from neutral, and no shift between the early and late halves appears. The generative, time-averaging-aware inference therefore reaches the same conclusion as the static convergence criterion, by modeling the transformation between states rather than the states alone, and with a posterior that is informative rather than merely broad.