# How reproducible is the ABC transmission posterior?

Produced by `analyses/55_abc_stability.py` (full: 800 particles, 6 rounds, 5 seeds). Chases F22 and settles the dependency for rule-18 item G, which draws population sizes from this posterior.

## Across-seed spread

| seed | mean b | 95% interval | posterior SD | P(b > 0) | N median |
|---|---|---|---|---|---|
| 101 | +0.0035 | [-0.0352, +0.0388] | 0.0171 | **0.624** | 222 |
| 202 | +0.0008 | [-0.0315, +0.0356] | 0.0168 | **0.525** | 139 |
| 303 | +0.0043 | [-0.0374, +0.0394] | 0.0179 | **0.633** | 91 |
| 404 | +0.0038 | [-0.0323, +0.0382] | 0.0165 | **0.589** | 104 |
| 505 | +0.0024 | [-0.0329, +0.0361] | 0.0160 | **0.602** | 137 |

- P(b > 0): range **0.525 to 0.633**, SD across seeds 0.038.
- mean b: range +0.0008 to +0.0043, SD 0.0013.
- posterior SD: range 0.0160 to 0.0179.

## Verdict

**P(b > 0) is not stable to the precision it is quoted at.** Across seeds it spans 0.525 to 0.633, a range of 0.108, against a gap of 0.098 between the supplement's 0.70 and the median of these runs. The supplement quotes a seed-dependent summary to two decimal places.

**The remedy is not to pick a seed.** P(b > 0) is a tail mass evaluated where this posterior is densest, so it is the least stable summary available and the one most sensitive to the particle sample. The mean and the interval, which are integrals over the whole posterior, are comparatively steady and reproduce the supplement closely. Report the interval and drop the point probability, or report the probability with its across-seed range attached.

## The dependency for item G

Item G draws the population size from this posterior. Its median across seeds runs 91 to 222, a spread of 94 percent of the mean. The neutrality posterior predictive check must therefore either pool draws across seeds or report its envelope's sensitivity to the seed, rather than conditioning on one ABC run as if it were the posterior.
