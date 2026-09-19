# How reproducible is the ABC transmission posterior?

Produced by `analyses/55_abc_stability.py` (full: 800 particles, 6 rounds, 5 seeds). Chases F22 and settles the dependency for rule-18 item G, which draws population sizes from this posterior.

## Across-seed spread

| seed | mean b | 95% interval | posterior SD | P(b > 0) | N median |
|---|---|---|---|---|---|
| 101 | +0.0101 | [-0.0195, +0.0518] | 0.0159 | **0.772** | 176 |
| 202 | +0.0109 | [-0.0152, +0.0409] | 0.0135 | **0.822** | 87 |
| 303 | +0.0030 | [-0.0273, +0.0356] | 0.0157 | **0.590** | 129 |
| 404 | +0.0023 | [-0.0234, +0.0329] | 0.0137 | **0.556** | 135 |
| 505 | -0.0023 | [-0.0290, +0.0308] | 0.0154 | **0.423** | 134 |

- P(b > 0): range **0.423 to 0.822**, SD across seeds 0.146.
- mean b: range -0.0023 to +0.0109, SD 0.0050.
- posterior SD: range 0.0135 to 0.0159.

## Verdict

**P(b > 0) is not stable to the precision it is quoted at.** Across seeds it spans 0.423 to 0.822, a range of 0.398, against a gap of 0.110 between the supplement's 0.70 and the median of these runs. The supplement quotes a seed-dependent summary to two decimal places.

**The remedy is not to pick a seed.** P(b > 0) is a tail mass evaluated where this posterior is densest, so it is the least stable summary available and the one most sensitive to the particle sample. The mean and the interval, which are integrals over the whole posterior, are comparatively steady and reproduce the supplement closely. Report the interval and drop the point probability, or report the probability with its across-seed range attached.

## The dependency for item G

Item G draws the population size from this posterior. Its median across seeds runs 87 to 176, a spread of 67 percent of the mean. The neutrality posterior predictive check must therefore either pool draws across seeds or report its envelope's sensitivity to the seed, rather than conditioning on one ABC run as if it were the posterior.
