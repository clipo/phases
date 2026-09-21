# How reproducible is the ABC transmission posterior?

Produced by `analyses/55_abc_stability.py` (full: 800 particles, 6 rounds, 5 seeds). Chases F22 and settles the dependency for rule-18 item G, which draws population sizes from this posterior.

## Across-seed spread

| seed | mean b | 95% interval | posterior SD | P(b > 0) | N median |
|---|---|---|---|---|---|
| 101 | +0.0008 | [-0.0259, +0.0294] | 0.0144 | **0.517** | 123 |
| 202 | +0.0072 | [-0.0191, +0.0402] | 0.0142 | **0.698** | 131 |
| 303 | +0.0100 | [-0.0122, +0.0420] | 0.0134 | **0.781** | 106 |
| 404 | +0.0142 | [-0.0109, +0.0430] | 0.0138 | **0.886** | 129 |
| 505 | +0.0126 | [-0.0138, +0.0459] | 0.0145 | **0.817** | 114 |

- P(b > 0): range **0.517 to 0.886**, SD across seeds 0.127.
- mean b: range +0.0008 to +0.0142, SD 0.0047.
- posterior SD: range 0.0134 to 0.0145.

## Verdict

**P(b > 0) is not stable to the precision it is quoted at.** Across seeds it spans 0.517 to 0.886, a range of 0.368, against a gap of 0.081 between the supplement's 0.70 and the median of these runs. The supplement quotes a seed-dependent summary to two decimal places.

**The remedy is not to pick a seed.** P(b > 0) is a tail mass evaluated where this posterior is densest, so it is the least stable summary available and the one most sensitive to the particle sample. The mean and the interval, which are integrals over the whole posterior, are comparatively steady and reproduce the supplement closely. Report the interval and drop the point probability, or report the probability with its across-seed range attached.

## The dependency for item G

Item G draws the population size from this posterior. Its median across seeds runs 106 to 131, a spread of 21 percent of the mean. The neutrality posterior predictive check must therefore either pool draws across seeds or report its envelope's sensitivity to the seed, rather than conditioning on one ABC run as if it were the posterior.
