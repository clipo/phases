# Hierarchical Bayesian convergence posterior

Four emergence slopes (3 CA-panel signatures + seriation fragmentation), partial pooling, measurement-error likelihood. full config: 2000 draws x 4 chains.

## Convergence

- **P(all four slopes > 0) = 0.000**.
- Convergence-score slope mean -0.512, P(>0) = 0.008.
- Shared emergence slope mu: mean -0.498, P(>0) = 0.046.

## Per-signature slopes (emergence-oriented, standardized)

| signature | posterior mean | 95% interval | P(>0) |
|---|---|---|---|
| neutral_departure | -0.618 | [-1.033, -0.244] | 0.001 |
| fst | -0.448 | [-0.786, -0.015] | 0.022 |
| spatial_boundary | -0.348 | [-1.113, +0.800] | 0.176 |
| seriation | -0.633 | [-0.739, -0.525] | 0.000 |

## MCMC diagnostics

- max R-hat = 1.0018 (want < 1.01); min ESS = 1979; divergences = 0.
