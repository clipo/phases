# Hierarchical Bayesian convergence posterior

Four emergence slopes (3 CA-panel signatures + seriation fragmentation), partial pooling, measurement-error likelihood. full config: 2000 draws x 4 chains.

## Convergence

- **P(all four slopes > 0) = 0.000**.
- Convergence-score slope mean -0.506, P(>0) = 0.008.
- Shared emergence slope mu: mean -0.485, P(>0) = 0.050.

## Per-signature slopes (emergence-oriented, standardized)

| signature | posterior mean | 95% interval | P(>0) |
|---|---|---|---|
| neutral_departure | -0.620 | [-1.009, -0.225] | 0.001 |
| fst | -0.437 | [-0.802, -0.021] | 0.028 |
| spatial_boundary | -0.334 | [-1.156, +0.812] | 0.192 |
| seriation | -0.633 | [-0.743, -0.525] | 0.000 |

## MCMC diagnostics

- max R-hat = 1.0025 (want < 1.01); min ESS = 2025; divergences = 0.
