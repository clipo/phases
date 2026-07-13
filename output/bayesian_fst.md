# Bayesian credible interval on the basin cultural F_ST

Observed St. Francis basin, 5 spatial clusters, 10 decorated types. Hierarchical Dirichlet-multinomial, 3000 draws x 4 chains (full).

## Cultural F_ST (Gini-Simpson estimator)

- Posterior mean F_ST = 0.0327, 95% credible interval [0.0302, 0.0353].
- Frequentist plug-in (variance.cultural_fst) = 0.0330.
- P(F_ST > 0) = 1.000.

The credible interval is the estimation uncertainty on the reported F_ST. Whether the value is distinguishable from stochastic drift is a separate test (analyses 21/33/35/37).

## MCMC diagnostics

- max R-hat = 1.0086 (want < 1.01); min ESS = 835; divergences = 0.
