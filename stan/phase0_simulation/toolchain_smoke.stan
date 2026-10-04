// Minimal model proving the C++ toolchain compiles and NUTS runs.
//
// Deliberately conjugate so the target is computed in R from the prior and the
// likelihood, never read back off this fit (Verification Regime rule 3).
// With prior mu ~ Normal(0, 10) and y_i ~ Normal(mu, 1) on N points,
//   posterior precision = 1/100 + N,  posterior mean = (N * ybar) / (1/100 + N).
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  real mu;
}
model {
  mu ~ normal(0, 10);
  y ~ normal(mu, 1);
}
