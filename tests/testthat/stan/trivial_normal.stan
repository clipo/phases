// Test fixture: a two-parameter normal model that samples cleanly.
data { int<lower=0> N; vector[N] y; }
parameters { real mu; real<lower=0> sigma; }
model { mu ~ normal(0, 5); sigma ~ normal(0, 2); y ~ normal(mu, sigma); }
