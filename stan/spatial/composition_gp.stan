// Decorated-ceramic composition as a latent spatial field.
//
// THE QUESTION. The manuscript asks whether the culture-historical phases are
// real units or a continuum. Every cultural F_ST it reports is measured between
// spatial clusters produced by k-means with k chosen by a silhouette score
// (docs/CODE_REVIEW_2026-08-31.md F15), so "between-group" is a partition we
// drew. This model removes the partition. It asks instead how assemblage
// composition varies over space, with two competing sources:
//
//   tau    an assemblage-level effect with NO spatial structure. Large tau means
//          each site has its own character regardless of its neighbours.
//   sigma  the amplitude of a SPATIALLY CORRELATED effect with length scale rho.
//          Large sigma with a resolved rho means nearby assemblages resemble
//          each other because they are near, not because they share a label.
//
// `spatial_share` = sigma^2 / (sigma^2 + tau^2) is the quantity the paper's
// argument turns on. Near 1, composition is organised by geography and smooth
// distance-decay suffices, so bounded phase units are unnecessary. Near 0,
// assemblages differ in ways geography does not explain, which is where a
// bounded-group account would have to live.
//
// STRUCTURE adapted from ../mataa stan/spatial/moai_site_gp.stan (2026-08-31).
// The change is the likelihood: mataa observes one categorical state per object,
// while an assemblage here is a vector of counts over K decorated classes, so
// the per-object categorical becomes one multinomial per assemblage. Class 1 is
// the reference category (eta = 0), which fixes the softmax translation
// redundancy.
//
// INHERITED WARNING, to be tested here rather than assumed away (rule 19).
// mataa's header records that on process-simulated moai data "the tau + field
// decomposition returned d50 of 7 to 13 km whatever the truth, the nugget-range
// trade-off". A nugget and a short length scale can mimic each other. The
// `use_tau` switch exists so that failure mode is measurable: fit with and
// without the non-spatial term and compare what rho does. If rho is pinned by
// the prior rather than the data, that must be found before any number is
// reported, not after.
functions {
#include kernels.stanfunctions
}
data {
  int<lower=2> N;                       // assemblages
  int<lower=2> K;                       // decorated classes
  array[N, K] int<lower=0> y;           // counts
  matrix[N, N] D;                       // pairwise distance, km

  // rho prior. 1 = inv_gamma(a, b); 2 = uniform on log(rho) in [rho_lo, rho_hi].
  // Both are offered because the choice is contested and rule 20 requires the
  // sensitivity refit to lean the other way; neither is hard-coded.
  int<lower=1, upper=2> prior_rho_type;
  real<lower=0> prior_rho_a;
  real<lower=0> prior_rho_b;
  real<lower=0> rho_lo;
  real<lower=rho_lo> rho_hi;

  int<lower=1, upper=3> kernel;         // 1 exp, 2 Matern 3/2, 3 Matern 5/2
  int<lower=0, upper=1> use_tau;        // 0 removes the non-spatial term
  real<lower=0> prior_scale_sd;         // half-normal sd, tau and sigma alike
  int<lower=0, upper=1> likelihood;     // 0 = prior predictive
}
transformed data {
  real jitter = 1e-8;
  real rc = range_const(kernel);
  real hc = half_const(kernel);
  int C = K - 1;                        // free categories; class 1 is reference
  array[N] int n_i;
  for (i in 1:N) n_i[i] = sum(y[i]);
}
parameters {
  vector[C] alpha;                      // class intercepts, relative to class 1
  vector<lower=0>[use_tau] tau;
  real<lower=0> sigma;
  real<lower=rho_lo, upper=rho_hi> rho;
  matrix[N, C] z_u;                     // non-centred iid assemblage effects
  matrix[N, C] z_f;                     // non-centred GP effects
}
transformed parameters {
  matrix[N, C] u = use_tau ? tau[1] * z_u : rep_matrix(0, N, C);
  matrix[N, C] f;
  {
    matrix[N, N] L = cholesky_decompose(k_matern(D, sigma, rho, N, jitter, kernel));
    f = L * z_f;
  }
}
model {
  // tau and sigma take the SAME prior deliberately. Giving the spatial term a
  // tighter or looser prior than the non-spatial one would tilt the very
  // comparison this model exists to make.
  if (use_tau) tau ~ normal(0, prior_scale_sd);
  sigma ~ normal(0, prior_scale_sd);
  if (prior_rho_type == 1) rho ~ inv_gamma(prior_rho_a, prior_rho_b);
  else target += -log(rho);             // p(rho) proportional to 1/rho
  alpha ~ normal(0, 2);
  to_vector(z_u) ~ std_normal();
  to_vector(z_f) ~ std_normal();

  if (likelihood) {
    for (i in 1:N) {
      vector[K] eta;
      eta[1] = 0;
      for (k in 2:K) eta[k] = alpha[k - 1] + u[i, k - 1] + f[i, k - 1];
      y[i] ~ multinomial_logit(eta);
    }
  }
}
generated quantities {
  real practical_range = rc * rho;      // correlation 0.05
  real half_distance = hc * rho;        // correlation 0.5; the primary scale
  real spatial_share = use_tau ? square(sigma) / (square(sigma) + square(tau[1])) : 1;
  array[N, K] int y_rep;
  {
    for (i in 1:N) {
      vector[K] eta;
      eta[1] = 0;
      for (k in 2:K) eta[k] = alpha[k - 1] + u[i, k - 1] + f[i, k - 1];
      y_rep[i] = multinomial_logit_rng(eta, n_i[i]);
    }
  }
}
