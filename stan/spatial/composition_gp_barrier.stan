// A3: A0 (composition_gp_marginal.stan) plus a region-level shift variance
// omega^2 on same-region pairs, with the placement a discrete latent variable
// over P candidate region assignments (design: docs/working/
// 2026-10-17-a3-barrier-covariance-spec.md, D23, 2026-10-17).
//
//   K_m = sigma^2 k(D; rho) + omega^2 S_m + tau^2 I,   S_m[i,j] = 1[region_m[i] == region_m[j]]
//   p(y | theta) = sum_m (1/P) prod_i Multinomial(y_i | softmax(0, alpha + f_m[i]))
//   f_m = chol(K_m) z_f  (one standard-normal matrix shared across placements)
//
// omega is the SD in logits of a shift shared by every site of a region, per
// class, beyond the distance-decay field. fix_omega = 1 replaces omega by 0,
// which is A0 exactly (code-equivalence test). Generated quantities: the
// placement responsibilities resp[m], one categorical draw m_draw, the field
// under m_draw, pointwise log_lik under the responsibilities, y_rep.
functions {
#include kernels.stanfunctions
}
data {
  int<lower=2> N;
  int<lower=2> K;
  array[N, K] int<lower=0> y;
  matrix[N, N] D;
  int<lower=1, upper=3> prior_rho_type;
  real<lower=0> prior_rho_a;
  real<lower=0> prior_rho_b;
  real<lower=0> rho_lo;
  real<lower=rho_lo> rho_hi;
  int<lower=1, upper=3> kernel;
  int<lower=0, upper=1> use_tau;
  real<lower=0> prior_scale_sd;
  // Intercept prior as data (2026-10-18, D25): the hard-coded Normal(0, 2) was too narrow for the
  // log-odds of rare classes against the reference, so the field's constant mode carried 0.3 to 5.6
  // logits of the class baselines (discussion log, 2026-10-18).
  vector[K - 1] prior_alpha_mean;
  real<lower=0> prior_alpha_sd;
  int<lower=0, upper=1> likelihood;
  // A3 additions
  int<lower=1> P;                          // candidate placements
  array[P, N] int<lower=1> region;         // region label of each site under each placement
  real<lower=0> prior_group_sd;            // half-normal sd on omega
  int<lower=0, upper=1> fix_omega;         // 1: omega = 0 (A0 exactly)
}
transformed data {
  real jitter = 1e-8;
  real rc = range_const(kernel);
  real hc = half_const(kernel);
  int C = K - 1;
  array[N] int n_i;
  array[P] matrix[N, N] S;
  for (i in 1:N) n_i[i] = sum(y[i]);
  for (m in 1:P) {
    for (i in 1:N) for (j in 1:N) S[m][i, j] = (region[m, i] == region[m, j]) ? 1 : 0;
  }
}
parameters {
  vector[C] alpha;
  vector<lower=0>[use_tau] tau;
  real<lower=0> sigma;
  real<lower=rho_lo, upper=rho_hi> rho;
  matrix[N, C] z_f;
  vector<lower=0>[1 - fix_omega] omega_free;
}
transformed parameters {
  real omega = fix_omega ? 0 : omega_free[1];
  array[P] matrix[N, C] f;
  vector[P] lp_m;
  {
    matrix[N, N] base = k_matern(D, sigma, rho, N, jitter, kernel);
    if (use_tau) {
      for (i in 1:N) base[i, i] += square(tau[1]);
    }
    for (m in 1:P) {
      f[m] = cholesky_decompose(base + square(omega) * S[m]) * z_f;
      lp_m[m] = -log(P);
      if (likelihood) {
        for (i in 1:N) {
          vector[K] eta;
          eta[1] = 0;
          for (k in 2:K) eta[k] = alpha[k - 1] + f[m][i, k - 1];
          lp_m[m] += multinomial_logit_lpmf(y[i] | eta);
        }
      }
    }
  }
}
model {
  if (use_tau) tau ~ normal(0, prior_scale_sd);
  sigma ~ normal(0, prior_scale_sd);
  if (prior_rho_type == 1) rho ~ inv_gamma(prior_rho_a, prior_rho_b);
  else if (prior_rho_type == 2) target += -log(rho);
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
  to_vector(z_f) ~ std_normal();
  if (!fix_omega) omega_free ~ normal(0, prior_group_sd);
  target += log_sum_exp(lp_m);
}
generated quantities {
  real practical_range = rc * rho;
  real half_distance = hc * rho;
  real spatial_share = use_tau ? square(sigma) / (square(sigma) + square(tau[1])) : 1;
  vector[P] resp = softmax(lp_m);
  int m_draw = categorical_rng(resp);
  matrix[N, C] f_draw = f[m_draw];
  vector[N] log_lik;
  array[N, K] int y_rep;
  for (i in 1:N) {
    vector[P] li;
    for (m in 1:P) {
      vector[K] eta;
      eta[1] = 0;
      for (k in 2:K) eta[k] = alpha[k - 1] + f[m][i, k - 1];
      li[m] = log(resp[m]) + multinomial_logit_lpmf(y[i] | eta);
    }
    log_lik[i] = log_sum_exp(li);
    {
      vector[K] eta;
      eta[1] = 0;
      for (k in 2:K) eta[k] = alpha[k - 1] + f_draw[i, k - 1];
      y_rep[i] = multinomial_logit_rng(eta, n_i[i]);
    }
  }
}
