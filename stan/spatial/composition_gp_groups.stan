// A1: A0 (composition_gp_marginal.stan) plus G latent composition groups.
//
// Written 2026-10-17 for the latent-group probe (spec revision 2.1, section 4;
// D13, D15, D16, D18). Everything A0 has is kept unchanged: the same data
// contract, priors, kernel include, nugget-in-covariance field. Added:
//
//   G latent groups with membership z_i ~ Categorical(pi), pi ~ Dirichlet(1),
//   summed out of the likelihood (one group per assemblage, never a per-sherd
//   mixture); group shifts delta (G x C) that SUM TO ZERO across groups within
//   each class, delta[, c] = sigma_g * Q * z_g[, c] with Q an orthonormal basis
//   of the subspace orthogonal to the ones vector (Helmert contrasts), so the
//   intercept prior is identical to A0's and each shift's marginal prior SD is
//   sigma_g * sqrt(1 - 1/G) (D16, declared prior change relative to iid shifts);
//   sigma_g ~ half-Normal(0, prior_group_sd). With fix_group_scale = 1 the group
//   scale is replaced by group_scale_fixed (0 for the code-equivalence test
//   against A0, spec 4.1c).
//
// Generated quantities, all per draw and label-invariant unless noted:
//   resp[i, g]        P(z_i = g | parameters, y_i)
//   comembership      sum_g resp[i,g] resp[j,g] for i != j, diagonal 1 (D13)
//   Delta             sqrt(sum_c sum_g pi_g (delta_gc - mean_c)^2), mean_c the
//                     pi-weighted mean (D15); Delta_w the same with the mean
//                     responsibility per group as weight (R2-5)
//   coherence         mean comembership over the graph edges passed as data
//                     minus mean over all pairs (D14, D18; conditional-mean form)
//   coherence_hard    the same on one categorical draw of z per iteration
//   log_lik[i]        the mixture log density, for held-out scoring
//   z_draw, y_rep     one membership draw and posterior predictive counts
//   delta, pi         label-sensitive; reported, never gated (spec 4.2)
functions {
#include kernels.stanfunctions
}
data {
  int<lower=2> N;                       // assemblages
  int<lower=2> K;                       // decorated classes
  array[N, K] int<lower=0> y;           // counts
  matrix[N, N] D;                       // pairwise distance, km
  int<lower=1, upper=3> prior_rho_type; // 1 inv_gamma, 2 proportional to 1/rho, 3 uniform
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
  // A1 additions
  int<lower=2> G;                       // latent groups
  real<lower=0> prior_group_sd;         // half-normal sd on sigma_g
  int<lower=0, upper=1> fix_group_scale;
  real<lower=0> group_scale_fixed;
  int<lower=0> M;                       // neighbor-graph edges for coherence
  array[M, 2] int<lower=1, upper=N> edges;
}
transformed data {
  real jitter = 1e-8;
  real rc = range_const(kernel);
  real hc = half_const(kernel);
  int C = K - 1;
  array[N] int n_i;
  int n_pairs = N * (N - 1) / 2;
  matrix[G, G - 1] Q;                   // orthonormal contrasts: columns sum to zero, Q'Q = I
  for (i in 1:N) n_i[i] = sum(y[i]);
  for (j in 1:(G - 1)) {
    for (g in 1:G) {
      if (g <= j) Q[g, j] = 1.0 / sqrt(j * (j + 1.0));
      else if (g == j + 1) Q[g, j] = -j / sqrt(j * (j + 1.0));
      else Q[g, j] = 0;
    }
  }
}
parameters {
  vector[C] alpha;
  vector<lower=0>[use_tau] tau;
  real<lower=0> sigma;
  real<lower=rho_lo, upper=rho_hi> rho;
  matrix[N, C] z_f;
  vector<lower=0>[1 - fix_group_scale] sigma_g_free;
  matrix[G - 1, C] z_g;                 // non-centered contrasts
  simplex[G] pi;
}
transformed parameters {
  matrix[N, C] f;
  real sigma_g;
  matrix[G, C] delta;
  {
    matrix[N, N] Sigma = k_matern(D, sigma, rho, N, jitter, kernel);
    if (use_tau) {
      for (i in 1:N) Sigma[i, i] += square(tau[1]);
    }
    f = cholesky_decompose(Sigma) * z_f;
  }
  sigma_g = fix_group_scale ? group_scale_fixed : sigma_g_free[1];
  delta = sigma_g * (Q * z_g);
}
model {
  if (use_tau) tau ~ normal(0, prior_scale_sd);
  sigma ~ normal(0, prior_scale_sd);
  if (prior_rho_type == 1) rho ~ inv_gamma(prior_rho_a, prior_rho_b);
  else if (prior_rho_type == 2) target += -log(rho);
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
  to_vector(z_f) ~ std_normal();
  if (!fix_group_scale) sigma_g_free ~ normal(0, prior_group_sd);
  to_vector(z_g) ~ std_normal();
  pi ~ dirichlet(rep_vector(1, G));
  if (likelihood) {
    for (i in 1:N) {
      vector[G] lp;
      for (g in 1:G) {
        vector[K] eta;
        eta[1] = 0;
        for (k in 2:K) eta[k] = alpha[k - 1] + f[i, k - 1] + delta[g, k - 1];
        lp[g] = log(pi[g]) + multinomial_logit_lpmf(y[i] | eta);
      }
      target += log_sum_exp(lp);
    }
  }
}
generated quantities {
  real practical_range = rc * rho;
  real half_distance = hc * rho;
  real spatial_share = use_tau ? square(sigma) / (square(sigma) + square(tau[1])) : 1;
  matrix[N, G] resp;
  vector[N] log_lik;
  matrix[N, N] comembership;
  real Delta;
  real Delta_w;
  real coherence;
  real coherence_hard;
  array[N] int z_draw;
  array[N, K] int y_rep;
  {
    for (i in 1:N) {
      vector[G] lp;
      for (g in 1:G) {
        vector[K] eta;
        eta[1] = 0;
        for (k in 2:K) eta[k] = alpha[k - 1] + f[i, k - 1] + delta[g, k - 1];
        lp[g] = log(pi[g]) + multinomial_logit_lpmf(y[i] | eta);
      }
      log_lik[i] = log_sum_exp(lp);
      resp[i] = softmax(lp)';
      z_draw[i] = categorical_rng(softmax(lp));
      {
        vector[K] eta;
        eta[1] = 0;
        for (k in 2:K) eta[k] = alpha[k - 1] + f[i, k - 1] + delta[z_draw[i], k - 1];
        y_rep[i] = multinomial_logit_rng(eta, n_i[i]);
      }
    }
    for (i in 1:N) {
      comembership[i, i] = 1;
      for (j in (i + 1):N) {
        real q = dot_product(resp[i], resp[j]);
        comembership[i, j] = q;
        comembership[j, i] = q;
      }
    }
    {
      real d2 = 0;
      real d2w = 0;
      vector[G] w;
      for (g in 1:G) w[g] = mean(col(resp, g));
      for (c in 1:C) {
        real m = dot_product(pi, col(delta, c));
        real mw = dot_product(w, col(delta, c));
        for (g in 1:G) {
          d2 += pi[g] * square(delta[g, c] - m);
          d2w += w[g] * square(delta[g, c] - mw);
        }
      }
      Delta = sqrt(d2);
      Delta_w = sqrt(d2w);
    }
    {
      real qe = 0;
      real qa = 0;
      real he = 0;
      real ha = 0;
      for (e in 1:M) {
        qe += comembership[edges[e, 1], edges[e, 2]];
        he += (z_draw[edges[e, 1]] == z_draw[edges[e, 2]]);
      }
      for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
          qa += comembership[i, j];
          ha += (z_draw[i] == z_draw[j]);
        }
      }
      coherence = (M > 0 ? qe / M : 0) - qa / n_pairs;
      coherence_hard = (M > 0 ? he / M : 0) - ha / n_pairs;
    }
  }
}
