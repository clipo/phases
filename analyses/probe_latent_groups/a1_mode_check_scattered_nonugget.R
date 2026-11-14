#!/usr/bin/env Rscript
# Does the posterior of A1 on a large-signal own-generator record have a
# "groups" mode at all, or does the long-range field absorb the planted
# area-level shifts? Pilot step 6 diagnostic (rule 17), 2026-10-17.
# Variant (2026-10-17): the no-nugget A1 (use_tau = 0) on a1_scattered_sg1, started at the truth, to tell a sampler failure from an identification failure for the no-nugget form.
#   (a) A0 alone, to compare in-sample pointwise log-likelihood with A1;
#   (b) A1 started at the truth (z_g at the generating contrasts, sigma_g at 2,
#       pi at the true sizes), to see whether chains stay in a groups mode and
#       what lp__ and log_lik they reach, against the default-start fit.
# Run: Rscript analyses/probe_latent_groups/a1_mode_check.R
suppressPackageStartupMessages({ library(here); library(rstan); library(posterior); library(jsonlite) })
here::i_am("pyproject.toml"); source(here("R", "fit_spatial.R")); options(mc.cores = 4L)
A0 <- here("stan", "spatial", "composition_gp_marginal.stan"); A1 <- here("stan", "spatial", "composition_gp_groups.stan")
sm0 <- compile_model(A0); sm1 <- compile_model(A1)
f <- here("data", "stan", "a1_scattered_sg1.json"); d <- read_stan_data(f); d$edges <- matrix(as.integer(unlist(d$edges)), nrow = d$M); d$use_tau <- 0L
tr <- attr(d, "provenance")$`_truth`
A1_KEYS <- c("G", "prior_group_sd", "fix_group_scale", "group_scale_fixed", "M", "edges")
d0 <- d[setdiff(names(d), A1_KEYS)]; attr(d0, "provenance") <- attr(d, "provenance"); attr(d0, "path") <- attr(d, "path")
GATED <- c("sigma", "tau", "rho", "alpha", "lp__", "f", "log_lik")
r0 <- fit_cached("modecheck_a0_scattered_sg1_nonug", A0, d0, gated = GATED, chains = 4L, seed = 20261017L, sm = sm0)
cat(diag_line(r0$diag, "A0 on sg2_r1"), "\n")
# Truth-start inits: delta_true = sigma_g Q z_g  ->  z_g = Q' delta_true / sigma_g (Q orthonormal, delta sums to zero).
G <- 3; Q <- matrix(0, G, G - 1)
for (j in 1:(G - 1)) for (g in 1:G) Q[g, j] <- if (g <= j) 1 / sqrt(j * (j + 1)) else if (g == j + 1) -j / sqrt(j * (j + 1)) else 0
delta_true <- matrix(unlist(tr$delta), nrow = G, byrow = TRUE)
z_g0 <- t(Q) %*% delta_true / tr$sigma_g
pi0 <- as.numeric(unlist(tr$pi_true))
init_fun <- function(chain_id) list(sigma_g_free = array(tr$sigma_g, 1), z_g = z_g0, pi = pi0,
                                    sigma = tr$sigma, rho = tr$rho, tau = numeric(0))
cfg <- list(iter = 3000L, warmup = 1500L, adapt_delta = 0.95, max_treedepth = 12L)
r1t <- fit_once(sm1, d, cfg, chains = 4L, seed = 20261017L, init = init_fun)
d1t <- gp_diagnostics(r1t$fit, gated = c(GATED, "sigma_g", "Delta", "Delta_w", "coherence", "comembership"),
                      exempt = c("delta", "pi", "resp", "z_g", "sigma_g_free", "z_draw", "coherence_hard", "y_rep"),
                      constants = "comembership", comembership = "comembership")
cat(diag_line(d1t, "A1 truth-start on sg2_r1"), " split", round(d1t$comembership_split$max_abs_diff, 3), "\n")
# Compare with the cached default-start A1 fit.
rds <- list.files(here("output", "gp_cache"), pattern = "^nonug_a1_scattered_sg1_.*\\.rds$", full.names = TRUE)
r1 <- readRDS(rds[1])
summ <- function(fit, pars) { dd <- posterior::subset_draws(posterior::as_draws_array(fit), variable = pars)
  as.data.frame(posterior::summarise_draws(dd, "median", ~posterior::quantile2(.x, probs = c(0.05, 0.95)))) }
ll <- function(fit) { x <- rstan::extract(fit, pars = "log_lik")$log_lik; c(mean_total = mean(rowSums(x)), sd_total = sd(rowSums(x))) }
margins <- function(fit, z) { q <- apply(rstan::extract(fit, pars = "comembership")$comembership, c(2, 3), mean)
  g <- sort(unique(z)); B <- matrix(NA_real_, length(g), length(g))
  for (a in seq_along(g)) for (b in seq_along(g)) { ia <- which(z == g[a]); ib <- which(z == g[b]); B[a, b] <- if (a == b) mean(q[ia, ia][upper.tri(q[ia, ia])]) else mean(q[ia, ib]) }
  sapply(seq_along(g), function(a) B[a, a] - max(B[a, -a])) }
per_chain_lp <- function(fit) { lp <- rstan::extract(fit, pars = "lp__", permuted = FALSE)[, , 1]; round(colMeans(lp), 1) }
out <- list(
  truth = list(sigma_g = tr$sigma_g, Delta_true = tr$Delta_true, realized_fst = tr$realized_fst_between_true_groups, sigma = tr$sigma, tau = tr$tau, rho = tr$rho),
  A0 = list(log_lik = ll(r0$fit), post = summ(r0$fit, c("sigma", "rho")), gate = r0$gate$pass),
  A1_default_start = list(log_lik = ll(r1$fit), post = summ(r1$fit, c("sigma", "rho", "sigma_g", "Delta", "coherence")), margins = margins(r1$fit, unlist(tr$z)), lp_per_chain = per_chain_lp(r1$fit), gate = r1$gate$pass, split = r1$diag$comembership_split$max_abs_diff),
  A1_truth_start = list(log_lik = ll(r1t$fit), post = summ(r1t$fit, c("sigma", "rho", "sigma_g", "Delta", "coherence")), margins = margins(r1t$fit, unlist(tr$z)), lp_per_chain = per_chain_lp(r1t$fit), gate = gate_fit(d1t)$pass, reasons = gate_fit(d1t)$reasons, split = d1t$comembership_split$max_abs_diff, seconds = r1t$seconds))
write_json(out, here("output", "probe_latent_groups", "a1_mode_check_scattered_nonugget.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
print(out[c("truth")]); cat("\nA0 log_lik total:", out$A0$log_lik, "\nA1 default:", out$A1_default_start$log_lik, " margins", round(out$A1_default_start$margins, 2), " lp per chain", out$A1_default_start$lp_per_chain,
    "\nA1 truth-start:", out$A1_truth_start$log_lik, " margins", round(out$A1_truth_start$margins, 2), " lp per chain", out$A1_truth_start$lp_per_chain, "\n")
print(out$A1_truth_start$post); print(out$A1_default_start$post)
