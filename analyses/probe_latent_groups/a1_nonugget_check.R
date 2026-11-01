#!/usr/bin/env Rscript
# A1 without the nugget (use_tau = 0): groups as the only non-spatial term.
# Step 6 follow-up (2026-10-17), testing the diagnosed cause of the
# non-identification (the per-assemblage latent term makes group shifts
# redundant). On five own-generator records (all generated WITH a nugget of
# 0.47): zero shift (false-detection risk when scatter has nowhere else to
# go), compact sigma_g 0.5 and 2, scattered sigma_g 0.5 and 1. A0 (with
# nugget) on each for the in-sample comparison. Everything through the layer
# and cache. Output: output/probe_latent_groups/a1_nonugget_check.json.
suppressPackageStartupMessages({ library(here); library(rstan); library(posterior); library(jsonlite) })
here::i_am("pyproject.toml"); source(here("R", "fit_spatial.R")); options(mc.cores = 4L)
A0 <- here("stan", "spatial", "composition_gp_marginal.stan"); A1 <- here("stan", "spatial", "composition_gp_groups.stan")
sm0 <- compile_model(A0); sm1 <- compile_model(A1)
A1_KEYS <- c("G", "prior_group_sd", "fix_group_scale", "group_scale_fixed", "M", "edges")
GATED0 <- c("sigma", "tau", "rho", "alpha", "lp__", "f", "log_lik")
GATED1 <- c("sigma", "rho", "alpha", "lp__", "f", "log_lik", "sigma_g", "Delta", "Delta_w", "coherence", "comembership")
EXEMPT <- c("delta", "pi", "resp", "z_g", "sigma_g_free", "z_draw", "coherence_hard", "y_rep", "tau")
read_a1 <- function(f) { d <- read_stan_data(f); d$edges <- matrix(as.integer(unlist(d$edges)), nrow = d$M, byrow = FALSE); d }
strip_a1 <- function(d) { p <- attr(d, "provenance"); pa <- attr(d, "path"); d <- d[setdiff(names(d), A1_KEYS)]; attr(d, "provenance") <- p; attr(d, "path") <- pa; d }
summ <- function(fit, pars) { dd <- posterior::subset_draws(posterior::as_draws_array(fit), variable = pars)
  as.data.frame(posterior::summarise_draws(dd, "median", ~posterior::quantile2(.x, probs = c(0.05, 0.95)))) }
margins <- function(fit, z) { q <- apply(rstan::extract(fit, pars = "comembership")$comembership, c(2, 3), mean)
  g <- sort(unique(z)); B <- matrix(NA_real_, length(g), length(g))
  for (a in seq_along(g)) for (b in seq_along(g)) { ia <- which(z == g[a]); ib <- which(z == g[b]); B[a, b] <- if (a == b) mean(q[ia, ia][upper.tri(q[ia, ia])]) else mean(q[ia, ib]) }
  list(margin = sapply(seq_along(g), function(a) B[a, a] - max(B[a, -a])), block = B) }
ll <- function(fit) { x <- rstan::extract(fit, pars = "log_lik")$log_lik; mean(rowSums(x)) }
out <- list()
for (tag in c("zero", "sg0.5_r1", "sg2_r1", "scattered_sg0.5", "scattered_sg1")) {
  d <- read_a1(here("data", "stan", paste0("a1_", tag, ".json"))); tr <- attr(d, "provenance")$`_truth`
  r0 <- fit_cached(paste0("nonug_a0_", tag), A0, strip_a1(d), gated = GATED0, chains = 4L, seed = 20261017L, sm = sm0)
  d1 <- d; d1$use_tau <- 0L
  r1 <- fit_cached(paste0("nonug_a1_", tag), A1, d1, gated = GATED1, exempt = EXEMPT, constants = "comembership", comembership = "comembership",
                   chains = 4L, seed = 20261017L, sm = sm1)
  cat(diag_line(r1$diag, paste("A1 no nugget,", tag)), "\n")
  m <- margins(r1$fit, unlist(tr$z))
  out[[tag]] <- list(truth = tr[intersect(names(tr), c("sigma_g", "Delta_true", "realized_fst_between_true_groups", "coherence_true_hard"))],
                     A0_log_lik = ll(r0$fit), A0_gate = r0$gate$pass, A1_log_lik = ll(r1$fit),
                     post = summ(r1$fit, c("sigma", "rho", "sigma_g", "Delta", "Delta_w", "coherence")),
                     margins = m$margin, block = m$block, gate = r1$gate, split = r1$diag$comembership_split$max_abs_diff, rung = r1$rung, seconds = r1$seconds)
  cat(sprintf("  %s: A0 ll %.1f | A1(no nugget) ll %.1f | margins %s | Delta %.2f | coherence %.3f | split %.2f\n", tag, out[[tag]]$A0_log_lik, out[[tag]]$A1_log_lik,
              paste(sprintf("%.2f", m$margin), collapse = " "), out[[tag]]$post$median[4], out[[tag]]$post$median[6], out[[tag]]$split))
}
write_json(out, here("output", "probe_latent_groups", "a1_nonugget_check.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat("done\n")
