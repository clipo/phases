#!/usr/bin/env Rscript
# Mechanism of the non-identification found in step 6 (rule 17), 2026-10-17.
#   (a) A1 on the sigma_g = 2 compact-group record with the field's range pinned
#       to 12 km (rho_lo 11.99, rho_hi 12.01, uniform prior): if the free range
#       is what lets the field absorb area-level steps, the groups reappear.
#   (b) A1, default priors, on records whose groups of the same sizes are
#       SCATTERED over the sites (sigma_g 0.5 and 1.0): groups the smooth field
#       cannot mimic should be recovered if A1 works as a group detector.
# Every fit through the layer and cache. Output: output/probe_latent_groups/a1_mechanism_check.json.
suppressPackageStartupMessages({ library(here); library(rstan); library(posterior); library(jsonlite) })
here::i_am("pyproject.toml"); source(here("R", "fit_spatial.R")); options(mc.cores = 4L)
A1 <- here("stan", "spatial", "composition_gp_groups.stan"); sm1 <- compile_model(A1)
GATED <- c("sigma", "tau", "rho", "alpha", "lp__", "f", "log_lik", "sigma_g", "Delta", "Delta_w", "coherence", "comembership")
EXEMPT <- c("delta", "pi", "resp", "z_g", "sigma_g_free", "z_draw", "coherence_hard", "y_rep")
read_a1 <- function(f) { d <- read_stan_data(f); d$edges <- matrix(as.integer(unlist(d$edges)), nrow = d$M, byrow = FALSE); d }
summ <- function(fit, pars) { dd <- posterior::subset_draws(posterior::as_draws_array(fit), variable = pars)
  as.data.frame(posterior::summarise_draws(dd, "median", ~posterior::quantile2(.x, probs = c(0.05, 0.95)))) }
margins <- function(fit, z) { q <- apply(rstan::extract(fit, pars = "comembership")$comembership, c(2, 3), mean)
  g <- sort(unique(z)); B <- matrix(NA_real_, length(g), length(g))
  for (a in seq_along(g)) for (b in seq_along(g)) { ia <- which(z == g[a]); ib <- which(z == g[b]); B[a, b] <- if (a == b) mean(q[ia, ia][upper.tri(q[ia, ia])]) else mean(q[ia, ib]) }
  list(margin = sapply(seq_along(g), function(a) B[a, a] - max(B[a, -a])), block = B) }
ll <- function(fit) { x <- rstan::extract(fit, pars = "log_lik")$log_lik; mean(rowSums(x)) }
out <- list()
# (a) pinned range on the compact sigma_g = 2 record
d <- read_a1(here("data", "stan", "a1_sg2_r1.json")); tr <- attr(d, "provenance")$`_truth`
dp <- d; dp$prior_rho_type <- 3L; dp$rho_lo <- 11.99; dp$rho_hi <- 12.01
r <- fit_cached("mech_a1_rho12_sg2_r1", A1, dp, gated = setdiff(GATED, "rho"), exempt = c(EXEMPT, "rho"), constants = "comembership",
                comembership = "comembership", chains = 4L, seed = 20261017L, sm = sm1)
cat(diag_line(r$diag, "A1 rho pinned 12 km, sg2_r1"), "\n")
m <- margins(r$fit, unlist(tr$z))
out$rho12_compact_sg2 <- list(truth = tr[c("sigma_g", "Delta_true", "realized_fst_between_true_groups")], post = summ(r$fit, c("sigma", "tau[1]", "sigma_g", "Delta", "Delta_w", "coherence")),
                             margins = m$margin, block = m$block, log_lik = ll(r$fit), gate = r$gate, split = r$diag$comembership_split$max_abs_diff, seconds = r$seconds, rung = r$rung)
cat("  margins", round(m$margin, 2), " Delta", round(out$rho12_compact_sg2$post$median[4], 2), " log_lik", round(out$rho12_compact_sg2$log_lik, 1), "\n")
# (b) scattered partitions, default priors
for (tag in c("scattered_sg0.5", "scattered_sg1")) {
  d <- read_a1(here("data", "stan", paste0("a1_", tag, ".json"))); tr <- attr(d, "provenance")$`_truth`
  r <- fit_cached(paste0("mech_a1_", tag), A1, d, gated = GATED, exempt = EXEMPT, constants = "comembership", comembership = "comembership",
                  chains = 4L, seed = 20261017L, sm = sm1)
  cat(diag_line(r$diag, tag), "\n")
  m <- margins(r$fit, unlist(tr$z))
  out[[tag]] <- list(truth = tr[c("sigma_g", "Delta_true", "coherence_true_hard", "realized_fst_between_true_groups")],
                     post = summ(r$fit, c("sigma", "tau[1]", "rho", "sigma_g", "Delta", "Delta_w", "coherence")),
                     margins = m$margin, block = m$block, log_lik = ll(r$fit), gate = r$gate, split = r$diag$comembership_split$max_abs_diff, seconds = r$seconds, rung = r$rung)
  cat("  margins", round(m$margin, 2), " Delta", round(out[[tag]]$post$median[5], 2), " coherence", round(out[[tag]]$post$median[7], 3), "\n")
}
write_json(out, here("output", "probe_latent_groups", "a1_mechanism_check.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat("done\n")
