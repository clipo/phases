#!/usr/bin/env Rscript
# A3 identification test with the placement supplied (design section 4, parts
# 1, 1b and 2; D23). Step 6b. Exploratory design validation (rule 1).
# Fits A3 (P = 1, the record's own placement) to every a3_*.json record, plus
# the code-equivalence pair (A0 against A3 with omega fixed at 0) on one
# zero-omega record. Reads omega against the truth, the identified field
# quantities, and the gate. Everything through the layer and the cache.
# Run: Rscript analyses/probe_latent_groups/validate_a3.R
suppressPackageStartupMessages({ library(here); library(rstan); library(posterior); library(jsonlite) })
here::i_am("pyproject.toml"); source(here("R", "fit_spatial.R")); source(here("R", "diagnostic_plots.R")); options(mc.cores = 4L)
A0 <- here("stan", "spatial", "composition_gp_marginal.stan"); A3 <- here("stan", "spatial", "composition_gp_barrier.stan")
sm0 <- compile_model(A0); sm3 <- compile_model(A3)
OUTD <- here("output", "probe_latent_groups"); FIG <- here("figures", "probe_latent_groups"); FIND <- here("output", "findings", "a3_validation.md")
SEED <- 20261017L
A3_KEYS <- c("P", "region", "prior_group_sd", "fix_omega")
GATED0 <- c("sigma", "tau", "rho", "alpha", "lp__", "f", "log_lik")
GATED3 <- c(GATED0, "omega")
EXEMPT3 <- c("resp", "m_draw", "lp_m", "f_draw", "y_rep", "omega_free")
read_a3 <- function(f) { d <- read_stan_data(f); d$region <- matrix(as.integer(unlist(d$region)), nrow = d$P, byrow = TRUE); d }
strip_a3 <- function(d) { p <- attr(d, "provenance"); pa <- attr(d, "path"); d <- d[setdiff(names(d), A3_KEYS)]; attr(d, "provenance") <- p; attr(d, "path") <- pa; d }
post <- function(fit, pars) { dd <- posterior::subset_draws(posterior::as_draws_array(fit), variable = pars)
  as.data.frame(posterior::summarise_draws(dd, "mean", "median", "mcse_mean", ~posterior::quantile2(.x, probs = c(0.05, 0.95)))) }
ll <- function(fit) mean(rowSums(rstan::extract(fit, pars = "log_lik")$log_lik))
timing <- list(); rows <- list()

# Code equivalence on a zero-omega record.
dz <- read_a3(here("data", "stan", "a3_own_G1_om0_r1.json"))
r0 <- fit_cached("a3val_a0_on_om0", A0, strip_a3(dz), gated = GATED0, chains = 4L, seed = SEED, sm = sm0)
dzf <- dz; dzf$fix_omega <- 1L
r3f <- fit_cached("a3val_a3_fixed0_on_om0", A3, dzf, gated = GATED0, exempt = c(EXEMPT3, "omega"), chains = 4L, seed = SEED, sm = sm3)
shared <- c("sigma", "tau[1]", "rho", paste0("alpha[", 1:9, "]"), "spatial_share", paste0("log_lik[", 1:28, "]"))
p0 <- post(r0$fit, shared); p3 <- post(r3f$fit, shared)
zsc <- (p3$mean - p0$mean) / sqrt(p0$mcse_mean^2 + p3$mcse_mean^2)
cat(sprintf("equivalence: max |z| over %d shared quantities %.2f (%s)\n", length(zsc), max(abs(zsc)), p0$variable[which.max(abs(zsc))]))

files <- sort(list.files(here("data", "stan"), pattern = "^a3_(own|slope|sim)_.*\\.json$", full.names = TRUE))
exemplar_done <- character()
for (f in files) {
  tag <- sub("^a3_(.*)\\.json$", "\\1", basename(f)); d <- read_a3(f); tr <- attr(d, "provenance")$`_truth`
  r <- fit_cached(paste0("a3val_", tag), A3, d, gated = GATED3, exempt = EXEMPT3, chains = 4L, seed = SEED, sm = sm3)
  timing[[tag]] <- list(rung = r$rung, seconds = r$seconds, pass = r$gate$pass, reasons = r$gate$reasons)
  cat(diag_line(r$diag, tag), "\n")
  pt <- post(r$fit, c("omega", "sigma", "tau[1]", "rho", "spatial_share"))
  kind <- sub("_r[0-9]+$", "", tag)
  if (!(kind %in% exemplar_done)) {
    plot_diagnostics(r$fit, pars = c("omega", "sigma", "rho", "tau[1]"), out_dir = FIG, prefix = paste0("a3_", kind), pairs_pars = c("omega", "sigma", "rho", "tau[1]"))
    exemplar_done <- c(exemplar_done, kind)
  }
  rows[[tag]] <- list(tag = tag, kind = tr$kind, placement = tr$placement, omega_true = if (is.null(tr$omega)) NA_real_ else tr$omega,
                      leak = if (is.null(tr$leak)) NA_real_ else tr$leak, rho_true = if (is.null(tr$rho)) NA_real_ else tr$rho,
                      fst = tr$realized_fst_between_regions, post = pt, log_lik = ll(r$fit), gate = r$gate$pass, reasons = paste(r$gate$reasons, collapse = ","),
                      rung = r$rung, seconds = r$seconds)
}
cell <- function(p, v) { i <- which(p$variable == v); sprintf("%.3f [%.3f, %.3f]", p$median[i], p$q5[i], p$q95[i]) }
L <- c("# A3 identification test with the placement supplied (step 6b; design section 4)", "",
       sprintf("Produced by `analyses/probe_latent_groups/validate_a3.R` on %s. Records from `a3_contracts.py`. Exploratory design validation, not evidence about the record.", format(Sys.Date())),
       "", "## Code equivalence", "",
       sprintf("A0 against A3 with omega fixed at 0 on `a3_own_G1_om0_r1`: largest |z| over %d shared posterior means %.2f at `%s` (threshold 3).", length(zsc), max(abs(zsc)), p0$variable[which.max(abs(zsc))]),
       diag_line(r0$diag, "A0"), diag_line(r3f$diag, "A3 fixed 0"), "",
       "## omega against the truth, placement supplied", "",
       "Prior half-Normal(0, 2): median 1.35, 90 percent below 3.3. A posterior that stays near that is the prior, not a detection.", "",
       "| record | kind | placement | truth omega / leak / rho | between-region F_ST | omega post | sigma | rho | tau | share | gate | rung | s |",
       "|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (x in rows) {
  tr_s <- if (!is.na(x$omega_true)) sprintf("omega %.2f, rho %.0f", x$omega_true, x$rho_true) else sprintf("leak %.2f", x$leak)
  L <- c(L, sprintf("| %s | %s | %s | %s | %.4f | %s | %s | %s | %s | %s | %s | %d | %.0f |", x$tag, x$kind, x$placement, tr_s, x$fst,
                    cell(x$post, "omega"), cell(x$post, "sigma"), cell(x$post, "rho"), cell(x$post, "tau[1]"), cell(x$post, "spatial_share"),
                    if (x$gate) "pass" else paste("fail:", x$reasons), x$rung, x$seconds))
}
writeLines(L, FIND)
write_json(list(equivalence = data.frame(variable = p0$variable, z = zsc), rows = rows, timing = timing),
           file.path(OUTD, "a3_validation_summary.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat("\nwrote", FIND, "\n")
