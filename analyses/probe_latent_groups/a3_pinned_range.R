#!/usr/bin/env Rscript
# A3 with the field's range pinned (the conditional form of the claim; model-forms
# note, "range profiled"; step 6b follow-up, 2026-10-18). Exploratory design
# validation (rule 1). Usage:
#   Rscript analyses/probe_latent_groups/a3_pinned_range.R <rho_km>
# Fits A3 (placement supplied) to every a3 record with rho pinned at <rho_km>
# (uniform prior on [rho - 0.01, rho + 0.01]; A0's exponential kernel), and reads
# omega against the truth. With the range pinned short the field is close to
# site-level noise, so omega becomes the only shared-structure term and is
# identified; the price is that a smooth gradient at a longer range (the
# step-against-slope records) is read as region shifts, which the table shows.
suppressPackageStartupMessages({ library(here); library(rstan); library(posterior); library(jsonlite) })
here::i_am("pyproject.toml"); source(here("R", "fit_spatial.R")); options(mc.cores = 4L)
args <- commandArgs(trailingOnly = TRUE); RHO <- as.numeric(args[1]); stopifnot(is.finite(RHO), RHO > 0)
A3 <- here("stan", "spatial", "composition_gp_barrier.stan"); sm3 <- compile_model(A3)
OUTD <- here("output", "probe_latent_groups"); FIND <- here("output", "findings", sprintf("a3_pinned_rho%g.md", RHO))
GATED <- c("sigma", "tau", "alpha", "lp__", "f", "log_lik", "omega"); EXEMPT <- c("resp", "m_draw", "lp_m", "f_draw", "y_rep", "omega_free", "rho")
read_a3 <- function(f) { d <- read_stan_data(f); d$region <- matrix(as.integer(unlist(d$region)), nrow = d$P, byrow = TRUE); d }
post <- function(fit, pars) { dd <- posterior::subset_draws(posterior::as_draws_array(fit), variable = pars)
  as.data.frame(posterior::summarise_draws(dd, "median", ~posterior::quantile2(.x, probs = c(0.05, 0.95)))) }
files <- sort(list.files(here("data", "stan"), pattern = "^a3_(own|slope|sim)_.*\\.json$", full.names = TRUE))
rows <- list()
for (f in files) {
  tag <- sub("^a3_(.*)\\.json$", "\\1", basename(f)); d <- read_a3(f); tr <- attr(d, "provenance")$`_truth`
  d$prior_rho_type <- 3L; d$rho_lo <- RHO - 0.01; d$rho_hi <- RHO + 0.01
  r <- fit_cached(sprintf("a3pin%g_%s", RHO, tag), A3, d, gated = GATED, exempt = EXEMPT, chains = 4L, seed = 20261018L, sm = sm3)
  cat(diag_line(r$diag, sprintf("rho %g, %s", RHO, tag)), "\n")
  pt <- post(r$fit, c("omega", "sigma", "tau[1]", "spatial_share"))
  rows[[tag]] <- list(tag = tag, kind = tr$kind, placement = tr$placement, omega_true = if (is.null(tr$omega)) NA_real_ else tr$omega,
                      leak = if (is.null(tr$leak)) NA_real_ else tr$leak, rho_true = if (is.null(tr$rho)) NA_real_ else tr$rho,
                      fst = tr$realized_fst_between_regions, post = pt, gate = r$gate$pass, reasons = paste(r$gate$reasons, collapse = ","), rung = r$rung, seconds = r$seconds)
}
cell <- function(p, v) { i <- which(p$variable == v); sprintf("%.3f [%.3f, %.3f]", p$median[i], p$q5[i], p$q95[i]) }
L <- c(sprintf("# A3 with the range pinned at %g km (conditional form; step 6b follow-up)", RHO), "",
       sprintf("Produced by `analyses/probe_latent_groups/a3_pinned_range.R %g` on %s. A0's exponential kernel with rho on [%g, %g]; placement supplied. Exploratory design validation, not evidence about the record.", RHO, format(Sys.Date()), RHO - 0.01, RHO + 0.01), "",
       "| record | kind | placement | truth | between-region F_ST | omega post | sigma | tau | share | gate | rung | s |", "|---|---|---|---|---|---|---|---|---|---|---|---|")
for (x in rows) {
  tr_s <- if (!is.na(x$omega_true)) sprintf("omega %.2f, rho %.0f", x$omega_true, x$rho_true) else sprintf("leak %.2f", x$leak)
  L <- c(L, sprintf("| %s | %s | %s | %s | %.4f | %s | %s | %s | %s | %s | %d | %.0f |", x$tag, x$kind, x$placement, tr_s, x$fst, cell(x$post, "omega"), cell(x$post, "sigma"),
                    cell(x$post, "tau[1]"), cell(x$post, "spatial_share"), if (x$gate) "pass" else paste("fail:", x$reasons), x$rung, x$seconds))
}
writeLines(L, FIND); write_json(rows, file.path(OUTD, sprintf("a3_pinned_rho%g.json", RHO)), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat("\nwrote", FIND, "\n")
