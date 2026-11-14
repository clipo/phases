# A1 under the corrected intercept prior (D25 follow-up; closing step 3, item 1 of D36).
#
# D24 closes A1 on recovery tables fitted under the inherited Normal(0, 2)
# intercept prior, which D25 replaced with Normal(0, 10). The A1 records in
# data/stan/ were regenerated under D25 (prior_alpha_sd = 10), so the old
# fits no longer match them and validate_a1.R would refit all thirteen. This
# script refits three: the zero-shift record (null calibration) and two
# planted records at group scale 0.5 (sg0.5_r1, which failed the gate under
# the old prior with divergences, and sg0.5_r5, which passed), and compares
# sigma_g, Delta, Delta_w and coherence with the Normal(0, 2) values stored in
# a1_validation_summary.json. Every fit carries the rule 16 diagnostics and
# the visual harness is written for the null fit. Exploratory design
# validation, not evidence about the record. Descriptive comparison only.
#
# Run: Rscript analyses/probe_latent_groups/a1_prior_check.R
# Output: output/findings/a1_prior_check.md, output/probe_latent_groups/a1_prior_check.json,
#         figures/probe_latent_groups/a1prior_null_{rank,trace,pairs,energy}.png
suppressPackageStartupMessages({ library(here); library(jsonlite) })
source(here("R", "fit_spatial.R")); source(here("R", "diagnostic_plots.R"))
options(mc.cores = 4L)
A1 <- here("stan", "spatial", "composition_gp_groups.stan")
OUTD <- here("output", "probe_latent_groups"); FIG <- here("figures", "probe_latent_groups")
FIND <- here("output", "findings", "a1_prior_check.md")
OLD <- fromJSON(file.path(OUTD, "a1_validation_summary.json"), simplifyVector = FALSE)
SEED <- 20261017L   # the seed validate_a1.R used, so the only change is the prior field in the records
GATED_A1 <- c("sigma", "tau", "rho", "alpha", "lp__", "f", "log_lik", "sigma_g", "Delta", "Delta_w", "coherence", "comembership")
EXEMPT_A1 <- c("delta", "pi", "resp", "z_g", "sigma_g_free", "z_draw", "coherence_hard", "y_rep")
PARS <- c("sigma_g", "Delta", "Delta_w", "coherence")
RECORDS <- c(zero = "a1_zero.json", sg0.5_r1 = "a1_sg0.5_r1.json", sg0.5_r5 = "a1_sg0.5_r5.json")
sm1 <- compile_model(A1)
read_a1 <- function(f) { d <- read_stan_data(f); d$edges <- matrix(as.integer(unlist(d$edges)), nrow = d$M, byrow = FALSE); d }
post <- function(fit, pars) {
  d <- posterior::subset_draws(posterior::as_draws_array(fit), variable = pars)
  as.data.frame(posterior::summarise_draws(d, "median", ~posterior::quantile2(.x, probs = c(0.05, 0.95))))
}
old_post <- function(tag) {
  p <- if (tag == "zero") OLD$null$post else OLD$recovery[[tag]]$post
  do.call(rbind, lapply(p, function(r) data.frame(variable = r$variable, median = r$median, q5 = r$q5, q95 = r$q95)))
}
cell <- function(p, v) { i <- which(p$variable == v); if (length(i) == 0) return("not stored"); sprintf("%.3f [%.3f, %.3f]", p$median[i], p$q5[i], p$q95[i]) }

rows <- list(); out <- list()
for (tag in names(RECORDS)) {
  d <- read_a1(here("data", "stan", RECORDS[[tag]]))
  if (!isTRUE(all.equal(d$prior_alpha_sd, 10))) stop(tag, ": the record does not carry prior_alpha_sd = 10")
  tr <- attr(d, "provenance")$`_truth`
  r <- fit_cached(paste0("a1prior_", tag), A1, d, gated = GATED_A1, exempt = EXEMPT_A1, constants = "comembership",
                  comembership = "comembership", chains = 4L, seed = SEED, sm = sm1)
  cat(diag_line(r$diag, tag), "\n")
  p_new <- post(r$fit, PARS); p_old <- old_post(tag)
  if (tag == "zero") plot_diagnostics(r$fit, pars = c("sigma", "rho", "tau[1]", "sigma_g", "Delta", "coherence", "alpha[1]", "alpha[2]"),
                                      out_dir = FIG, prefix = "a1prior_null", pairs_pars = c("sigma", "rho", "tau[1]", "sigma_g"))
  rows[[tag]] <- c(list(record = tag, truth = if (is.null(tr$sigma_g)) 0 else tr$sigma_g),
                   setNames(lapply(PARS, function(v) cell(p_old, v)), paste0("old_", PARS)),
                   setNames(lapply(PARS, function(v) cell(p_new, v)), paste0("new_", PARS)),
                   list(gate = if (r$gate$pass) "pass" else paste("fail:", paste(r$gate$reasons, collapse = ",")),
                        rung = r$rung, seconds = round(r$seconds), diag = diag_line(r$diag, tag)))
  out[[tag]] <- list(truth_sigma_g = rows[[tag]]$truth, post_new = p_new, post_old = p_old, gate = r$gate, rung = r$rung, seconds = r$seconds, diag = r$diag)
}
write_json(out, file.path(OUTD, "a1_prior_check.json"), digits = NA, auto_unbox = TRUE, pretty = TRUE)

L <- c("# A1 under the corrected intercept prior (D25 follow-up)", "",
       paste0("Produced by `analyses/probe_latent_groups/a1_prior_check.R` on ", format(Sys.Date()), ". Three A1 fits on the D25 records ",
              "(prior_alpha_sd = 10; seed ", SEED, ", the seed of validate_a1.R), against the Normal(0, 2) posteriors stored in ",
              "`a1_validation_summary.json`. Median [5%, 95%]. The gate and its reasons are validate_a1.R's. Exploratory design ",
              "validation, not evidence about the record; a descriptive comparison, no threshold."), "",
       "| record | true sigma_g | sigma_g, Normal(0, 2) | sigma_g, Normal(0, 10) | Delta, Normal(0, 2) | Delta, Normal(0, 10) | Delta_w, Normal(0, 2) | Delta_w, Normal(0, 10) | coherence, Normal(0, 2) | coherence, Normal(0, 10) | gate (new) | rung | seconds |",
       "|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (tag in names(rows)) {
  x <- rows[[tag]]
  L <- c(L, sprintf("| %s | %.2f | %s | %s | %s | %s | %s | %s | %s | %s | %s | %d | %d |", x$record, x$truth,
                    x$old_sigma_g, x$new_sigma_g, x$old_Delta, x$new_Delta, x$old_Delta_w, x$new_Delta_w, x$old_coherence, x$new_coherence,
                    x$gate, as.integer(x$rung), as.integer(x$seconds)))
}
L <- c(L, "", "## Diagnostics (rule 16)", "")
for (tag in names(rows)) L <- c(L, rows[[tag]]$diag)
L <- c(L, "", "Visual set for the null fit: `figures/probe_latent_groups/a1prior_null_{rank,trace,pairs,energy}.png`.", "")
writeLines(L, FIND); cat(paste(L, collapse = "\n"), "\n"); cat("wrote", FIND, "\n")
