#!/usr/bin/env Rscript
# A0 refit on the real St. Francis basin data with the corrected diagnostics
# layer, its visual diagnostic set, and the rule-20 program for A0's priors.
# Latent-group probe, pilot step 4 (spec revision 2.1, sections 4.1 and 6.4).
#
# Exploratory design validation (rule 1: committed even though exploratory).
# Every fit goes through R/fit_spatial.R, so every posterior below carries its
# diagnostics, its source hashes and the data fingerprint; fits are cached by
# content in output/gp_cache, so a rerun with unchanged inputs reads them back.
#
# Run from the repository root:
#   Rscript analyses/probe_latent_groups/refit_a0_basin.R
# Inputs: data/stan/basin_composition.json and data/stan/a0_recovery_*.json
#   (analyses/49_write_stan_data.py). Outputs: output/findings/a0_basin_refit.md,
#   output/probe_latent_groups/*.{json,csv}, figures/probe_latent_groups/a0_*.png.
# The prior- and posterior-predictive counts written here are summarized by
#   analyses/probe_latent_groups/a0_predictive_summaries.py (Python owns the
#   resemblance summaries, spec 3.3).

suppressPackageStartupMessages({ library(here); library(rstan); library(posterior); library(jsonlite) })
here::i_am("pyproject.toml")
source(here("R", "fit_spatial.R"))
source(here("R", "diagnostic_plots.R"))
options(mc.cores = 4L)

MODEL <- here("stan", "spatial", "composition_gp_marginal.stan")
BASE_COMMIT <- "bb2fd21"   # last commit before the 2026-10-16 additions to the model file
FIG  <- here("figures", "probe_latent_groups")
OUTD <- here("output", "probe_latent_groups")
FIND <- here("output", "findings", "a0_basin_refit.md")
for (d in c(FIG, OUTD, dirname(FIND))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

GATED <- c("sigma", "tau", "rho", "alpha", "lp__", "f", "log_lik")   # spec 4.2
SEED  <- 20261016L
N_PRED_DRAWS <- 200L

sm  <- compile_model(MODEL)
dat <- read_stan_data(here("data", "stan", "basin_composition.json"))
prov <- attr(dat, "provenance")
cat(sprintf("contract: N %d, K %d, sherds %d, max distance %.2f km, counts sha256 %s\n",
            dat$N, dat$K, sum(dat$y), prov$max_pairwise_km, substr(prov$sha256_counts, 1, 12)))

post_table <- function(fit, pars) {
  d <- posterior::subset_draws(posterior::as_draws_array(fit), variable = pars)
  s <- posterior::summarise_draws(d, "median", ~posterior::quantile2(.x, probs = c(0.05, 0.95)))
  s <- as.data.frame(s)
  stopifnot(all(c("variable", "median", "q5", "q95") %in% names(s)))
  s
}
fmt_row <- function(s, i) sprintf("| %s | %.3f | %.3f | %.3f |", s$variable[i], s$median[i], s$q5[i], s$q95[i])

write_yrep <- function(fit, path, n_draws = N_PRED_DRAWS, seed = 1L) {
  yr <- rstan::extract(fit, pars = "y_rep")$y_rep        # draws x N x K
  set.seed(seed); idx <- sort(sample(dim(yr)[1], min(n_draws, dim(yr)[1])))
  rows <- do.call(rbind, lapply(seq_along(idx), function(j) {
    m <- yr[idx[j], , ]
    data.frame(draw = j, i = rep(seq_len(nrow(m)), ncol(m)), k = rep(seq_len(ncol(m)), each = nrow(m)),
               count = as.integer(m))
  }))
  utils::write.csv(rows, path, row.names = FALSE)
  invisible(path)
}

timing <- list()
note <- function(label, res) timing[[label]] <<- list(rung = res$rung, seconds = res$seconds, cached = isTRUE(res$cached),
                                                     pass = res$gate$pass, reasons = res$gate$reasons)

# ---------------------------------------------------------------------------
# 1. The refit on the real data
# ---------------------------------------------------------------------------
cat("\n[1] A0 on the real basin data\n")
res <- fit_cached("a0_basin", MODEL, dat, gated = GATED, chains = 4L, seed = SEED, sm = sm)
note("a0_basin", res)
cat(diag_line(res$diag, "a0_basin"), "\n")
main_pars <- c("sigma", "rho", "tau[1]", paste0("alpha[", 1:9, "]"))
figs <- plot_diagnostics(res$fit, pars = main_pars, out_dir = FIG, prefix = "a0_basin",
                         pairs_pars = c("sigma", "rho", "tau[1]"))
# A second pairs plot on the intercepts alone (spec 4.2 lists alpha in the pairs set).
alpha_pairs <- file.path(FIG, "a0_basin_pairs_alpha.png")
save_png(bayesplot::mcmc_pairs(posterior::as_draws_array(res$fit), pars = paste0("alpha[", 1:9, "]"),
                               np = bayesplot::nuts_params(res$fit),
                               off_diag_args = list(size = 0.3, alpha = 0.3)), alpha_pairs, height = 7)
write_yrep(res$fit, file.path(OUTD, "a0_posterior_predictive_yrep.csv"))
real_post <- post_table(res$fit, c("sigma", "tau[1]", "rho", "spatial_share", "half_distance", "practical_range"))

# ---------------------------------------------------------------------------
# 2a. Prior predictive (likelihood switched off in the model's own data)
# ---------------------------------------------------------------------------
cat("\n[2a] prior predictive\n")
dat0 <- dat; dat0$likelihood <- 0L
res0 <- fit_cached("a0_prior_predictive", MODEL, dat0, gated = setdiff(GATED, "log_lik"), exempt = "log_lik",
                   chains = 4L, seed = SEED,
                   ladder = list(list(iter = 2000L, warmup = 1000L, adapt_delta = 0.9, max_treedepth = 10L)), sm = sm)
note("a0_prior_predictive", res0)
cat(diag_line(res0$diag, "a0_prior_predictive"), "\n")
write_yrep(res0$fit, file.path(OUTD, "a0_prior_predictive_yrep.csv"))

# ---------------------------------------------------------------------------
# 2b. Sensitivity refits, each leaning against a default (spec 4.1)
# ---------------------------------------------------------------------------
cat("\n[2b] sensitivity refits\n")
dat_s <- dat; dat_s$prior_scale_sd <- 0.5
res_s <- fit_cached("a0_sens_scale05", MODEL, dat_s, gated = GATED, chains = 4L, seed = SEED, sm = sm)
note("a0_sens_scale05", res_s); cat(diag_line(res_s$diag, "a0_sens_scale05"), "\n")
dat_r <- dat; dat_r$prior_rho_type <- 3L
res_r <- fit_cached("a0_sens_rho_uniform", MODEL, dat_r, gated = GATED, chains = 4L, seed = SEED, sm = sm)
note("a0_sens_rho_uniform", res_r); cat(diag_line(res_r$diag, "a0_sens_rho_uniform"), "\n")
sens <- lapply(list(default = res, scale05 = res_s, rho_uniform = res_r),
               function(r) post_table(r$fit, c("sigma", "tau[1]", "rho", "spatial_share", "half_distance")))

# Identified quantities (D20): the data see between-site contrasts, whose
# variance is 2 sigma^2 (1 - exp(-d/rho)); on 28 sites within 79 km only the
# contrast scale sigma^2/rho and the contrast SD at distances inside the
# basin are informed, while sigma, rho and alpha individually are not.
contrast_table <- function(fit, dists = c(4.3, 36.6, 79.1)) {
  d <- posterior::as_draws_matrix(fit); sg <- d[, "sigma"]; rho <- d[, "rho"]; tau <- d[, "tau[1]"]
  q <- function(x) sprintf("%.2f [%.2f, %.2f]", stats::median(x), stats::quantile(x, 0.05), stats::quantile(x, 0.95))
  rows <- c(`sigma^2/rho (logit^2 per km)` = q(sg^2 / rho))
  for (dd in dists) {
    rows[sprintf("field contrast SD at %.1f km (logits)", dd)] <- q(sqrt(2 * sg^2 * (1 - exp(-dd / rho))))
    rows[sprintf("contrast SD with nugget at %.1f km", dd)] <- q(sqrt(2 * sg^2 * (1 - exp(-dd / rho)) + 2 * tau^2))
  }
  rows
}
ident <- lapply(list(default = res, scale05 = res_s, rho_uniform = res_r), function(r) contrast_table(r$fit))

# ---------------------------------------------------------------------------
# 2c. Recovery against the prior (spec 4.1c): four own-generator records
# ---------------------------------------------------------------------------
cat("\n[2c] recovery fits\n")
rec_files <- sort(list.files(here("data", "stan"), pattern = "^a0_recovery_share.*\\.json$", full.names = TRUE))
recovery <- list()
for (f in rec_files) {
  tag <- sub("^a0_recovery_(.*)\\.json$", "\\1", basename(f))
  d <- read_stan_data(f); truth <- attr(d, "provenance")$`_truth`
  r <- fit_cached(paste0("a0_recovery_", tag), MODEL, d, gated = GATED, chains = 4L, seed = SEED, sm = sm)
  note(paste0("a0_recovery_", tag), r); cat(diag_line(r$diag, tag), "\n")
  pt <- post_table(r$fit, c("spatial_share", "rho", "sigma", "tau[1]"))
  tv <- c(spatial_share = truth$spatial_share, rho = truth$rho, sigma = truth$sigma, `tau[1]` = truth$tau)
  pt$truth <- tv[pt$variable]
  pt$covered <- pt$truth >= pt$q5 & pt$truth <= pt$q95
  recovery[[tag]] <- list(table = pt, diag = diag_line(r$diag, tag), gate = r$gate)
}

# ---------------------------------------------------------------------------
# Findings file
# ---------------------------------------------------------------------------
diff_txt <- tryCatch(system2("git", c("diff", BASE_COMMIT, "--", "stan/spatial/composition_gp_marginal.stan"),
                             stdout = TRUE, stderr = TRUE), error = function(e) "git diff unavailable")
ver <- software_versions()
L <- c(
  "# A0 refit on the real basin data, with corrected diagnostics and the rule-20 program",
  "",
  sprintf("Produced by `analyses/probe_latent_groups/refit_a0_basin.R` on %s (pilot step 4 of the latent-group probe; exploratory design validation, not evidence about the record).", format(Sys.Date())),
  "",
  "## Data contract and provenance",
  "",
  sprintf("`data/stan/basin_composition.json`: %d assemblages, %d classes, %d sherds. Counts sha256 `%s`; IDs `%s`; coordinates `%s`; distance matrix `%s`; raw workbook `%s`.",
          dat$N, dat$K, sum(dat$y), prov$sha256_counts, prov$sha256_ids, prov$sha256_coords, prov$sha256_distance_km,
          if (is.null(prov$sha256_raw_pfg_xlsx)) "absent" else prov$sha256_raw_pfg_xlsx),
  sprintf("Geometry (recomputed by the fingerprint): maximum pairwise %.2f km, median pairwise %.2f km, nearest-neighbor min %.2f, median %.2f, max %.2f km. Provisional coordinates: %s.",
          prov$max_pairwise_km, prov$median_pairwise_km, prov$nearest_neighbor_min_km, prov$nearest_neighbor_median_km,
          prov$nearest_neighbor_max_km, paste(unlist(prov$provisional_coordinates), collapse = ", ")),
  sprintf("Prior in the contract: `prior_rho_type = %d` (1 inverse-gamma, 2 proportional to 1/rho, 3 uniform), rho on [%.1f, %.0f] km, half-Normal(0, %.1f) on sigma and tau, kernel %d. This settles the lineage question of inventory row 14: the contract Carl's basin fit read carries the 1/rho prior; only the two repair scripts switched to the inverse-gamma.",
          dat$prior_rho_type, dat$rho_lo, dat$rho_hi, dat$prior_scale_sd, dat$kernel),
  sprintf("Software: %s; rstan %s; StanHeaders %s; posterior %s; loo %s; bayesplot %s. Stan source sha256 `%s`; fit layer sha256 `%s`.",
          ver$R, ver$rstan, ver$StanHeaders, ver$posterior, ver$loo, ver$bayesplot, stan_source_hash(MODEL)$sha256, layer_hash()),
  "",
  sprintf("## Changes to the model file since %s (model block unchanged except the third rho-prior case)", BASE_COMMIT),
  "", "```diff", diff_txt, "```", "",
  "## 1. The refit",
  "",
  sprintf("Gate: %s%s. Rung %d (iter %d, warmup %d, adapt_delta %.2f, max_treedepth %d); %.0f s wall time (4 chains in parallel); %s.",
          if (res$gate$pass) "PASS" else "FAIL", if (res$gate$pass) "" else paste0(" (", paste(res$gate$reasons, collapse = ", "), ")"),
          res$rung, res$rung_cfg$iter, res$rung_cfg$warmup, res$rung_cfg$adapt_delta, res$rung_cfg$max_treedepth,
          res$seconds, if (isTRUE(res$cached)) "read from cache" else "fitted in this run"),
  diag_line(res$diag, "a0_basin"),
  sprintf("Treedepth, measured against the configured maximum for the first time: %d of %d post-warmup iterations (%.1f%%) at depth %d; largest depth observed %d. Gated variables: %d (sigma, tau, rho, 9 alpha, lp__, 252 field coordinates, 28 pointwise log-likelihoods); exempt %d; stuck %d.",
          res$diag$treedepth$n_at_max, res$diag$treedepth$n_total, res$diag$treedepth$pct_at_max,
          res$diag$treedepth$max_treedepth, res$diag$treedepth$max_observed, res$diag$n_gated, res$diag$n_exempt, res$diag$n_stuck),
  "",
  "| quantity | median | 5% | 95% |", "|---|---|---|---|",
  vapply(seq_len(nrow(real_post)), function(i) fmt_row(real_post, i), character(1)),
  "",
  sprintf("Visual diagnostics (%s, 300 dpi, 7 in): %s; plus `%s` (pairs on the nine intercepts).",
          figs$font[1], paste(sprintf("`%s`", basename(figs$file)), collapse = ", "), basename(alpha_pairs)),
  "Posterior-predictive counts for 200 draws: `output/probe_latent_groups/a0_posterior_predictive_yrep.csv` (summarized by the Python script).",
  "",
  "## 2a. Prior predictive (rule 20a)",
  "",
  sprintf("Sampled with the model's own `likelihood = 0` switch. %s", diag_line(res0$diag, "a0_prior_predictive")),
  "Counts for 200 draws: `output/probe_latent_groups/a0_prior_predictive_yrep.csv`. The resemblance summaries (mean within-assemblage Gini-Simpson, mean richness, pooled Gini-Simpson, two-cluster plug-in F_ST) are computed beside the observed values by `analyses/probe_latent_groups/a0_predictive_summaries.py`, which owns those definitions (spec 3.3); its output is `output/findings/a0_predictive_summaries.md`.",
  "",
  "## 2b. Prior sensitivity (rule 20b)",
  "",
  "Direction of each default and what the alternative leans toward: half-Normal(0, 2) on sigma and tau admits large variance of both kinds; half-Normal(0, 0.5) leans toward little of either. The 1/rho prior puts half its mass below 12.2 km (the geometric mean of 0.5 and 300); the uniform prior on [0.5, 300] puts 90 percent of its mass above 30 km and leans toward long ranges.",
  "",
  "| quantity | default median [5%, 95%] | half-Normal(0, 0.5) on sigma, tau | uniform rho |", "|---|---|---|---|",
  vapply(seq_len(nrow(sens$default)), function(i) {
    cell <- function(s) sprintf("%.3f [%.3f, %.3f]", s$median[i], s$q5[i], s$q95[i])
    sprintf("| %s | %s | %s | %s |", sens$default$variable[i], cell(sens$default), cell(sens$scale05), cell(sens$rho_uniform))
  }, character(1)),
  "",
  diag_line(res_s$diag, "a0_sens_scale05"), diag_line(res_r$diag, "a0_sens_rho_uniform"),
  "",
  "### What is identified (D20)",
  "",
  "The data see between-site contrasts, with variance 2 sigma^2 (1 - exp(-d/rho)). On 28 sites within 79 km the contrast scale sigma^2/rho and the contrast SD at distances inside the basin are informed; sigma, rho and the intercepts individually are not (the intercepts are anti-correlated with the field's per-class mean, correlations -0.96 to -1.00 on the default fit). Distances below are the nearest-neighbor median, the median pair and the maximum pair. Length scales in km and sigma as an amplitude are not reportable as estimates; these are.",
  "",
  "| quantity | default prior | half-Normal(0, 0.5) on sigma, tau | uniform rho |", "|---|---|---|---|",
  vapply(names(ident$default), function(k) sprintf("| %s | %s | %s | %s |", k, ident$default[[k]], ident$scale05[[k]], ident$rho_uniform[[k]]), character(1)),
  "",
  "## 2c. Recovery against the prior (rule 20c)",
  "",
  "Data simulated from the model by `analyses/49_write_stan_data.simulate_from_model` (numpy, no Stan) at the real geometry and sherd totals; spatial share 0.2 and 0.9 crossed with rho 5 and 60 km; total latent variance 1 logit^2. Covered = truth inside the 5 to 95 percent interval.",
  "",
  "| record | quantity | truth | median | 5% | 95% | covered |", "|---|---|---|---|---|---|---|",
  unlist(lapply(names(recovery), function(tag) {
    pt <- recovery[[tag]]$table
    vapply(seq_len(nrow(pt)), function(i) sprintf("| %s | %s | %.3f | %.3f | %.3f | %.3f | %s |", tag, pt$variable[i], pt$truth[i],
                                                 pt$median[i], pt$q5[i], pt$q95[i], if (pt$covered[i]) "yes" else "no"), character(1))
  })),
  "",
  vapply(names(recovery), function(tag) recovery[[tag]]$diag, character(1)),
  "",
  "## Timing and gates",
  "",
  "| fit | rung | seconds | gate | reasons | cached |", "|---|---|---|---|---|---|",
  vapply(names(timing), function(nm) { t <- timing[[nm]]
    sprintf("| %s | %d | %.0f | %s | %s | %s |", nm, t$rung, t$seconds, if (t$pass) "pass" else "fail",
            if (length(t$reasons)) paste(t$reasons, collapse = ", ") else "", if (t$cached) "yes" else "no") }, character(1)),
  "",
  "Fit records (JSON, with every diagnostic, hash and version) are beside the cached fits in `output/gp_cache/`."
)
writeLines(L, FIND)
jsonlite::write_json(list(timing = timing, real = real_post, sensitivity = sens,
                          recovery = lapply(recovery, `[[`, "table"), identified = ident, figures = figs,
                          treedepth = res$diag$treedepth, gate = res$gate[c("pass", "reasons")]),
                     file.path(OUTD, "a0_refit_summary.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat("\nwrote", FIND, "\n")
