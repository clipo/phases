#!/usr/bin/env Rscript
#' Does a resolvable-range prior on rho remove the ridge?
#'
#' THE PROBLEM. analyses/02_spatial/01_recovery.R found the full model hitting
#' 100% treedepth saturation with ZERO divergences at every length scale: a long
#' flat ridge, not a step-size failure. Measured cause: the uniform-on-log-rho
#' prior over [0.5, 300] km puts 36% of its mass BELOW the basin's median
#' nearest-neighbour spacing (5.03 km) and 13% ABOVE the maximum pairwise
#' distance (127.8 km). Half the prior sits where the likelihood is flat,
#' because no pair of sites can inform a correlation scale shorter than the
#' closest pair or longer than the farthest. The ridge is manufactured by the
#' prior, not discovered in the data.
#'
#' THE REPAIR. An inverse-gamma with vanishing density outside the resolvable
#' band, the standard construction for GP length scales. Solved for 1% tail mass
#' at each end: inv_gamma(2.5874, 38.7141), median 17.1 km, 90% interval
#' [6.8, 63.1] km.
#'
#' RULE 20, ALL THREE PARTS.
#'  Direction. This prior leans toward MODERATE length scales, with its median
#'    at 17 km against the manuscript's reported ~24 km. It is therefore
#'    SYMPATHETIC to the paper's operating point, which makes the burden heavy
#'    rather than light. Said plainly because rule 20 says a prior whose
#'    direction has not been stated has not been defended.
#'  (b) Sensitivity. Every dataset is fitted under BOTH priors and the movement
#'    is reported whether or not it is small.
#'  (c) Recovery against the prior. The grid includes rho = 5 km and rho = 80 km,
#'    which sit in the 1% tails the new prior disfavours. If those are still
#'    recovered, the prior is not doing the work. This is the decisive check.
#'
#' Usage: Rscript analyses/02_spatial/02_prior_repair.R

suppressPackageStartupMessages({ library(rstan); library(here); library(jsonlite) })
source(here("R", "fit_spatial.R"))
options(mc.cores = as.integer(Sys.getenv("MLS_MC_CORES", "2")))

IG_A <- 2.5874; IG_B <- 38.7141
NN_MEDIAN <- 5.03; MAX_D <- 127.8
rhos <- c(5, 10, 20, 40, 80)

cat(strrep("=", 80), "\n"); cat("Prior repair: uniform-on-log vs resolvable-range inverse-gamma\n")
cat(strrep("=", 80), "\n")

sm <- stan_model(here("stan", "spatial", "composition_gp.stan"))

set_prior <- function(dat, which) {
  if (which == "uniform") {
    dat$prior_rho_type <- 2L; dat$prior_rho_a <- 0; dat$prior_rho_b <- 0
    dat$rho_lo <- 0.5; dat$rho_hi <- 300
  } else {
    dat$prior_rho_type <- 1L; dat$prior_rho_a <- IG_A; dat$prior_rho_b <- IG_B
    dat$rho_lo <- 0.5; dat$rho_hi <- 300
  }
  dat
}

rows <- list()
for (r in rhos) {
  dat0 <- read_stan_data(here("data", "stan", sprintf("sim_rho%d.json", r)))
  for (pr in c("uniform", "invgamma")) {
    cat(sprintf("\ntrue rho = %d km, prior = %s\n", r, pr))
    res <- fit_composition_gp_ladder(sm, set_prior(dat0, pr), use_tau = 1,
                                     cache_key = sprintf("prior_%s_rho%d", pr, r))
    rows[[length(rows) + 1]] <- list(rho_true = r, prior = pr, res = res)
    cat(sprintf("  rho %.1f [%.1f, %.1f] width %.1fx | share %.2f [%.2f, %.2f] | treedepth %.0f%%\n",
                res$rho[1], res$rho[2], res$rho[3], res$rho[3] / res$rho[2],
                res$spatial_share[1], res$spatial_share[2], res$spatial_share[3],
                100 * res$diag$n_treedepth / res$diag$n_total))
  }
}

get <- function(pr, f) vapply(rows[vapply(rows, function(x) x$prior == pr, TRUE)], f, 0)
wid_u <- get("uniform",  function(x) x$res$rho[3] / x$res$rho[2])
wid_i <- get("invgamma", function(x) x$res$rho[3] / x$res$rho[2])
td_u  <- get("uniform",  function(x) 100 * x$res$diag$n_treedepth / x$res$diag$n_total)
td_i  <- get("invgamma", function(x) 100 * x$res$diag$n_treedepth / x$res$diag$n_total)
sw_u  <- get("uniform",  function(x) x$res$spatial_share[3] - x$res$spatial_share[2])
sw_i  <- get("invgamma", function(x) x$res$spatial_share[3] - x$res$spatial_share[2])
med_i <- get("invgamma", function(x) x$res$rho[1])
cov_i <- vapply(rows[vapply(rows, function(x) x$prior == "invgamma", TRUE)],
                function(x) x$rho_true >= x$res$rho[2] && x$rho_true <= x$res$rho[3], TRUE)

L <- c("# Does a resolvable-range prior remove the ridge?", "",
       sprintf("Produced by `analyses/02_spatial/02_prior_repair.R` on %s.", format(Sys.Date())),
       "", "## Why the uniform prior was the problem", "",
       sprintf("The basin's median nearest-neighbour spacing is %.2f km and its maximum pairwise distance is %.1f km. No pair of sites can inform a correlation length shorter than the closest pair or longer than the farthest, so the likelihood is flat outside that band. Uniform on log(rho) over [0.5, 300] km places **%.0f%% of its mass below %.2f km and %.0f%% above %.1f km**, so about half the prior sits where the data say nothing. That is what the sampler was exploring for 100%% of its trajectories.",
               NN_MEDIAN, MAX_D, 36, NN_MEDIAN, 13, MAX_D),
       "",
       sprintf("The replacement is inv_gamma(%.4f, %.4f), solved for 1%% tail mass at each end of that band: median %.1f km, 90%% interval [%.1f, %.1f] km.",
               IG_A, IG_B, 17.1, 6.8, 63.1),
       "", "## Side by side", "",
       "| true rho | prior | rho median | 95% CI | width | share width | treedepth sat. |",
       "|---|---|---|---|---|---|---|")
for (x in rows) {
  rr <- x$res
  L <- c(L, sprintf("| %d | %s | %.1f | %.1f to %.1f | %.1fx | %.2f | %.0f%% |",
                    x$rho_true, x$prior, rr$rho[1], rr$rho[2], rr$rho[3],
                    rr$rho[3] / rr$rho[2],
                    rr$spatial_share[3] - rr$spatial_share[2],
                    100 * rr$diag$n_treedepth / rr$diag$n_total))
}
L <- c(L, "", "## Diagnostics (rule 16, every fit)", "")
for (x in rows) L <- c(L, diag_line(x$res$diag, sprintf("rho %d, %s (rung %d%s)",
                                                        x$rho_true, x$prior, x$res$rung,
                                                        ifelse(x$res$passed, "", ", NOT PASSING"))))
L <- c(L, "", "## Movement", "",
       sprintf("- Treedepth saturation: %.0f%% mean under uniform, **%.0f%% under inverse-gamma**.", mean(td_u), mean(td_i)),
       sprintf("- rho interval width: %.1fx mean under uniform, **%.1fx under inverse-gamma**.", mean(wid_u), mean(wid_i)),
       sprintf("- spatial_share interval width: %.2f mean under uniform, **%.2f under inverse-gamma**.", mean(sw_u), mean(sw_i)),
       "",
       "## Rule 20(c): recovery at values the prior disfavours", "",
       sprintf("rho = 5 km and rho = 80 km sit in the prior's 1%% tails. Under the inverse-gamma they return %.1f and %.1f, covering truth: %s and %s. %s",
               med_i[1], med_i[length(med_i)],
               ifelse(cov_i[1], "yes", "**NO**"),
               ifelse(cov_i[length(cov_i)], "yes", "**NO**"),
               if (cov_i[1] && cov_i[length(cov_i)])
                 "The prior is not dragging the estimate to its own median; the data move it into both tails."
               else
                 "**At least one disfavoured truth is NOT recovered. The prior is doing work the data should be doing, and this prior is not reportable under rule 20(c).**"),
       "")
out <- here("output", "findings", "spatial_prior_repair.md")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
writeLines(L, out)
cat("\nwrote", out, "\n")
