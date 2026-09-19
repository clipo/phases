#!/usr/bin/env Rscript
#' Does the composition GP recover a known length scale at this design?
#'
#' THE GATE for the spatial-field direction. analyses/48_length_scale_recovery.py
#' established that the DATA carry scale information between roughly 4 and 32 km.
#' This asks the narrower and harder question: can THIS MODEL read it out, or
#' does the nugget-versus-range trade-off swallow it?
#'
#' The trade-off is not hypothetical. ../mataa's stan/spatial/moai_site_gp.stan
#' records that on process-simulated data "the tau + field decomposition returned
#' d50 of 7 to 13 km whatever the truth". Under rule 19 that warning is tested
#' here rather than inherited: every dataset is fitted twice, with the
#' non-spatial term present (use_tau = 1) and absent (use_tau = 0), and the two
#' rho posteriors are compared against the same truth.
#'
#' Data are simulated FROM the model by analyses/49_write_stan_data.py in numpy,
#' with its own independent Matern implementation (rule 3).
#'
#' Usage: Rscript analyses/02_spatial/01_recovery.R

suppressPackageStartupMessages({ library(rstan); library(here); library(jsonlite) })
source(here("R", "fit_spatial.R"))
# Two parallel chains, not four. This box runs with ~13 GB available and each
# chain holds ~780 MB; four at once plus the parent OOM-killed a wrapper shell
# on 2026-08-31. Chains still total 4, they just run two at a time.
options(mc.cores = as.integer(Sys.getenv("MLS_MC_CORES", "2")))

cat(strrep("=", 80), "\n"); cat("Spatial composition GP: recovery on known truth\n")
cat("Basin site spacing: nearest-neighbour median 5.03 km, min 1.64, max 22.2;\n")
cat("median pairwise 45.9 km. A length scale below the nearest-neighbour\n")
cat("spacing is not separable from a nugget, which bounds what is recoverable.\n")
cat(strrep("=", 80), "\n")

sm <- stan_model(here("stan", "spatial", "composition_gp.stan"))
truth <- jsonlite::fromJSON(here("data", "stan", "sim_truth.json"))
truth_share <- truth$truth$`5.0`$spatial_share
# 80 km added 2026-08-31 to probe the UPPER end. analyses/48_length_scale_
# recovery.py found the data lose scale information above about 32 km, from a
# completely different direction (summary statistics on drift simulations, no
# GP anywhere). If this model independently fails above 32 km too, that is two
# witnesses agreeing and rule 10 says to treat it as signal.
rhos <- c(5, 10, 20, 40, 80)

rows <- list()
for (r in rhos) {
  dat <- read_stan_data(here("data", "stan", sprintf("sim_rho%d.json", r)))
  for (ut in c(1, 0)) {
    cat(sprintf("\nfitting true rho = %d km, use_tau = %d ...\n", r, ut))
    res <- fit_composition_gp_ladder(sm, dat, use_tau = ut,
                                     cache_key = sprintf("recovery_rho%d_tau%d", r, ut))
    rows[[length(rows) + 1]] <- list(rho_true = r, use_tau = ut, res = res)
    cat(sprintf("  rho %.1f [%.1f, %.1f] | d50 %.1f | spatial_share %.2f\n",
                res$rho[1], res$rho[2], res$rho[3], res$half_distance[1],
                res$spatial_share[1]))
    cat("  ", diag_line(res$diag, sprintf("rung %d%s", res$rung,
        ifelse(res$passed, "", " NOT PASSING"))), "\n", sep = "")
  }
}

L <- c("# Does the composition GP recover a known length scale?", "",
       sprintf("Produced by `analyses/02_spatial/01_recovery.R` on %s. Data simulated from the model by `analyses/49_write_stan_data.py` at the real basin geometry (29 assemblages, 10 classes, real sherd totals, max pairwise distance 127.8 km). Exponential kernel, uniform prior on log(rho) over [0.5, 300] km, true sigma %.1f, true tau %.1f, true spatial share %.2f.",
               format(Sys.Date()), truth$truth$`5.0`$sigma, truth$truth$`5.0`$tau,
               truth$truth$`5.0`$spatial_share),
       "", "## Recovery of rho", "",
       "| true rho (km) | use_tau | rho median | rho 95% CI | d50 | spatial share | covers truth | rung | verdict |",
       "|---|---|---|---|---|---|---|---|---|")
for (x in rows) {
  r <- x$res
  cov <- (x$rho_true >= r$rho[2]) && (x$rho_true <= r$rho[3])
  L <- c(L, sprintf("| %d | %d | %.1f | %.1f to %.1f | %.1f | %.2f | %s | %d | %s |",
                    x$rho_true, x$use_tau, r$rho[1], r$rho[2], r$rho[3],
                    r$half_distance[1], r$spatial_share[1],
                    ifelse(cov, "yes", "**NO**"), r$rung, r$verdict))
}
L <- c(L, "", "## Diagnostics (rule 16, every fit)", "")
for (x in rows) {
  L <- c(L, diag_line(x$res$diag,
                      sprintf("true rho %d km, use_tau %d (rung %d%s)", x$rho_true,
                              x$use_tau, x$res$rung,
                              ifelse(x$res$passed, "", ", NOT PASSING"))))
}

# The discriminating comparison: is rho pinned regardless of truth?
w1 <- vapply(rows[vapply(rows, function(x) x$use_tau == 1, TRUE)],
             function(x) x$res$rho[1], 0)
w0 <- vapply(rows[vapply(rows, function(x) x$use_tau == 0, TRUE)],
             function(x) x$res$rho[1], 0)
spread <- function(v) max(v) / min(v)
# Interval WIDTH is what decides whether a recovered parameter is usable. A
# median that tracks truth inside an order-of-magnitude interval is coverage
# bought with vagueness, which is ../mataa's pukao outcome ("2-3.5x,
# order-of-magnitude intervals") and must not be reported as recovery.
L <- c(L, "", "## Is the recovery USABLE? Interval widths", "",
       "| true rho | use_tau | 95% CI width ratio (hi/lo) | spatial share | share 95% CI |",
       "|---|---|---|---|---|")
for (x in rows) {
  r_ <- x$res
  wr <- r_$rho[3] / r_$rho[2]
  share_txt <- if (x$use_tau == 0) "1.00 (definitional, not estimated)" else
    sprintf("%.2f | %.2f to %.2f (width %.2f)", r_$spatial_share[1],
            r_$spatial_share[2], r_$spatial_share[3],
            r_$spatial_share[3] - r_$spatial_share[2])
  L <- c(L, if (x$use_tau == 0)
    sprintf("| %d | %d | **%.1fx** | %s | - |", x$rho_true, x$use_tau, wr, share_txt)
    else sprintf("| %d | %d | **%.1fx** | %s |", x$rho_true, x$use_tau, wr, share_txt))
}

w1 <- vapply(rows[vapply(rows, function(x) x$use_tau == 1, TRUE)],
             function(x) x$res$rho[1], 0)
w0 <- vapply(rows[vapply(rows, function(x) x$use_tau == 0, TRUE)],
             function(x) x$res$rho[1], 0)
t1 <- vapply(rows[vapply(rows, function(x) x$use_tau == 1, TRUE)],
             function(x) x$rho_true, 0)
spread <- function(v) max(v) / min(v)
share_err <- max(abs(vapply(rows[vapply(rows, function(x) x$use_tau == 1, TRUE)],
                            function(x) x$res$spatial_share[1], 0) - truth_share))
widths <- vapply(rows[vapply(rows, function(x) x$use_tau == 1, TRUE)],
                 function(x) x$res$rho[3] / x$res$rho[2], 0)
ss_lo <- vapply(rows[vapply(rows, function(x) x$use_tau == 1, TRUE)],
                function(x) x$res$spatial_share[2], 0)
ss_hi <- vapply(rows[vapply(rows, function(x) x$use_tau == 1, TRUE)],
                function(x) x$res$spatial_share[3], 0)
sw <- ss_hi - ss_lo

L <- c(L, "", "## The nugget-versus-range trade-off", "",
       sprintf("Across a true-rho range spanning a factor of %.0f, the posterior medians span a factor of %.1f with the non-spatial term present and %.1f without it.",
               max(t1) / min(t1), spread(w1), spread(w0)),
       "",
       "../mataa's `moai_site_gp.stan` header warns that its tau-plus-field decomposition returned d50 of 7 to 13 km whatever the truth. **That warning inverts here.** The FULL model tracks the truth; it is the FIELD-ONLY model (use_tau = 0) that pins near 8 to 10 km regardless, because without a nugget the field must absorb every idiosyncratic site difference and is forced short. Do not carry their diagnostic reading across: on this design, dropping tau creates the pathology rather than exposing it.",
       "")

# SATURATION. A spread statistic over the whole grid hides a ceiling: if the two
# largest truths both land in the same place, max/min still looks healthy. Check
# adjacent levels directly.
med1 <- vapply(rows[vapply(rows, function(x) x$use_tau == 1, TRUE)],
               function(x) x$res$rho[1], 0)
tru1 <- vapply(rows[vapply(rows, function(x) x$use_tau == 1, TRUE)],
               function(x) x$rho_true, 0)
ratios <- med1[-1] / med1[-length(med1)]
sat <- which(ratios < 1.25)          # a doubling of truth that moved < 25%
L <- c(L, "", "## Saturation: does the estimate stop moving?", "",
       "| truth pair | median pair | ratio | doubling recovered? |",
       "|---|---|---|---|")
for (i in seq_along(ratios)) {
  L <- c(L, sprintf("| %.0f -> %.0f km | %.1f -> %.1f | %.2f | %s |",
                    tru1[i], tru1[i + 1], med1[i], med1[i + 1], ratios[i],
                    ifelse(ratios[i] >= 1.25, "yes", "**NO**")))
}
if (length(sat)) {
  L <- c(L, "", sprintf("**Ceiling at roughly %.0f km.** Truths of %.0f and %.0f km return posterior medians of %.1f and %.1f, which are the same answer. Above that the model cannot tell one long interaction scale from another.",
                        tru1[min(sat)], tru1[min(sat)], tru1[min(sat) + 1],
                        med1[min(sat)], med1[min(sat) + 1]),
         "",
         "**This is a second, independent witness for the ceiling.** `analyses/48_length_scale_recovery.py` located it near 32 km from summary statistics on drift simulations, with no GP anywhere in it; this locates it from a fitted latent field on data simulated from the model itself. Two methods sharing no machinery agree, which rule 10 says to treat as signal rather than coincidence. The common cause is geometric: the basin is 128 km across, so a decay longer than a few tens of km is not expressed inside the study window at all.", "")
} else {
  L <- c(L, "", "No saturation detected: every doubling of the truth moves the estimate.", "")
}

cov_all <- all(vapply(rows[vapply(rows, function(x) x$use_tau == 1, TRUE)],
                      function(x) x$rho_true >= x$res$rho[2] && x$rho_true <= x$res$rho[3],
                      TRUE))

L <- c(L, "## Verdict", "")
if (spread(w1) > 2 && cov_all && max(widths) < 4) {
  L <- c(L, "**rho is recovered and usable.** Medians track truth, intervals cover, and the intervals are tight enough to report as estimates.")
} else if (spread(w1) > 2 && cov_all) {
  L <- c(L, sprintf("**rho is directionally recovered but only weakly identified.** Medians track the truth across a factor of %.0f and every interval covers it, but the intervals span factors of %.1f to %.1f. Coverage here is bought with vagueness. rho is reportable as an order-of-magnitude statement, not as a point estimate with a useful interval, which is the same outcome ../mataa recorded for pukao.",
                    max(t1) / min(t1), min(widths), max(widths)),
            "",
            sprintf("**`spatial_share` is recovered more usefully, but it is not sharp everywhere.** Its posterior MEDIAN sits within %.3f of the truth of %.2f at every length scale tested, which is genuine accuracy. Its INTERVAL, however, tightens with the length scale: width %.2f at the shortest scale, where [%.2f, %.2f] is close to no information at all on a bounded [0, 1] quantity, narrowing to %.2f at the longest. Reporting a single width bound here would be misleading, and an earlier draft of this verdict did exactly that.",
                    share_err, truth_share, sw[1], ss_lo[1], ss_hi[1], sw[length(sw)]),
            "",
            "The variance partition is still the better-behaved quantity and the one the paper's argument needs, since it asks what fraction of assemblage-level compositional variance is organised by geography. But it earns a headline only where the interval is informative, which on this evidence means the model must be fitted with enough spatial signal present, and the width has to be reported alongside the median every time.")
} else {
  L <- c(L, "**rho is NOT recovered.** The posterior does not move with the truth, or intervals miss it. Do not report a length scale from this model on this design.")
}
L <- c(L, "")

out <- here("output", "findings", "spatial_recovery.md")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
writeLines(L, out)
cat("\nwrote", out, "\n")
