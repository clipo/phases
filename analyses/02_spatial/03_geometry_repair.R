#!/usr/bin/env Rscript
#' Does marginalising the nugget remove the ridge?
#'
#' The two-additive-fields form (composition_gp.stan) saturated max_treedepth on
#' 100% of iterations with ZERO divergences, on simulated and real data, under
#' both the uniform and the resolvable-range prior on rho. Zero divergences with
#' total treedepth saturation is a long flat direction, not a step-size problem,
#' and raising the prior's information content did not remove it.
#'
#' Diagnosis: u (iid, N x C) and f (GP, N x C) are additive at every site, so
#' 522 latent parameters compete to explain 261 effective observations. Small
#' rho makes f iid and therefore identical to u; large rho makes f constant and
#' therefore absorbed by alpha.
#'
#' Repair: composition_gp_marginal.stan puts the nugget in the covariance,
#' Cholesky(K + tau^2 I), leaving one latent field.
#'
#' This script measures the repair rather than assuming it, fitting BOTH models
#' to the SAME simulated data at known rho and comparing geometry and recovery.
#'
#' Usage: Rscript analyses/02_spatial/03_geometry_repair.R

suppressPackageStartupMessages({library(rstan); library(here); library(jsonlite)})
source(here("R", "fit_spatial.R"))
options(mc.cores = as.integer(Sys.getenv("MLS_MC_CORES", "2")))

IG_A <- 2.5874; IG_B <- 38.7141
rhos <- c(10, 20, 40)

cat(strrep("=", 80), "\n"); cat("GP geometry repair: two additive fields vs marginalised nugget\n")
cat(strrep("=", 80), "\n")

models <- list(
  additive   = stan_model(here("stan", "spatial", "composition_gp.stan")),
  marginal   = stan_model(here("stan", "spatial", "composition_gp_marginal.stan"))
)

rows <- list()
for (r in rhos) {
  d <- read_stan_data(here("data", "stan", sprintf("sim_rho%d.json", r)))
  d$prior_rho_type <- 1L; d$prior_rho_a <- IG_A; d$prior_rho_b <- IG_B
  d$rho_lo <- 0.5; d$rho_hi <- 300
  for (m in names(models)) {
    cat(sprintf("\ntrue rho = %d km, model = %s\n", r, m))
    res <- fit_composition_gp_ladder(models[[m]], d, use_tau = 1,
                                     cache_key = sprintf("repair_%s_rho%d", m, r))
    rows[[length(rows) + 1]] <- list(rho_true = r, model = m, res = res)
    cat(sprintf("  rho %.1f [%.1f, %.1f] width %.1fx | share %.2f [%.2f, %.2f] | treedepth %.0f%%\n",
                res$rho[1], res$rho[2], res$rho[3], res$rho[3]/res$rho[2],
                res$spatial_share[1], res$spatial_share[2], res$spatial_share[3],
                100*res$diag$n_treedepth/res$diag$n_total))
    cat("  ", diag_line(res$diag, sprintf("rung %d%s", res$rung,
        ifelse(res$passed, "", "  NOT PASSING"))), "\n", sep = "")
  }
}

g <- function(mm, f) vapply(rows[vapply(rows, function(x) x$model == mm, TRUE)], f, 0)
td_a <- g("additive", function(x) 100*x$res$diag$n_treedepth/x$res$diag$n_total)
td_m <- g("marginal", function(x) 100*x$res$diag$n_treedepth/x$res$diag$n_total)
w_a  <- g("additive", function(x) x$res$rho[3]/x$res$rho[2])
w_m  <- g("marginal", function(x) x$res$rho[3]/x$res$rho[2])
e_a  <- g("additive", function(x) x$res$diag$ess_bulk)
e_m  <- g("marginal", function(x) x$res$diag$ess_bulk)

L <- c("# Does marginalising the nugget remove the ridge?", "",
       sprintf("Produced by `analyses/02_spatial/03_geometry_repair.R` on %s. Both models fitted to the SAME simulated data at known rho, under the resolvable-range inverse-gamma prior, so the only difference is the latent structure.", format(Sys.Date())),
       "", "## Side by side", "",
       "| true rho | model | rho median | 95% CI | width | spatial share | treedepth sat. | bulk ESS |",
       "|---|---|---|---|---|---|---|---|")
for (x in rows) {
  rr <- x$res
  L <- c(L, sprintf("| %d | %s | %.1f | %.1f to %.1f | %.1fx | %.2f | **%.0f%%** | %.0f |",
                    x$rho_true, x$model, rr$rho[1], rr$rho[2], rr$rho[3],
                    rr$rho[3]/rr$rho[2], rr$spatial_share[1],
                    100*rr$diag$n_treedepth/rr$diag$n_total, rr$diag$ess_bulk))
}
L <- c(L, "", "## Diagnostics (rule 16, every fit)", "")
for (x in rows) L <- c(L, diag_line(x$res$diag,
  sprintf("rho %d, %s (rung %d%s)", x$rho_true, x$model, x$res$rung,
          ifelse(x$res$passed, "", ", NOT PASSING"))))
L <- c(L, "", "## Movement", "",
       sprintf("- Treedepth saturation: **%.0f%% additive -> %.0f%% marginalised** (mean).", mean(td_a), mean(td_m)),
       sprintf("- rho interval width: %.1fx -> %.1fx (mean).", mean(w_a), mean(w_m)),
       sprintf("- bulk ESS: %.0f -> %.0f (mean).", mean(e_a), mean(e_m)),
       "",
       if (mean(td_m) < 5 && mean(td_a) > 20)
         "**The repair works.** Removing the additive redundancy removes the flat direction. The two-field form should not be used again, and the geometry defect was structural rather than a matter of prior information or sampler settings."
       else if (mean(td_m) < mean(td_a) / 2)
         "**Partial.** The marginalised form is materially better but still saturating. Something beyond the u/f redundancy is contributing; the sigma-rho trade-off intrinsic to a GP on 29 points is the next candidate."
       else
         "**The repair does not work.** The flat direction is not the u/f redundancy. Do not report a length scale from either form until the cause is found; on current evidence the honest statement is that rho is not identified at this design.",
       "")
out <- here("output", "findings", "gp_geometry_repair.md")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
writeLines(L, out); cat("\nwrote", out, "\n")
