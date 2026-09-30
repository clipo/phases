# 05_basin_fit.R -- the basin GP fit the paper's headline number comes from.
#
# WHY THIS EXISTS. The abstract reports that 96 percent of assemblage-level
# compositional variation is spatially structured. Until this script, no
# committed code produced that fit. `output/gp_real_*.rds` existed in the tree
# and `04_posterior_predictive.R` read it, but nothing in the repository wrote
# it: it came from an ad-hoc interactive call. That violates rule 1 (no number
# without a committed procedure) on the single most prominent quantity in the
# paper, and it meant the fit could not be re-derived, re-seeded or audited.
#
# It also violated rule 16 in effect. The diagnostics are attached to every fit
# by construction in R/fit_spatial.R, so they existed, but nothing reported them
# alongside the number, and the marginalised fit does NOT pass: it saturates
# treedepth because the length scale is not identified above this design's
# resolving ceiling. That is a real property of the data, argued in
# output/findings/gp_posterior_predictive.md and gp_geometry_repair.md, and it
# is the reason the paper reports a variance SHARE and refuses to report a
# length-scale NUMBER. Reporting the share while silently omitting the
# diagnostics is not acceptable; reporting both is.
#
# Writes output/gp_real_<model>.rds (consumed by 04_posterior_predictive.R) and
# output/findings/gp_basin_fit.md.
#
# Usage: Rscript analyses/02_spatial/05_basin_fit.R

suppressPackageStartupMessages({
  library(here)
  library(rstan)
})

source(here("R", "fit_spatial.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))

DATA <- here("data", "stan", "basin_composition.json")
stopifnot(file.exists(DATA))
dat <- read_stan_data(DATA)

cat(sprintf("basin contract: N = %d assemblages, K = %d classes\n",
            dat$N, dat$K))

models <- list(
  composition_gp          = here("stan", "spatial", "composition_gp.stan"),
  composition_gp_marginal = here("stan", "spatial", "composition_gp_marginal.stan")
)

L <- c(
  "# The basin spatial GP: the fit behind the reported variance share",
  "",
  sprintf("Produced by `analyses/02_spatial/05_basin_fit.R` on %s.",
          format(Sys.Date())),
  "",
  sprintf("Basin contract `data/stan/basin_composition.json`: %d assemblages, %d decorated classes.",
          dat$N, dat$K),
  "",
  "`spatial_share` is the fraction of assemblage-level compositional variance",
  "carried by the spatial field rather than by the assemblage-specific term. It",
  "is the quantity the abstract reports. The length scale `rho` is reported here",
  "for completeness and is NOT reportable as a number: see the verdict below.",
  "",
  "| model | rung | passed | verdict | spatial_share | 95% CI | rho median | rho 95% CI |",
  "|---|---|---|---|---|---|---|---|"
)
diags <- c()
results <- list()

for (mdl in names(models)) {
  cat(sprintf("\n== %s ==\n", mdl))
  sm <- stan_model(models[[mdl]])
  res <- fit_composition_gp_ladder(sm, dat, use_tau = 1, likelihood = 1,
                                   cache_key = paste0("basin_", mdl))
  results[[mdl]] <- res
  saveRDS(res, here("output", paste0("gp_real_", mdl, ".rds")))
  L <- c(L, sprintf("| %s | %d | %s | %s | %.3f | [%.3f, %.3f] | %.1f | [%.1f, %.1f] |",
                    mdl, res$rung, ifelse(isTRUE(res$passed), "yes", "**no**"),
                    res$verdict,
                    res$spatial_share[1], res$spatial_share[2], res$spatial_share[3],
                    res$rho[1], res$rho[2], res$rho[3]))
  diags <- c(diags, diag_line(res$diag, sprintf("%s (rung %d)", mdl, res$rung)))
}

L <- c(L, "", "## Diagnostics (rule 16: reported for every fit, without exception)",
       "", diags, "")

marg <- results[["composition_gp_marginal"]]
sat <- 100 * marg$diag$n_treedepth / marg$diag$n_total

L <- c(
  L, "## Reading", "",
  sprintf(paste0("The marginalised fit is the one the paper reports. It does not pass the ",
                 "rule-16 gate: treedepth saturation is %.0f percent with %d divergences. ",
                 "Under `diag_verdict` that is a **%s**, not a step-size problem, and the ",
                 "ladder therefore stops escalating rather than spending compute on a ",
                 "geometry that will not respond to it (rule 17 distinguishes an ",
                 "identifiability limit from a fixable parameterisation)."),
          sat, marg$diag$n_div, marg$verdict),
  "",
  paste0("The cause is established elsewhere rather than assumed here. ",
         "`gp_geometry_repair.md` shows the saturation is a property of the ridge and ",
         "not of the additive-versus-marginal parameterisation, and ",
         "`gp_posterior_predictive.md` shows the model reproduces cultural F_ST, the ",
         "distance decay and the spread of within-assemblage diversity, so ",
         "misspecification is excluded as the explanation. What remains is that `rho` ",
         "is flat above the design's resolving ceiling: the basin spans about 128 km ",
         "and the posterior piles up against the top of its range, which is the data ",
         "saying no decay in composition is detectable across the drainage."),
  "",
  paste0("The consequence for the manuscript is the split it already makes. ",
         "`spatial_share` is identified and is reported with its interval. A ",
         "length-scale NUMBER is not identified and is reported as unresolved, never ",
         "as an estimate. A reader is entitled to both the share and the diagnostics ",
         "above, which is why this file exists."),
  ""
)

out <- here("output", "findings", "gp_basin_fit.md")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
writeLines(L, out)
cat("\nwrote", out, "\n")
cat(paste(L[grepl("^\\||^- ", L)], collapse = "\n"), "\n")
