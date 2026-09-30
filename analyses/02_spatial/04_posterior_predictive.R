#!/usr/bin/env Rscript
#' Posterior predictive check for the spatial composition model.
#'
#' THE QUESTION THIS SETTLES. The real-basin fit puts rho far above what the
#' design can resolve and reports spatial_share = 0.96. Every diagnostic is
#' healthy except treedepth saturation, which is consistent with a genuinely
#' flat likelihood in rho above the resolving ceiling. But all of that assumes
#' the model can reproduce the observed compositions in the first place. The
#' simulated data were generated FROM this model and behave well; the real data
#' are not, and misspecification produces exactly this signature.
#'
#' If the checks below fail, no length scale and no variance share from this
#' model is reportable, whatever the diagnostics say.
#'
#' Rule 18 explicitly retains posterior predictive checks, so this is a
#' compliant instrument, not a loophole.
#'
#' THE STATISTICS ARE CHOSEN TO DISCRIMINATE. They are quantities the model was
#' NOT fitted to directly and that the paper's claims depend on:
#'   1. cultural F_ST across the three spatial clusters -- the paper's own
#'      headline quantity;
#'   2. the distance-decay of Brainerd-Robinson similarity -- the spatial
#'      signature the GP is meant to replace the Mantel test for;
#'   3. per-assemblage Gini-Simpson diversity, spread across assemblages.
#' A model that reproduces the raw counts but not these would pass a naive check
#' and still be unusable for the claims made from it.
#'
#' Usage: Rscript analyses/02_spatial/04_posterior_predictive.R

suppressPackageStartupMessages({library(rstan); library(here); library(jsonlite)})
source(here("R", "fit_spatial.R"))

dat <- read_stan_data(here("data", "stan", "basin_composition.json"))
y <- dat$y; D <- dat$D; N <- dat$N; K <- dat$K

gini <- function(v) { p <- v / sum(v); if (sum(v) == 0) NA else 1 - sum(p^2) }

fst_gs <- function(Y, lab) {
  gc <- t(sapply(sort(unique(lab)), function(c) colSums(Y[lab == c, , drop = FALSE])))
  gc <- gc[rowSums(gc) > 0, , drop = FALSE]
  if (nrow(gc) < 2) return(NA_real_)
  w <- rowSums(gc) / sum(gc)
  hs <- sum(w * apply(gc, 1, gini))
  ht <- gini(colSums(gc))
  if (is.na(ht) || ht == 0) NA_real_ else (ht - hs) / ht
}

br <- function(Y) {                       # Brainerd-Robinson similarity matrix
  P <- Y / rowSums(Y)
  S <- matrix(0, nrow(P), nrow(P))
  for (i in 1:nrow(P)) for (j in 1:nrow(P)) S[i, j] <- 200 - 100 * sum(abs(P[i, ] - P[j, ]))
  S
}
decay <- function(Y) {                    # rank correlation of similarity vs distance
  S <- br(Y); iu <- upper.tri(S)
  suppressWarnings(cor(S[iu], D[iu], method = "spearman"))
}

cl <- kmeans(cmdscale(as.dist(D), k = 2), centers = 3, nstart = 25)$cluster

stats_of <- function(Y) c(fst = fst_gs(Y, cl), decay = decay(Y),
                          div_sd = sd(apply(Y, 1, gini), na.rm = TRUE),
                          div_mean = mean(apply(Y, 1, gini), na.rm = TRUE))
obs <- stats_of(y)

L <- c("# Posterior predictive check: can the model reproduce the basin?", "",
       sprintf("Produced by `analyses/02_spatial/04_posterior_predictive.R` on %s.", format(Sys.Date())),
       "", "Statistics chosen to discriminate: quantities the model was not fitted to directly and that the paper's claims rest on.", "",
       "| model | statistic | observed | predictive median | 95% predictive interval | Bayesian p | verdict |",
       "|---|---|---|---|---|---|---|")
verdicts <- c()
for (mdl in c("composition_gp", "composition_gp_marginal")) {
  f <- here("output", paste0("gp_real_", mdl, ".rds"))
  if (!file.exists(f)) { cat("missing", f, "\n"); next }
  res <- readRDS(f)
  yr <- rstan::extract(res$fit, pars = "y_rep")$y_rep      # draws x N x K
  nd <- min(400, dim(yr)[1])
  idx <- round(seq(1, dim(yr)[1], length.out = nd))
  sim <- t(sapply(idx, function(s) stats_of(matrix(yr[s, , ], N, K))))
  for (st in colnames(sim)) {
    v <- sim[, st]; v <- v[is.finite(v)]
    pb <- mean(v >= obs[[st]])
    ok <- pb > 0.05 && pb < 0.95
    verdicts <- c(verdicts, ok)
    L <- c(L, sprintf("| %s | %s | %.4f | %.4f | [%.4f, %.4f] | %.3f | %s |",
                      mdl, st, obs[[st]], median(v),
                      quantile(v, 0.025), quantile(v, 0.975), pb,
                      ifelse(ok, "pass", "**FAIL**")))
  }
}
L <- c(L, "", "## Verdict", "",
       if (all(verdicts))
         "**The model reproduces every discriminating statistic.** Observed cultural F_ST, distance decay and diversity spread all fall inside the posterior predictive interval. Misspecification is therefore NOT the explanation for the treedepth saturation on real data, which leaves the flat-likelihood reading: rho is not identified above the design's resolving ceiling, and the posterior piling up against the top of its range is the data saying there is no decay within the basin. `spatial_share` is reportable; a length-scale NUMBER is not."
       else
         "**At least one discriminating statistic falls outside the posterior predictive interval.** The model does not reproduce the basin on a quantity the paper's claims depend on. No length scale and no variance share from this model is reportable until that is understood, whatever the convergence diagnostics say.",
       "")
out <- here("output", "findings", "gp_posterior_predictive.md")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
writeLines(L, out)
cat(paste(L[grepl("^\\||^\\*\\*", L)], collapse = "\n"), "\n")
cat("\nwrote", out, "\n")
