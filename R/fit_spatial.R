# Fitting layer for the spatial composition model. Sourced by analyses/02_spatial/.
#
# Rule 16 is enforced HERE, by construction: fit_composition_gp() returns
# diagnostics attached to every fit it produces, so no caller can report a
# posterior from this model without them. That is the whole reason this is a
# function rather than a call to sampling() in each script.

suppressPackageStartupMessages({ library(rstan); library(jsonlite); library(here); library(digest) })

#' Read a Stan data contract written by analyses/49_write_stan_data.py.
#' Keys beginning with "_" are provenance for humans and are dropped.
read_stan_data <- function(path) {
  # jsonlite::fromJSON treats a non-existent path as a JSON STRING and fails
  # later with a lexical error pointing at the path itself. Refuse up front
  # with a reason (rule 5); a missing contract is a missing upstream step, not
  # a parse problem.
  if (!file.exists(path)) {
    stop(sprintf(paste0("Stan data contract not found: %s\n",
                        "  Generate it with: .venv/bin/python analyses/49_write_stan_data.py"),
                 path), call. = FALSE)
  }
  d <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  d <- d[!startsWith(names(d), "_")]
  d$y <- matrix(as.integer(unlist(d$y)), nrow = d$N, byrow = FALSE)
  d$D <- matrix(as.numeric(unlist(d$D)), nrow = d$N, byrow = FALSE)
  d
}

#' Full rule-16 diagnostic set for one stanfit.
#' Returned for EVERY fit; never optional, never computed by the caller.
gp_diagnostics <- function(fit, pars = c("rho", "sigma", "tau", "spatial_share")) {
  have <- intersect(pars, fit@model_pars)
  s <- rstan::summary(fit, pars = have)$summary
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  ndiv <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  ntot <- sum(sapply(sp, nrow))
  ntree <- sum(sapply(sp, function(x) sum(x[, "treedepth__"] >= 10)))
  # E-BFMI per chain, from the definition; min across chains is what is reported.
  ebfmi <- min(sapply(sp, function(x) {
    e <- x[, "energy__"]; mean(diff(e)^2) / stats::var(e)
  }))
  dr <- rstan::extract(fit, pars = "rho", permuted = FALSE)
  list(rhat = max(s[, "Rhat"], na.rm = TRUE),
       ess_bulk = min(s[, "n_eff"], na.rm = TRUE),
       ess_tail = min(apply(dr, 3, function(m) posterior::ess_tail(m))),
       n_div = ndiv, n_total = ntot, pct_div = 100 * ndiv / ntot,
       n_treedepth = ntree, ebfmi = ebfmi)
}

#' Is a fit acceptable under rule 16 thresholds?
# Thresholds. Divergences are a VALIDITY signal and are not tolerated at all.
# Treedepth saturation is an EFFICIENCY signal in Stan's own guidance, so a
# handful is acceptable and a large fraction is not: at 100% every trajectory
# is truncated and the tail of a long-tailed posterior goes unexplored, which
# is how a weakly identified parameter comes back looking tidier than it is.
# The earlier n_treedepth == 0 rule was stricter than Stan's guidance and
# escalated fits that were mixing perfectly well (R-hat 1.004, ESS 1195).
diag_ok <- function(d) {
  d$rhat < 1.01 && d$ess_bulk >= 400 && d$n_div == 0 &&
    (d$n_treedepth / d$n_total) < 0.02 && d$ebfmi > 0.3
}

#' ESCALATION LADDER. Rungs of increasing sampler effort, tried in order until
#' the diagnostics pass. Rule 17: a model that does not converge is a defect to
#' fix, not a reason to report the number anyway. The rung reached is reported,
#' because needing rung 3 is itself information about the geometry.
# Two rungs, not three, and the reason is rule 17 rather than compute budget.
# The failure mode measured on 2026-08-31 is TREEDEPTH SATURATION WITHOUT
# DIVERGENCES, which is the signature of a long ridge, and escalating made it
# worse (25% of iterations saturating at max_treedepth 12, 46% at 14). More
# compute does not identify a parameter the data do not identify; it only buys
# longer trajectories along the same ridge. A third rung would blur the line
# between "cannot be sampled" and "is not identified", which is exactly the
# distinction rule 17 requires be drawn.
GP_LADDER <- list(
  list(iter = 3000, warmup = 1500, adapt_delta = 0.95, max_treedepth = 12),
  list(iter = 6000, warmup = 3000, adapt_delta = 0.99, max_treedepth = 13)
)

#' Classify a fit that did not pass. Distinguishes a ridge (not identified)
#' from ordinary sampling trouble, because they call for different responses.
diag_verdict <- function(d, passed) {
  if (passed) return("recovered")
  if (d$n_treedepth / d$n_total > 0.20 && d$n_div == 0) return("ridge")
  if (d$n_div > 0) return("divergent")
  "underpowered"
}

#' Fit with the ladder, returning the first rung whose diagnostics pass (or the
#' last rung tried, flagged as not passing).
fit_composition_gp_ladder <- function(sm, dat, use_tau = 1, likelihood = 1,
                                      chains = 4, seed = 20260831,
                                      ladder = GP_LADDER, cache_key = NULL,
                                      cache_dir = here::here("output", "gp_cache")) {
  # Cache on the CONTENT of the fit request, not on a label, so changing the
  # data or the ladder invalidates it automatically.
  ck <- NULL
  if (!is.null(cache_key)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    sig <- digest::digest(list(cache_key, dat, use_tau, likelihood, chains,
                               seed, ladder))
    ck <- file.path(cache_dir, paste0(sig, ".rds"))
    if (file.exists(ck)) {
      cat("    [cached]\n")
      return(readRDS(ck))
    }
  }
  last <- NULL
  for (i in seq_along(ladder)) {
    cfg <- ladder[[i]]
    res <- fit_composition_gp(sm, dat, use_tau = use_tau, likelihood = likelihood,
                              chains = chains, iter = cfg$iter, warmup = cfg$warmup,
                              adapt_delta = cfg$adapt_delta,
                              max_treedepth = cfg$max_treedepth, seed = seed)
    res$rung <- i
    res$rung_cfg <- cfg
    res$passed <- diag_ok(res$diag)
    res$verdict <- diag_verdict(res$diag, res$passed)
    last <- res
    if (res$passed) { if (!is.null(ck)) saveRDS(res, ck); return(res) }
    # A ridge does not respond to compute. Stop escalating and say so.
    if (res$verdict == "ridge" && i >= 1) {
      cat(sprintf("    rung %d: treedepth saturation %.0f%% with zero divergences -> ridge, not a step-size problem\n",
                  i, 100 * res$diag$n_treedepth / res$diag$n_total))
      if (i >= 2) { if (!is.null(ck)) saveRDS(res, ck); return(res) }
    }
    cat(sprintf("    rung %d insufficient (R-hat %.4f, bulk ESS %.0f, div %d, treedepth %d); escalating\n",
                i, res$diag$rhat, res$diag$ess_bulk, res$diag$n_div, res$diag$n_treedepth))
  }
  if (!is.null(ck)) saveRDS(last, ck)
  last
}

#' Fit composition_gp.stan and return draws plus diagnostics, always together.
fit_composition_gp <- function(sm, dat, use_tau = 1, likelihood = 1,
                               chains = 4, iter = 2000, warmup = 1000,
                               adapt_delta = 0.95, max_treedepth = 12,
                               seed = 20260831) {
  dat$use_tau <- as.integer(use_tau)
  dat$likelihood <- as.integer(likelihood)
  fit <- rstan::sampling(
    sm, data = dat, chains = chains, iter = iter, warmup = warmup, seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    pars = c("rho", "sigma", "tau", "spatial_share", "half_distance",
             "practical_range", "alpha", "y_rep"),
    refresh = 0)
  q <- function(nm) {
    if (!nm %in% fit@model_pars) return(rep(NA_real_, 3))
    x <- as.numeric(rstan::extract(fit, pars = nm)[[nm]])
    c(stats::median(x), stats::quantile(x, 0.025), stats::quantile(x, 0.975))
  }
  list(fit = fit,
       rho = q("rho"), sigma = q("sigma"), tau = q("tau"),
       spatial_share = q("spatial_share"), half_distance = q("half_distance"),
       diag = gp_diagnostics(fit))
}

#' One-line diagnostic string, for reports.
diag_line <- function(d, label) {
  sprintf(paste0("- %s: R-hat %.4f | bulk ESS %.0f | tail ESS %.0f | ",
                 "divergences %d/%d (%.2f%%) | treedepth>=10 %d | E-BFMI %.3f"),
          label, d$rhat, d$ess_bulk, d$ess_tail, d$n_div, d$n_total,
          d$pct_div, d$n_treedepth, d$ebfmi)
}
