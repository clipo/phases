#!/usr/bin/env Rscript
#' Verify the C++ toolchain compiles a Stan model and NUTS recovers a known
#' posterior.
#'
#' This script is NOT optional. The program plan names the aarch64 rstan build
#' as its one platform risk; this is where the risk is retired. It fails in
#' seconds rather than forty minutes into a real fit.
#'
#' The target is computed here from the prior and the likelihood in closed form.
#' It is never read off the fit being tested (Verification Regime rule 3).
#'
#' Usage: Rscript analyses/00_setup/01_verify_toolchain.R

suppressPackageStartupMessages({
  library(rstan)
  library(here)
})

cat(strrep("=", 80), "\n")
cat("mls-emergence: Stan toolchain verification\n")
cat(strrep("=", 80), "\n")
cat("platform:    ", R.version$platform, "\n")
cat("R:           ", R.version.string, "\n")
cat("rstan:       ", as.character(packageVersion("rstan")), "\n")
cat("StanHeaders: ", as.character(packageVersion("StanHeaders")), "\n")
cat("compiler:    ", system("g++ -dumpversion", intern = TRUE), "\n\n")

SEED <- 20260831
y <- c(0.5, 1.0, 1.5)          # ybar is exactly 1.0
N <- length(y)

# Analytic posterior, from the prior and likelihood. Nothing below touches Stan.
prior_prec <- 1 / 100
post_prec  <- prior_prec + N
post_mean  <- (N * mean(y)) / post_prec
post_sd    <- 1 / sqrt(post_prec)
cat(sprintf("analytic posterior: mean = %.6f, sd = %.6f\n\n", post_mean, post_sd))

cat("compiling stan/phase0_simulation/toolchain_smoke.stan ...\n")
t0 <- Sys.time()
sm <- stan_model(here("stan", "phase0_simulation", "toolchain_smoke.stan"))
cat(sprintf("compiled in %.1f s\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

fit <- sampling(sm, data = list(N = N, y = y), seed = SEED,
                chains = 4, iter = 4000, warmup = 2000, refresh = 0)

draws <- as.matrix(fit, pars = "mu")[, 1]
est_mean <- mean(draws)
est_sd   <- sd(draws)

# Full rule-16 diagnostic set, even here. If it is not attached to the smoke
# test it will not be attached to anything.
sm_summary <- summary(fit, pars = "mu")$summary
rhat  <- sm_summary[1, "Rhat"]
ess   <- sm_summary[1, "n_eff"]
sp    <- get_sampler_params(fit, inc_warmup = FALSE)
ndiv  <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
ntot  <- sum(sapply(sp, nrow))
ntree <- sum(sapply(sp, function(x) sum(x[, "treedepth__"] >= 10)))

cat(sprintf("posterior mean: %.6f  (analytic %.6f, diff %.2e)\n",
            est_mean, post_mean, abs(est_mean - post_mean)))
cat(sprintf("posterior sd:   %.6f  (analytic %.6f, diff %.2e)\n",
            est_sd, post_sd, abs(est_sd - post_sd)))
cat(sprintf("diagnostics: R-hat %.4f | ESS %.0f | divergences %d/%d (%.2f%%) | treedepth>=10 %d\n\n",
            rhat, ess, ndiv, ntot, 100 * ndiv / ntot, ntree))

# Monte Carlo tolerance: the standard error of the mean of the retained draws,
# times a generous factor. Derived, not chosen by eye.
mcse <- est_sd / sqrt(ess)
tol  <- 10 * mcse

ok <- abs(est_mean - post_mean) < tol &&
      abs(est_sd - post_sd) < 0.02 &&
      rhat < 1.01 && ndiv == 0

if (!ok) {
  cat("FAIL: the toolchain compiled but did not recover the analytic posterior.\n")
  cat(sprintf("      tolerance was %.2e (10 x MCSE)\n", tol))
  quit(status = 1)
}

cat(strrep("=", 80), "\n")
cat("PASS. rstan compiles and samples correctly on this platform.\n")
cat(strrep("=", 80), "\n")
