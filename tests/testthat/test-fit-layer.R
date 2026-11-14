# Tests for R/fit_spatial.R and R/diagnostic_plots.R (probe spec 4.2; pilot
# step 3). The pure functions are checked on synthetic input with analytic
# answers; the two real fits are discriminating probes: a model that must pass
# the gate and one that must fail it for a stated reason. The real fits need
# rstan and run only with NOT_CRAN=true, as test-kernels-stan.R does:
#   NOT_CRAN=true Rscript -e \
#     'testthat::test_file("tests/testthat/test-fit-layer.R")'
library(testthat)
library(here)
source(here("R", "fit_spatial.R"))

fake_sp <- function(depths_per_chain, divergent = NULL, energy = NULL) {
  lapply(seq_along(depths_per_chain), function(c) {
    n <- length(depths_per_chain[[c]])
    cbind(accept_stat__ = rep(0.9, n), stepsize__ = rep(0.1, n),
          treedepth__ = depths_per_chain[[c]], n_leapfrog__ = 2^depths_per_chain[[c]],
          divergent__ = if (is.null(divergent)) rep(0, n) else divergent[[c]],
          energy__ = if (is.null(energy)) rnorm(n) else energy[[c]])
  })
}

test_that("treedepth saturation counts the configured maximum, not depth 10", {
  # The spec's own case: every draw at depth 10 under a maximum of 13 is zero
  # saturation. The old layer would have reported 100 percent here.
  sp <- fake_sp(list(rep(10, 500), rep(10, 500), rep(10, 500), rep(10, 500)))
  td <- treedepth_saturation(sp, 13)
  expect_equal(td$n_at_max, 0L); expect_equal(td$pct_at_max, 0); expect_equal(td$n_total, 2000L)
  expect_equal(td$max_observed, 10L)
  # Every draw at the maximum is 100 percent.
  td <- treedepth_saturation(fake_sp(list(rep(13, 100), rep(13, 100))), 13)
  expect_equal(td$pct_at_max, 100)
  # A mixed case has the exact fraction.
  td <- treedepth_saturation(fake_sp(list(c(rep(12, 30), rep(9, 70)), rep(8, 100))), 12)
  expect_equal(td$n_at_max, 30L); expect_equal(td$pct_at_max, 15)
  # A depth above the stated maximum means the maximum passed is not this fit's.
  expect_error(treedepth_saturation(fake_sp(list(rep(11, 10))), 10), "exceeds the configured maximum")
  # Missing column is refused, not treated as zero.
  bad <- list(cbind(divergent__ = rep(0, 5)))
  expect_error(treedepth_saturation(bad, 10), "treedepth__")
})

test_that("the configured maximum treedepth is read from the fit's arguments, default 10", {
  expect_equal(configured_max_treedepth(list(list(control = list(max_treedepth = 13)))), 13L)
  expect_equal(configured_max_treedepth(list(list(control = list(adapt_delta = 0.99)))), 10L)
  expect_equal(configured_max_treedepth(list(list())), 10L)
})

test_that("divergences and E-BFMI follow their definitions on synthetic input", {
  sp <- fake_sp(list(rep(5, 100), rep(5, 100)), divergent = list(c(rep(1, 3), rep(0, 97)), rep(0, 100)))
  dv <- divergences(sp)
  expect_equal(dv$n_div, 3L); expect_equal(dv$n_total, 200L); expect_equal(dv$pct_div, 1.5)
  # Energies 1..n: squared differences are all 1, variance is n(n+1)/12, so
  # E-BFMI is 12/(n(n+1)); at n = 100 that is 0.0011881 (computed by hand, not
  # through the function). An iid sequence would give about 2, so this probe
  # point separates "energy changes track its spread" from "they do not".
  n <- 100
  sp <- fake_sp(list(rep(5, n)), energy = list(as.numeric(1:n)))
  expect_equal(ebfmi_per_chain(sp), 12 / (n * (n + 1)), tolerance = 1e-12)
})

synthetic_draws <- function() {
  set.seed(7)
  it <- 400; ch <- 4
  arr <- array(NA_real_, dim = c(it, ch, 6),
               dimnames = list(NULL, NULL, c("sigma", "f[1,1]", "f[2,1]", "comembership[1,1]",
                                             "comembership[1,2]", "delta[1,1]")))
  arr[, , "sigma"] <- rnorm(it * ch)
  arr[, , "f[1,1]"] <- 0.5                                   # stuck: zero variance, not a named constant
  arr[, , "f[2,1]"] <- rnorm(it * ch)
  arr[, , "comembership[1,1]"] <- 1                          # structural constant (diagonal)
  arr[, , "comembership[1,2]"] <- 0.3 + 0.01 * rnorm(it * ch)
  for (c in 1:ch) arr[, c, "delta[1,1]"] <- rnorm(it, mean = 5 * c)   # label-sensitive, far apart by chain
  posterior::as_draws_array(arr)
}

test_that("monitored variables are classed gated, exempt, constant or stuck, and the gate reads them", {
  d <- synthetic_draws()
  s <- monitored_summary(d, gated = c("sigma", "f", "comembership"), exempt = "delta", constants = "comembership")
  cls <- setNames(s$class, s$variable)
  expect_equal(cls[["sigma"]], "gated")
  expect_equal(cls[["f[2,1]"]], "gated")
  expect_equal(cls[["f[1,1]"]], "stuck")
  expect_equal(cls[["comembership[1,1]"]], "constant")
  expect_equal(cls[["comembership[1,2]"]], "gated")
  expect_equal(cls[["delta[1,1]"]], "exempt")
  # The exempt variable has a very large R-hat (chains 5 apart) and is reported, not gated.
  expect_gt(s$rhat[s$variable == "delta[1,1]"], 1.5)
  # Base-name matching does not bleed: "f" must not match "foo".
  expect_equal(.vars_for(c("f", "f[1]", "foo", "ff[2]"), "f"), c("f", "f[1]"))
  # Assemble a diagnostics list by hand and gate it.
  gt <- s[s$class == "gated", ]
  diag <- list(rhat_max = max(gt$rhat), ess_bulk_min = min(gt$ess_bulk), ess_tail_min = min(gt$ess_tail),
               n_div = 0L, ebfmi_min = 0.9, n_stuck = sum(s$class == "stuck"), comembership_split = NULL)
  g <- gate_fit(diag)
  expect_false(g$pass); expect_equal(g$reasons, "stuck")
  diag$n_stuck <- 0L; expect_true(gate_fit(diag)$pass)
  diag$n_div <- 1L; expect_equal(gate_fit(diag)$reasons, "divergences")
  diag$n_div <- 0L; diag$ebfmi_min <- 0.2; expect_equal(gate_fit(diag)$reasons, "ebfmi")
  diag$ebfmi_min <- 0.9; diag$rhat_max <- 1.02; expect_equal(gate_fit(diag)$reasons, "rhat")
})

test_that("chain_split finds chains that disagree on co-membership and ignores the diagonal", {
  it <- 200
  arr <- array(NA_real_, dim = c(it, 2, 3),
               dimnames = list(NULL, NULL, c("comembership[1,1]", "comembership[1,2]", "comembership[2,1]")))
  arr[, , "comembership[1,1]"] <- 1
  arr[, 1, "comembership[1,2]"] <- 0.9; arr[, 2, "comembership[1,2]"] <- 0.1
  arr[, , "comembership[2,1]"] <- 0.5
  cs <- chain_split(posterior::as_draws_array(arr), "comembership")
  expect_equal(cs$max_abs_diff, 0.8, tolerance = 1e-12)
  expect_equal(cs$n_entries, 2L)
  g <- gate_fit(list(rhat_max = 1, ess_bulk_min = 1000, ess_tail_min = 1000, n_div = 0L, ebfmi_min = 1,
                     n_stuck = 0L, comembership_split = cs))
  expect_equal(g$reasons, "comembership_split")
})

test_that("the source hash covers nested includes and changes when an include changes", {
  d <- tempfile("stanmod"); dir.create(d)
  writeLines("real inner(real x) { return x; }", file.path(d, "inner.stanfunctions"))
  writeLines(c("#include inner.stanfunctions", "real outer(real x) { return inner(x); }"),
             file.path(d, "outer.stanfunctions"))
  writeLines(c("functions {", "#include outer.stanfunctions", "}", "parameters { real z; }", "model { z ~ normal(0, 1); }"),
             file.path(d, "m.stan"))
  files <- stan_include_files(file.path(d, "m.stan"))
  expect_setequal(basename(files), c("m.stan", "outer.stanfunctions", "inner.stanfunctions"))
  h1 <- stan_source_hash(file.path(d, "m.stan"))$sha256
  writeLines("real inner(real x) { return 2 * x; }", file.path(d, "inner.stanfunctions"))
  h2 <- stan_source_hash(file.path(d, "m.stan"))$sha256
  expect_false(identical(h1, h2))
  writeLines(c("#include missing.stanfunctions"), file.path(d, "bad.stan"))
  expect_error(stan_include_files(file.path(d, "bad.stan")), "include not found")
  # The project's A0 model resolves its kernel include.
  a0 <- stan_include_files(here("stan", "spatial", "composition_gp_marginal.stan"))
  expect_true(any(basename(a0) == "kernels.stanfunctions"))
})

test_that("a clean model passes the gate and its record and cache behave; a funnel fails for divergences", {
  skip_on_cran()
  cache <- tempfile("stancache"); fitcache <- tempfile("fitcache")
  set.seed(1)
  dat <- list(N = 60L, y = rnorm(60, 1, 2))
  attr(dat, "provenance") <- list(`_note` = "synthetic test data"); attr(dat, "path") <- "synthetic"
  mf <- here("tests", "testthat", "stan", "trivial_normal.stan")
  sm <- compile_model(mf, cache_dir = cache)
  expect_s4_class(sm, "stanmodel")
  expect_equal(length(list.files(cache, pattern = "\\.rds$")), 1L)
  expect_s4_class(compile_model(mf, cache_dir = cache), "stanmodel")   # cache hit, no recompile
  small <- list(list(iter = 1000L, warmup = 500L, adapt_delta = 0.9, max_treedepth = 8L))
  res <- fit_cached("trivial", mf, dat, gated = c("mu", "sigma", "lp__"), chains = 4L, seed = 11L,
                    ladder = small, cache_dir = fitcache, sm = sm, verbose = FALSE)
  expect_true(res$gate$pass)
  expect_equal(res$diag$treedepth$max_treedepth, 8L)
  expect_lte(res$diag$treedepth$max_observed, 8L)
  expect_equal(res$diag$n_gated, 3L)
  expect_false(res$cached)
  rec <- jsonlite::fromJSON(res$record_path)
  expect_true(all(c("stan_source_sha256", "fit_layer_sha256", "software", "data_provenance",
                    "sampler", "gate", "diagnostics", "variables") %in% names(rec)))
  expect_equal(rec$sampler$rung, 1L)
  expect_true(rec$gate$pass)
  # Second call with identical inputs is a cache hit; a changed seed is not.
  res2 <- fit_cached("trivial", mf, dat, gated = c("mu", "sigma", "lp__"), chains = 4L, seed = 11L,
                     ladder = small, cache_dir = fitcache, sm = sm, verbose = FALSE)
  expect_true(res2$cached)
  res3 <- fit_cached("trivial", mf, dat, gated = c("mu", "sigma", "lp__"), chains = 4L, seed = 12L,
                     ladder = small, cache_dir = fitcache, sm = sm, verbose = FALSE)
  expect_false(res3$cached)
  # Codex WB-3 (2026-10-19): a changed diagnostic policy on a cache hit is re-read from the saved draws, not
  # returned as the old gate. A threshold no finite chain can meet must fail the gate on the hit; a changed
  # draw-generation option (init) must miss the cache. The old behavior returned res$gate$pass TRUE for both.
  res4 <- fit_cached("trivial", mf, dat, gated = c("mu", "sigma", "lp__"), chains = 4L, seed = 11L,
                     ladder = small, cache_dir = fitcache, sm = sm, verbose = FALSE,
                     thresholds = modifyList(GATE, list(rhat_max = 0.5)))
  expect_true(res4$cached); expect_true(res4$policy_recomputed)
  expect_false(res4$gate$pass); expect_true("rhat" %in% res4$gate$reasons)
  expect_equal(res4$gate$thresholds$rhat_max, 0.5)
  res5 <- fit_cached("trivial", mf, dat, gated = c("mu", "sigma"), chains = 4L, seed = 11L,
                     ladder = small, cache_dir = fitcache, sm = sm, verbose = FALSE)
  expect_true(res5$cached); expect_true(res5$policy_recomputed); expect_equal(res5$diag$n_gated, 2L)
  res6 <- fit_cached("trivial", mf, dat, gated = c("mu", "sigma", "lp__"), chains = 4L, seed = 11L,
                     ladder = small, cache_dir = fitcache, sm = sm, verbose = FALSE, init = 0)
  expect_false(res6$cached)
  expect_error(fit_cached("trivial", mf, dat, gated = "mu", chains = 4L, seed = 11L, ladder = small,
                          cache_dir = fitcache, sm = sm, bogus = 1), "unknown option")
  expect_match(diag_line(res$diag, "trivial"), "treedepth at max 8:")
  # Plots are produced as PNG files.
  source(here("R", "diagnostic_plots.R"))
  figs <- tempfile("figs")
  tab <- plot_diagnostics(res$fit, pars = c("mu", "sigma"), out_dir = figs, prefix = "trivial")
  expect_setequal(tab$kind, c("rank", "trace", "pairs", "energy"))
  expect_true(all(file.exists(tab$file))); expect_true(all(file.size(tab$file) > 10000))
  expect_true(tab$font[1] %in% c("Arial", "Liberation Sans"))

  # The probe point where the gate must fail: Neal's funnel, centered, at a
  # modest adapt_delta, produces divergent transitions.
  fm <- here("tests", "testthat", "stan", "funnel.stan")
  smf <- compile_model(fm, cache_dir = cache)
  resf <- fit_with_ladder(smf, list(), gated = c("v", "x", "lp__"), chains = 4L, seed = 3L,
                          ladder = list(list(iter = 1000L, warmup = 500L, adapt_delta = 0.8, max_treedepth = 10L)),
                          verbose = FALSE)
  expect_false(resf$gate$pass)
  expect_true("divergences" %in% resf$gate$reasons)
  expect_gt(resf$diag$n_div, 0L)
})
