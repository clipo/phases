# Does kernels.stanfunctions' MATRIX code match the R kernels it was copied
# from? Separated from test-kernels.R because it needs rstan, and rstan makes
# this session hang on exit on this toolchain.
#
# STATUS: skipped. See the skip() reason below.

library(testthat)
library(rstan)
library(here)

corr_fun <- function(kernel) {
  switch(kernel,
         function(s) exp(-s),
         function(s) (1 + sqrt(3) * s) * exp(-sqrt(3) * s),
         function(s) { x <- sqrt(5) * s; (1 + x + x^2 / 3) * exp(-x) })
}

test_that("the Stan implementation agrees with the R kernel it was copied from", {
  # This is the part that touches Stan. Distances and rho chosen so the
  # correlation spans the full range rather than clustering near 1.
  skip_on_cran()
  fns <- tempfile(fileext = ".stan")
  writeLines(c("functions {",
               paste0("#include ", here("stan", "spatial", "kernels.stanfunctions")),
               "}", "model {}"), fns)
  expose_stan_functions(stanc(file = fns))

  D <- matrix(c(0, 5, 20, 5, 0, 12, 20, 12, 0), 3, 3)
  for (k in 1:3) {
    K <- k_matern(D, 1.0, 10.0, 3L, 0.0, as.integer(k))
    f <- corr_fun(k)
    expect_equal(K[1, 2], f(5 / 10), tolerance = 1e-10)
    expect_equal(K[1, 3], f(20 / 10), tolerance = 1e-10)
    expect_equal(K[2, 3], f(12 / 10), tolerance = 1e-10)
    expect_equal(K[1, 1], 1.0, tolerance = 1e-12)   # sigma^2 on the diagonal
  }
  # The constants, as corrected, must actually deliver their stated correlation.
  # This is the assertion that would have caught the inherited error: it checks
  # the CONSEQUENCE (correlation at the quoted distance) rather than the value.
  for (k in 1:3) {
    rc <- range_const(as.integer(k)); hc <- half_const(as.integer(k))
    f <- corr_fun(k)
    expect_equal(f(rc), 0.05, tolerance = 1e-3,
                 info = sprintf("range_const(%d) must give correlation 0.05", k))
    expect_equal(f(hc), 0.50, tolerance = 1e-5,
                 info = sprintf("half_const(%d) must give correlation 0.50", k))
  }
})
