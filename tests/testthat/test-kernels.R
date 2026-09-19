# Re-pin the Matern constants that arrived from ../mataa (rule 19: a
# collaborator's committed result is re-derived before anything depends on it).
#
# The constants are the distances, in units of rho, at which each kernel's
# correlation falls to 0.05 (practical range) and to 0.5 (half-distance). They
# are verified here by solving k(s) = target with uniroot on the kernel written
# out from its mathematical definition in R. Nothing below reads the Stan
# implementation, so agreement is evidence and not tautology.

library(testthat)
library(here)

# Kernel correlation functions, from the definitions. Independent of Stan.
corr_fun <- function(kernel) {
  switch(kernel,
         function(s) exp(-s),
         function(s) (1 + sqrt(3) * s) * exp(-sqrt(3) * s),
         function(s) {
           x <- sqrt(5) * s
           (1 + x + x^2 / 3) * exp(-x)
         })
}

solve_const <- function(kernel, target) {
  f <- corr_fun(kernel)
  uniroot(function(s) f(s) - target, interval = c(1e-8, 50), tol = 1e-12)$root
}

# The constants under test are distances in units of rho. The roots below are
# solved in the kernel's own argument and then divided by the sqrt factor that
# relates argument to distance -- which is precisely the step the inherited
# constants omitted (output/findings/kernel_constants.md).
arg_scale <- function(kernel) c(1, sqrt(3), sqrt(5))[kernel]

test_that("half-distance constants match an independent solve", {
  # k(s) = 0.5. These are the PRIMARY reported scale, so they carry the tight
  # tolerance. The rival account -- that the constants were transcribed from a
  # different kernel parameterisation -- would show up as a mismatch of order
  # sqrt(3) or sqrt(5), far outside 1e-6.
  expect_equal(solve_const(1, 0.5), 0.6931472, tolerance = 1e-6)
  expect_equal(solve_const(2, 0.5), 0.968994, tolerance = 1e-6)
  expect_equal(solve_const(3, 0.5), 1.042122, tolerance = 1e-6)
  # The discriminating probe: the INHERITED values are the same roots before
  # the sqrt division. Assert they are wrong, so a future re-import of the
  # uncorrected file fails here loudly instead of silently.
  expect_equal(solve_const(2, 0.5) * sqrt(3), 1.678347, tolerance = 1e-6)
  expect_equal(solve_const(3, 0.5) * sqrt(5), 2.330256, tolerance = 1e-6)
})

test_that("practical-range constants match, and kernel 1 is a DOCUMENTED rounding", {
  # Kernels 2 and 3 are exact roots.
  expect_equal(solve_const(2, 0.05), 2.738871, tolerance = 1e-6)
  expect_equal(solve_const(3, 0.05), 2.646900, tolerance = 1e-6)
  expect_equal(solve_const(2, 0.05) * sqrt(3), 4.743865, tolerance = 1e-6)
  expect_equal(solve_const(3, 0.05) * sqrt(5), 5.918649, tolerance = 1e-6)
  # Kernel 1 is NOT the exact root: -log(0.05) = 2.9957, and the file uses 3
  # because exp(-3) = 0.0498. Asserting exact equality here would fail, and
  # asserting a loose tolerance everywhere would hide a real transcription
  # error in the other two. So the approximation is pinned explicitly.
  exact <- solve_const(1, 0.05)
  expect_equal(exact, 2.995732, tolerance = 1e-6)
  expect_lt(abs(3 - exact), 0.005)          # the rounding, bounded
  expect_equal(exp(-3), 0.0498, tolerance = 1e-4)
})

# The Stan-side comparison lives in test-kernels-stan.R, deliberately in its own
# file: loading rstan here makes the session hang on exit on this toolchain
# (aarch64, rstan 2.32.7), and the constants below must be runnable without it.
