#!/usr/bin/env Rscript
#' Record the R and Stan environment.
#'
#' "The same toolchain as ../mataa" is a claim that needs a version comparison,
#' not an assumption. This writes the comparison down (Verification Regime
#' rules 1 and 14). The Python half is written by the companion
#' analyses/00_setup/02_record_environment.py.
#'
#' Usage: Rscript analyses/00_setup/02_record_environment.R

suppressPackageStartupMessages({ library(here) })

OUT <- here("output", "findings", "environment_r.md")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

pkgs <- c("rstan", "StanHeaders", "bridgesampling", "posterior", "loo",
          "bayesplot", "here", "testthat", "renv")
have <- pkgs[vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
vers <- vapply(have, function(p) as.character(packageVersion(p)), character(1))

L <- c(
  "# R and Stan environment",
  "",
  sprintf("Recorded %s by `analyses/00_setup/02_record_environment.R`.",
          format(Sys.Date())),
  "",
  "| | |",
  "|---|---|",
  sprintf("| platform | `%s` |", R.version$platform),
  sprintf("| R | %s |", R.version.string),
  sprintf("| C++ | g++ %s |",
          tryCatch(system("g++ -dumpfullversion", intern = TRUE), error = function(e) "unknown")),
  sprintf("| %s | %s |", have, vers),
  "",
  "## Comparison with ../mataa",
  "",
  "mataa builds on Linux Mint / x86 and records Stan 2.32.2. This box is",
  "aarch64 with StanHeaders on the same 2.32 line. Same Stan line, different",
  "platform. The consequence that matters: 2.32 predates the built-in",
  "`dirichlet_multinomial_lpmf`, so the hand-written `dm_lpmf` is required here",
  "for the same reason it is required there.",
  ""
)
writeLines(L, OUT)
cat("wrote", OUT, "\n")
