#!/usr/bin/env Rscript
#' Install R dependencies for the mls-emergence Stan rebuild.
#'
#' Stan needs a C++ toolchain. gcc/g++ 13 and make are already present on this
#' box; `build-essential` is the package that supplies them on Ubuntu if they
#' are not. Nothing here needs sudo: packages install into R_LIBS_USER.
#'
#' PLATFORM NOTE. This machine is aarch64 (ARM64) Linux running R 4.3.3, which
#' is NOT the platform ../mataa builds on. rstan is compiled from source on
#' aarch64 rather than installed from a binary, so this script is the point
#' where the build risk is settled. analyses/00_setup/01_verify_toolchain.R
#' confirms the result by compiling and sampling a model with a known posterior.
#'
#' Usage: Rscript analyses/00_setup/00_install_deps.R

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Parallel compile jobs. Held at 4 rather than the 20 available cores: each
# Stan translation unit peaks around 2-4 GB and this box runs with ~13 GB free.
NCPUS <- as.integer(Sys.getenv("MLS_NCPUS", "4"))

cat(strrep("=", 80), "\n")
cat("mls-emergence: R dependency install\n")
cat(strrep("=", 80), "\n")
cat("R:        ", R.version.string, "\n")
cat("platform: ", R.version$platform, "\n")
cat("lib:      ", .libPaths()[1], "\n")
cat("Ncpus:    ", NCPUS, "\n\n")

# Ordered so the heavy Stan chain comes first and fails early if it is going to.
PKGS <- c(
  # Stan
  "StanHeaders", "rstan",
  # marginal likelihoods, the stated reason for rstan over cmdstanr
  "bridgesampling",
  # posterior handling and diagnostics (rule 16)
  "posterior", "loo", "bayesplot",
  # project plumbing
  "here", "testthat", "renv",
  # data handling. tidyverse in full pulls systemfonts/ragg, which need
  # freetype2 dev headers this box does not have; take the pieces we use.
  "dplyr", "tidyr", "readr", "readxl", "jsonlite"
)

for (p in PKGS) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("[have] %-16s %s\n", p, as.character(packageVersion(p))))
    next
  }
  cat(sprintf("[install] %s ...\n", p))
  install.packages(p, Ncpus = NCPUS)
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("\n[FAIL] %s did not install.\n", p))
    if (p == "rstan") {
      cat("rstan is the one package with a real build risk on this platform.\n",
          "Fallback, in order: (1) retry with Ncpus=1, (2) cmdstanr, which\n",
          "downloads a prebuilt CmdStan toolchain and avoids compiling the\n",
          "rstan R/C++ interface at all. Falling back to cmdstanr costs\n",
          "bridgesampling, which operates natively on stanfit objects; the\n",
          "replacement is bridge sampling on extracted draws or a grid\n",
          "quadrature reference. Record the decision in\n",
          "output/findings/toolchain_decision.md before proceeding.\n", sep = "")
    }
    quit(status = 1)
  }
  cat(sprintf("[ok] %-16s %s\n", p, as.character(packageVersion(p))))
}

cat("\n", strrep("=", 80), "\n", sep = "")
cat("All dependencies present. Next: Rscript analyses/00_setup/01_verify_toolchain.R\n")
cat(strrep("=", 80), "\n")
