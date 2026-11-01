#!/usr/bin/env Rscript
# Identified nuisance values from the cached A0 basin refit, for pinning the
# own-generator records of the A1 validation (spec 4.1c; D20).
# The intercept alone is not identified (it trades one for one with the field's
# per-class mean), so the generator takes the posterior mean of the identified
# sum alpha_c + mean_i f_ic as its intercept and generates a field centered per
# class; sigma, tau and rho are posterior medians (sigma and rho sit on the
# ridge, but the field's contrast structure they imply is identified).
# Run: Rscript analyses/probe_latent_groups/a0_pinned_nuisances.R
suppressPackageStartupMessages({ library(here); library(posterior); library(jsonlite) })
here::i_am("pyproject.toml")
rds <- list.files(here("output", "gp_cache"), pattern = "^a0_basin_.*\\.rds$", full.names = TRUE)
stopifnot(length(rds) == 1)
r <- readRDS(rds); d <- as_draws_matrix(r$fit)
N <- 28; C <- 9
fm <- sapply(1:C, function(c) rowMeans(d[, paste0("f[", 1:N, ",", c, "]")]))
a <- d[, paste0("alpha[", 1:C, "]")]
out <- list(source = basename(rds), gate_pass = r$gate$pass,
            sigma = median(d[, "sigma"]), tau = median(d[, "tau[1]"]), rho = median(d[, "rho"]),
            alpha_plus_field_mean = as.numeric(colMeans(a + fm)),
            alpha_plus_field_mean_sd = as.numeric(apply(a + fm, 2, sd)),
            alpha_alone_mean = as.numeric(colMeans(a)), alpha_alone_sd = as.numeric(apply(a, 2, sd)))
dir.create(here("output", "probe_latent_groups"), showWarnings = FALSE, recursive = TRUE)
write_json(out, here("output", "probe_latent_groups", "a0_pinned_nuisances.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE)
cat(toJSON(out[c("sigma", "tau", "rho", "alpha_plus_field_mean")], auto_unbox = TRUE, digits = 4), "\n")
