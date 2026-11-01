#!/usr/bin/env Rscript
# A1 validated on its own generator (spec 4.1c, 6.6; R2-1). Pilot step 6.
# Exploratory design validation (rule 1). Three parts, every fit through the
# diagnostics layer and the content-keyed cache:
#   (i)  code equivalence at zero group scale: A0 and A1 with fix_group_scale = 1
#        on the same zero-shift record must agree on the shared parameters and
#        the pointwise log-likelihood within Monte Carlo error;
#   (ii) null calibration: unconstrained A1 on the zero-shift record; Delta,
#        Delta_w, coherence and the separation margins against the generating
#        labels are the measured "what no groups looks like" reference;
#   (iii) recovery: A1 on five records at sigma_g = 0.25 and five at 0.5, read
#        as Delta, Delta_w, margins and sigma_g against the truth and each
#        record's realized between-group contrast.
# Twenty draws of one fit are exported for the Stan-versus-numpy readout
# cross-check (a1_readout_crosscheck.py).
# Run: Rscript analyses/probe_latent_groups/validate_a1.R
suppressPackageStartupMessages({ library(here); library(rstan); library(posterior); library(jsonlite) })
here::i_am("pyproject.toml")
source(here("R", "fit_spatial.R")); source(here("R", "diagnostic_plots.R"))
options(mc.cores = 4L)
A0 <- here("stan", "spatial", "composition_gp_marginal.stan")
A1 <- here("stan", "spatial", "composition_gp_groups.stan")
OUTD <- here("output", "probe_latent_groups"); FIG <- here("figures", "probe_latent_groups")
FIND <- here("output", "findings", "a1_validation.md")
SEED <- 20261017L
GATED_A0 <- c("sigma", "tau", "rho", "alpha", "lp__", "f", "log_lik")
GATED_A1 <- c(GATED_A0, "sigma_g", "Delta", "Delta_w", "coherence", "comembership")
EXEMPT_A1 <- c("delta", "pi", "resp", "z_g", "sigma_g_free", "z_draw", "coherence_hard", "y_rep")
A1_KEYS <- c("G", "prior_group_sd", "fix_group_scale", "group_scale_fixed", "M", "edges")
sm0 <- compile_model(A0); sm1 <- compile_model(A1)
read_a1 <- function(f) { d <- read_stan_data(f); d$edges <- matrix(as.integer(unlist(d$edges)), nrow = d$M, byrow = FALSE); d }
strip_a1 <- function(d) { p <- attr(d, "provenance"); pa <- attr(d, "path"); d <- d[setdiff(names(d), A1_KEYS)]; attr(d, "provenance") <- p; attr(d, "path") <- pa; d }
post <- function(fit, pars) {
  d <- posterior::subset_draws(posterior::as_draws_array(fit), variable = pars)
  as.data.frame(posterior::summarise_draws(d, "mean", "median", "mcse_mean", ~posterior::quantile2(.x, probs = c(0.05, 0.95))))
}
margins <- function(fit, z) {
  q <- apply(rstan::extract(fit, pars = "comembership")$comembership, c(2, 3), mean)
  g <- sort(unique(z)); B <- matrix(NA_real_, length(g), length(g))
  for (a in seq_along(g)) for (b in seq_along(g)) {
    ia <- which(z == g[a]); ib <- which(z == g[b])
    B[a, b] <- if (a == b) mean(q[ia, ia][upper.tri(q[ia, ia])]) else mean(q[ia, ib])
  }
  list(block = B, margin = sapply(seq_along(g), function(a) B[a, a] - max(B[a, -a])))
}
timing <- list(); note <- function(l, r) timing[[l]] <<- list(rung = r$rung, seconds = r$seconds, cached = isTRUE(r$cached), pass = r$gate$pass, reasons = r$gate$reasons)

# ---------------------------------------------------------------------------
cat("\n[i] code equivalence at zero group scale\n")
dz <- read_a1(here("data", "stan", "a1_zero.json")); truth0 <- attr(dz, "provenance")$`_truth`
r0 <- fit_cached("a1val_a0_on_zero", A0, strip_a1(dz), gated = GATED_A0, chains = 4L, seed = SEED, sm = sm0); note("A0 on zero", r0)
dzf <- dz; dzf$fix_group_scale <- 1L; dzf$group_scale_fixed <- 0
r1f <- fit_cached("a1val_a1_fixed0_on_zero", A1, dzf, gated = GATED_A0,
                  exempt = c(EXEMPT_A1, "sigma_g", "Delta", "Delta_w", "coherence", "comembership"),
                  chains = 4L, seed = SEED, sm = sm1); note("A1 fixed 0 on zero", r1f)
shared <- c("sigma", "tau[1]", "rho", paste0("alpha[", 1:9, "]"), "spatial_share", paste0("log_lik[", 1:28, "]"))
p0 <- post(r0$fit, shared); p1 <- post(r1f$fit, shared)
zsc <- (p1$mean - p0$mean) / sqrt(p0$mcse_mean^2 + p1$mcse_mean^2)
equiv <- data.frame(variable = p0$variable, mean_A0 = p0$mean, mean_A1fixed = p1$mean, z = zsc)
cat(sprintf("  max |z| over %d shared quantities: %.2f (%s)\n", nrow(equiv), max(abs(zsc)), equiv$variable[which.max(abs(zsc))]))

# ---------------------------------------------------------------------------
cat("\n[ii] null calibration: unconstrained A1 on the zero-shift record\n")
r1 <- fit_cached("a1val_a1_on_zero", A1, dz, gated = GATED_A1, exempt = EXEMPT_A1, constants = "comembership",
                 comembership = "comembership", chains = 4L, seed = SEED, sm = sm1); note("A1 on zero", r1)
cat(diag_line(r1$diag, "A1 on zero"), "\n")
null_post <- post(r1$fit, c("sigma_g", "Delta", "Delta_w", "coherence", "coherence_hard", "sigma", "tau[1]", "rho"))
null_marg <- margins(r1$fit, truth0$z)
figs <- plot_diagnostics(r1$fit, pars = c("sigma", "rho", "tau[1]", "sigma_g", "Delta", "coherence", "alpha[1]", "alpha[2]"),
                         out_dir = FIG, prefix = "a1_null", pairs_pars = c("sigma", "rho", "tau[1]", "sigma_g"))
# Export 20 draws for the numpy cross-check.
ex <- rstan::extract(r1$fit, pars = c("delta", "pi", "resp", "z_draw", "Delta", "Delta_w", "coherence", "coherence_hard"))
idx <- round(seq(1, length(ex$Delta), length.out = 20))
write_json(list(edges = dz$edges, draws = lapply(idx, function(s) list(
  delta = ex$delta[s, , ], pi = ex$pi[s, ], resp = ex$resp[s, , ], z_draw = ex$z_draw[s, ],
  Delta = ex$Delta[s], Delta_w = ex$Delta_w[s], coherence = ex$coherence[s], coherence_hard = ex$coherence_hard[s]))),
  file.path(OUTD, "a1_crosscheck_draws.json"), digits = NA, auto_unbox = TRUE)

# ---------------------------------------------------------------------------
cat("\n[iii] recovery on ten nonzero records\n")
rec_files <- sort(list.files(here("data", "stan"), pattern = "^a1_sg.*\\.json$", full.names = TRUE))
recovery <- list()
for (f in rec_files) {
  tag <- sub("^a1_(.*)\\.json$", "\\1", basename(f)); d <- read_a1(f); tr <- attr(d, "provenance")$`_truth`
  r <- fit_cached(paste0("a1val_", tag), A1, d, gated = GATED_A1, exempt = EXEMPT_A1, constants = "comembership",
                  comembership = "comembership", chains = 4L, seed = SEED, sm = sm1); note(tag, r)
  cat(diag_line(r$diag, tag), "\n")
  pt <- post(r$fit, c("sigma_g", "Delta", "Delta_w", "coherence", "coherence_hard"))
  mg <- margins(r$fit, tr$z)
  recovery[[tag]] <- list(truth_sigma_g = tr$sigma_g, Delta_true = tr$Delta_true, realized_fst = tr$realized_fst_between_true_groups,
                          post = pt, margin = mg$margin, block = mg$block, gate = r$gate, split = r$diag$comembership_split$max_abs_diff)
}

# ---------------------------------------------------------------------------
cell <- function(p, v) { i <- which(p$variable == v); sprintf("%.3f [%.3f, %.3f]", p$median[i], p$q5[i], p$q95[i]) }
L <- c("# A1 validated on its own generator (pilot step 6; spec 4.1c)", "",
       sprintf("Produced by `analyses/probe_latent_groups/validate_a1.R` on %s. Records from `a1_contracts.py` (nuisances pinned to the A0 refit's identified quantities; field centered per class; membership = Phillips's partition 12/11/5; Gabriel edges for coherence). Exploratory design validation, not evidence about the record.", format(Sys.Date())),
       "", "## (i) Code equivalence at zero group scale", "",
       sprintf("A0 and A1 with the group scale fixed at zero, same record, same seed. Largest |z| over %d shared posterior means (difference over combined MCSE): **%.2f** at `%s`. Values above 3 would mean the two programs disagree.", nrow(equiv), max(abs(zsc)), equiv$variable[which.max(abs(zsc))]),
       "", "| quantity | A0 mean | A1 (fixed 0) mean | z |", "|---|---|---|---|",
       sprintf("| %s | %.4f | %.4f | %.2f |", equiv$variable[1:13], equiv$mean_A0[1:13], equiv$mean_A1fixed[1:13], equiv$z[1:13]),
       sprintf("| log_lik[1..28] | | | max abs z %.2f |", max(abs(equiv$z[14:41]))), "",
       diag_line(r0$diag, "A0 on zero"), diag_line(r1f$diag, "A1 fixed 0 on zero"), "",
       "## (ii) Null calibration: unconstrained A1 on the zero-shift record", "",
       sprintf("Truth: sigma_g = 0, no groups; the generating labels (Phillips's partition) are used only to read the margins. Realized between-'group' F_ST of the record from the field alone: %.4f (the record's is 0.0123; see the A0 adequacy note below).", truth0$realized_fst_between_true_groups),
       "", "| quantity | median [5%, 95%] |", "|---|---|",
       sprintf("| %s | %s |", c("sigma_g", "Delta", "Delta_w", "coherence (soft)", "coherence (hard)"), sapply(c("sigma_g", "Delta", "Delta_w", "coherence", "coherence_hard"), function(v) cell(null_post, v))),
       sprintf("| separation margins vs generating labels (12, 11, 5) | %s |", paste(sprintf("%.3f", null_marg$margin), collapse = ", ")),
       sprintf("| per-chain co-membership split (max abs diff) | %.3f |", r1$diag$comembership_split$max_abs_diff), "",
       diag_line(r1$diag, "A1 on zero"),
       sprintf("Visual set: %s.", paste(sprintf("`%s`", basename(figs$file)), collapse = ", ")), "",
       "## (iii) Recovery against the prior (sigma_g = 0.25 and 0.5, five records each)", "",
       "| record | true sigma_g | realized between-group F_ST | Delta true | sigma_g post | Delta post | Delta_w post | coherence | margins (12, 11, 5) | gate | chain split |",
       "|---|---|---|---|---|---|---|---|---|---|---|",
       sapply(names(recovery), function(tg) { x <- recovery[[tg]]
         sprintf("| %s | %.2f | %.4f | %.3f | %s | %s | %s | %s | %s | %s | %.3f |", tg, x$truth_sigma_g, x$realized_fst, x$Delta_true,
                 cell(x$post, "sigma_g"), cell(x$post, "Delta"), cell(x$post, "Delta_w"), cell(x$post, "coherence"),
                 paste(sprintf("%.2f", x$margin), collapse = ", "), if (x$gate$pass) "pass" else paste("fail:", paste(x$gate$reasons, collapse = ",")), x$split) }),
       "", "## A0 adequacy note (found while pinning the generator, 2026-10-17)", "",
       "The A0 posterior predictive with the fitted field reproduces the record's between-area F_ST (0.0122 against 0.0123), but fresh field draws at the same hyperparameters give a median of 0.085 (10th to 90th percentile 0.02 to 0.26). Between Phillips's Parkin and Kent areas the pooled logit share differs by +0.47 for Parkin_Punctated, -0.12 for Barton/Kent/MPI, -0.90 for Painted, and from -2.2 to +0.83 [corrected 2026-10-19 after the Codex audit of the regional-drift spec: Fortune_Noded is +0.83 in area_logits.md; first written as -0.17 to -2.2, which omitted it] for the rare classes, two of which (Walls_Engraved, Hull_Engraved) are absent from every Parkin-area site; per-class SDs of the logit share across sites are 0.8 to 1.8 (`analyses/probe_latent_groups/area_logits.py`, `output/findings/area_logits.md`). The reference class is not the flattest on this scale. A0 has one field scale for all classes and holds the record's pattern only in the realized field. A per-class field scale is the extension this points to. The own-generator records above therefore carry far more between-area contrast than the record (realized F_ST column), which is acceptable for a test of the software and of identification under the model (rule 20c) but is not a resemblance claim.",
       "", "## Timing and gates", "", "| fit | rung | seconds | gate | reasons | cached |", "|---|---|---|---|---|---|",
       sapply(names(timing), function(nm) { t <- timing[[nm]]; sprintf("| %s | %d | %.0f | %s | %s | %s |", nm, t$rung, t$seconds, if (t$pass) "pass" else "fail", paste(t$reasons, collapse = ", "), if (t$cached) "yes" else "no") }))
writeLines(L, FIND)
write_json(list(equivalence = equiv, null = list(post = null_post, margin = null_marg$margin, block = null_marg$block),
                recovery = lapply(recovery, function(x) x[c("truth_sigma_g", "Delta_true", "realized_fst", "post", "margin", "block", "split")]),
                timing = timing), file.path(OUTD, "a1_validation_summary.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat("\nwrote", FIND, "\n")
