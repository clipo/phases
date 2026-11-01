#!/usr/bin/env Rscript
# A0 on the real basin data under three intercept priors (D25, 2026-10-18).
# Exploratory design validation (rule 1). The hard-coded Normal(0, 2) intercept
# prior left the intercepts near their prior while the identified baselines
# alpha_c + mean_i f_ic ran to -7.27 logits, so the field's constant mode
# carried 0.3 to 5.6 logits of them (discussion log, 2026-10-18). Fits:
#   P2   alpha ~ Normal(0, 2)         the inherited prior (reproduces the step-4 refit)
#   P10  alpha ~ Normal(0, 10)        wide, centred on zero, no use of the data
#   PC   alpha ~ Normal(m, 2)         centred on the pooled log-odds m_c = log(n_c / n_1)
# For each: the gate; sigma, tau, rho; the spatial share sigma^2/(sigma^2+tau^2);
# a contrast-based share (Codex technical notes, section 5) using the field's
# variance after removing its basin-wide mean, s_c = sigma^2 c(rho) with
# c(rho) = (N - sum_ij R_ij / N) / (N - 1), R = exp(-D/rho), against tau^2;
# the identified contrast SDs at 4.3, 36.6 and 79.1 km; and how much of each
# class baseline the field's mean still carries. Visual set for P10.
# Run: Rscript analyses/probe_latent_groups/a0_intercept_prior.R
suppressPackageStartupMessages({ library(here); library(rstan); library(posterior); library(jsonlite) })
here::i_am("pyproject.toml"); source(here("R", "fit_spatial.R")); source(here("R", "diagnostic_plots.R")); options(mc.cores = 4L)
A0 <- here("stan", "spatial", "composition_gp_marginal.stan"); sm <- compile_model(A0)
OUTD <- here("output", "probe_latent_groups"); FIG <- here("figures", "probe_latent_groups"); FIND <- here("output", "findings", "a0_intercept_prior.md")
GATED <- c("sigma", "tau", "rho", "alpha", "lp__", "f", "log_lik"); SEED <- 20261016L
dat <- read_stan_data(here("data", "stan", "basin_composition.json"))
pooled <- colSums(dat$y); m <- log(pooled[-1] / pooled[1])
variants <- list(P2 = list(mean = rep(0, dat$K - 1), sd = 2, label = "Normal(0, 2), inherited"),
                 P10 = list(mean = rep(0, dat$K - 1), sd = 10, label = "Normal(0, 10)"),
                 PC = list(mean = as.numeric(m), sd = 2, label = "Normal(pooled log-odds, 2)"))
N <- dat$N; D <- dat$D
qs <- function(x) c(median = stats::median(x), q5 = unname(stats::quantile(x, 0.05)), q95 = unname(stats::quantile(x, 0.95)))
fmt <- function(v) sprintf("%.3f [%.3f, %.3f]", v[1], v[2], v[3])
res <- list(); summ <- list()
for (nm in names(variants)) {
  v <- variants[[nm]]; d <- dat; d$prior_alpha_mean <- v$mean; d$prior_alpha_sd <- v$sd
  r <- fit_cached(paste0("a0_alpha_", nm), A0, d, gated = GATED, chains = 4L, seed = SEED, sm = sm)
  cat(diag_line(r$diag, paste("A0", nm)), "\n")
  dr <- as_draws_matrix(r$fit); sg <- as.numeric(dr[, "sigma"]); ta <- as.numeric(dr[, "tau[1]"]); rh <- as.numeric(dr[, "rho"])
  cfac <- vapply(rh, function(p) (N - sum(exp(-D / p)) / N) / (N - 1), numeric(1))
  fm <- sapply(1:(d$K - 1), function(c) rowMeans(dr[, paste0("f[", 1:N, ",", c, "]")]))
  al <- dr[, paste0("alpha[", 1:(d$K - 1), "]")]
  ctr <- function(dist) sqrt(2 * sg^2 * (1 - exp(-dist / rh)))
  summ[[nm]] <- list(label = v$label, gate = r$gate$pass, reasons = r$gate$reasons, rung = r$rung, seconds = r$seconds,
    diag = diag_line(r$diag, nm), td = r$diag$treedepth$pct_at_max,
    sigma = qs(sg), tau = qs(ta), rho = qs(rh), share = qs(sg^2 / (sg^2 + ta^2)),
    contrast_share = qs(sg^2 * cfac / (sg^2 * cfac + ta^2)), field_contrast_var = qs(sg^2 * cfac),
    c4 = qs(ctr(4.3)), c37 = qs(ctr(36.6)), c79 = qs(ctr(79.1)), ratio = qs(sg^2 / rh),
    alpha_mean = colMeans(al), alpha_sd = apply(al, 2, sd), field_mean = colMeans(fm), sum_mean = colMeans(al + fm), sum_sd = apply(al + fm, 2, sd),
    cor_alpha_fm = sapply(1:(d$K - 1), function(c) cor(al[, c], fm[, c])))
  res[[nm]] <- r
}
figs <- plot_diagnostics(res$P10$fit, pars = c("sigma", "rho", "tau[1]", paste0("alpha[", 1:9, "]")), out_dir = FIG, prefix = "a0_alphaP10",
                         pairs_pars = c("sigma", "rho", "tau[1]", "alpha[9]"))
row <- function(k, lab) sprintf("| %s | %s | %s | %s |", lab, fmt(summ$P2[[k]]), fmt(summ$P10[[k]]), fmt(summ$PC[[k]]))
L <- c("# A0 under three intercept priors (D25)", "",
       sprintf("Produced by `analyses/probe_latent_groups/a0_intercept_prior.R` on %s, real basin data (28 assemblages, counts sha256 `%s`). Exploratory design validation; the reported A0 fit is decided from this table.", format(Sys.Date()), attr(dat, "provenance")$`_provenance`$sha256_counts),
       "", sprintf("Pooled log-odds against the reference class, m_c: %s.", paste(sprintf("%.2f", m), collapse = ", ")), "",
       "| quantity | P2: Normal(0, 2), inherited | P10: Normal(0, 10) | PC: Normal(pooled log-odds, 2) |", "|---|---|---|---|",
       row("sigma", "sigma"), row("tau", "tau"), row("rho", "rho (km)"), row("ratio", "sigma^2/rho"),
       row("share", "spatial share sigma^2/(sigma^2+tau^2)"), row("field_contrast_var", "field variance after removing its basin mean, sigma^2 c(rho)"),
       row("contrast_share", "contrast-based share sigma^2 c(rho) / (sigma^2 c(rho) + tau^2)"),
       row("c4", "field contrast SD at 4.3 km"), row("c37", "field contrast SD at 36.6 km"), row("c79", "field contrast SD at 79.1 km"), "",
       "## How the class baselines are carried", "",
       "| class | baseline alpha + mean f (P2) | field mean, P2 | field mean, P10 | field mean, PC | corr(alpha, field mean) P2 / P10 / PC |", "|---|---|---|---|---|---|",
       sprintf("| %d | %.2f (sd %.2f) | %.2f | %.2f | %.2f | %.2f / %.2f / %.2f |", 2:dat$K, summ$P2$sum_mean, summ$P2$sum_sd, summ$P2$field_mean, summ$P10$field_mean, summ$PC$field_mean,
               summ$P2$cor_alpha_fm, summ$P10$cor_alpha_fm, summ$PC$cor_alpha_fm), "",
       "## Diagnostics", "", vapply(names(summ), function(k) sprintf("%s (rung %d, %.0f s, treedepth at max %.1f%%)", summ[[k]]$diag, summ[[k]]$rung, summ[[k]]$seconds, summ[[k]]$td), character(1)),
       "", sprintf("Visual set for P10 (%s): %s.", figs$font[1], paste(sprintf("`%s`", basename(figs$file)), collapse = ", ")))
writeLines(L, FIND); write_json(summ, file.path(OUTD, "a0_intercept_prior.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat(paste(L[7:20], collapse = "\n"), "\n\nwrote", FIND, "\n")
