# Visual diagnostics harness (shared conventions, "Bayesian Reporting"; probe
# spec 4.2). Produces, for one stanfit, the rank, trace, pairs and energy
# plots as PNG at 300 dpi and 7 inches wide, sans-serif, no titles. Called by
# committed exemplar scripts that refit at fixed seeds; never from a persisted
# fit that the script did not itself produce.

suppressPackageStartupMessages({
  library(rstan); library(bayesplot); library(ggplot2); library(posterior); library(systemfonts); library(ragg)
})

#' Arial when the machine has it, else the metric-compatible Liberation Sans.
figure_font <- function() {
  fam <- unique(systemfonts::system_fonts()$family)
  if ("Arial" %in% fam) return("Arial")
  if ("Liberation Sans" %in% fam) return("Liberation Sans")
  "sans"
}

save_png <- function(p, path, width = 7, height = 4.5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ragg::agg_png(path, width = width, height = height, units = "in", res = 300)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  invisible(path)
}

#' Rank, trace, pairs and energy plots for `fit`.
#' `pars`: variables for rank and trace; `pairs_pars`: a short list for the
#' pairs plot (divergent transitions marked). Returns a data.frame of files.
plot_diagnostics <- function(fit, pars, out_dir, prefix, pairs_pars = pars) {
  font <- figure_font()
  bayesplot::bayesplot_theme_set(bayesplot::theme_default(base_family = font, base_size = 9))
  bayesplot::color_scheme_set("viridis")
  draws <- posterior::as_draws_array(fit)
  have <- posterior::variables(draws)
  pars <- intersect(pars, have); pairs_pars <- intersect(pairs_pars, have)
  if (length(pars) == 0) stop("none of the requested variables is in the fit", call. = FALSE)
  np <- bayesplot::nuts_params(fit)
  files <- character()
  n <- length(pars); h <- max(3, 1.6 * ceiling(n / 3))
  files["rank"]  <- save_png(bayesplot::mcmc_rank_overlay(draws, pars = pars) + ggplot2::theme(legend.position = "bottom"),
                             file.path(out_dir, paste0(prefix, "_rank.png")), height = h)
  files["trace"] <- save_png(bayesplot::mcmc_trace(draws, pars = pars, np = np) + ggplot2::theme(legend.position = "bottom"),
                             file.path(out_dir, paste0(prefix, "_trace.png")), height = h)
  if (length(pairs_pars) >= 2) {
    files["pairs"] <- save_png(bayesplot::mcmc_pairs(draws, pars = pairs_pars, np = np,
                                                     off_diag_args = list(size = 0.4, alpha = 0.4)),
                               file.path(out_dir, paste0(prefix, "_pairs.png")), height = 7)
  }
  files["energy"] <- save_png(bayesplot::mcmc_nuts_energy(np) + ggplot2::theme(legend.position = "bottom"),
                              file.path(out_dir, paste0(prefix, "_energy.png")), height = 4)
  data.frame(kind = names(files), file = unname(files), font = font, dpi = 300, width_in = 7,
             stringsAsFactors = FALSE)
}
