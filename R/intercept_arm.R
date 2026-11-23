# Intercept-prior arm of the validation fits.
#
# Kept outside R/fit_spatial.R on purpose: the fit cache key hashes that file
# (layer_hash), so a helper added there would discard every cached fit.

#' The validation scripts fit their records under one of two class-intercept
#' priors: Normal(0, 10), the paper's prior and the default, or Normal(0, 2),
#' the prior-sensitivity arm, selected with the argument --alpha-sd=2. Each
#' arm has its own records (simulated at that prior's basin-fit values by
#' a0_pinned_nuisances.R and the contracts scripts) and its own findings
#' files. Returns the prior SD, the suffix of the arm's output files and the
#' directory of its records. Any other argument or value raises.
intercept_arm <- function(args = character()) {
  hit <- grepl("^--alpha-sd=", args)
  if (sum(hit) > 1 || any(!hit)) stop("unrecognised argument(s): ", paste(args, collapse = " "), "; the only option is --alpha-sd=2", call. = FALSE)
  sd <- if (any(hit)) suppressWarnings(as.numeric(sub("^--alpha-sd=", "", args[hit]))) else 10
  if (is.na(sd) || !(sd %in% c(2, 10))) stop("--alpha-sd takes 2 or 10 (the two intercept priors the paper reports)", call. = FALSE)
  if (sd == 10) return(list(sd = 10, tag = "", dir = here::here("data", "stan"), flag = ""))
  list(sd = 2, tag = "_alphaN0-2", dir = here::here("data", "stan", "alphaN0-2"), flag = " --alpha-sd=2")
}

#' Stop unless a record carries the arm's intercept prior (rule 4: the arm is
#' never applied silently to records written for the other prior).
check_arm <- function(d, arm, path) {
  if (!isTRUE(all.equal(as.numeric(d$prior_alpha_sd), arm$sd)))
    stop(sprintf("%s has prior_alpha_sd %s, but this run is the Normal(0, %g) arm", path, format(d$prior_alpha_sd), arm$sd), call. = FALSE)
  invisible(d)
}
