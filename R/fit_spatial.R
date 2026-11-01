# Fitting and diagnostics layer for the spatial composition models (A0 and A1).
#
# Rewritten 2026-10-16 (latent-group probe spec, revision 2.1, section 4.2;
# inventory row 13; pilot step 3). The previous layer (first commit 32292d9,
# 2026-09-02) counted treedepth__ >= 10 as saturation under configured maxima
# of 12 and 13, computed tail ESS for rho only, used rstan's legacy n_eff for
# bulk ESS, and did not diagnose alpha or the field; its "ridge" verdict and
# the ladder's stop rule rested on that count. None of it is kept.
#
# Rule 16 is enforced here by construction: every fit returns through
# gp_diagnostics() and gate_fit(), and fit_record() writes the diagnostics,
# the source hashes, the software versions and the data provenance beside
# every posterior. Diagnostic quantities are pure functions of draws and
# sampler parameters so that tests/testthat/test-fit-layer.R can check them
# on synthetic input with analytic answers (rule 3).
#
# Environment: native R (docs/working/2026-10-16-environment-baseline.md, section
# "Native toolchain"). From 2026-10-16 to 2026-10-20 the project ran in a
# micromamba R with R_MAKEVARS_USER=R/Makevars.conda; neither is needed now.

suppressPackageStartupMessages({
  library(rstan); library(posterior); library(jsonlite); library(here); library(digest)
})

STAN_CACHE_DIR <- here::here("output", "stan_cache")   # compiled models; gitignored
FIT_CACHE_DIR  <- here::here("output", "gp_cache")     # fits and records; gitignored

# Gate thresholds (spec 4.2). Treedepth saturation is reported, not gated,
# until the A0 refit measures it correctly and the rule is logged.
GATE <- list(rhat_max = 1.01, ess_min = 400, n_div_max = 0L, ebfmi_min = 0.3,
             comembership_split_max = 0.5)

# Two rungs (spec 4.2). Escalation happens only on a failed gate.
GP_LADDER <- list(
  list(iter = 3000L, warmup = 1500L, adapt_delta = 0.95, max_treedepth = 12L),
  list(iter = 6000L, warmup = 3000L, adapt_delta = 0.99, max_treedepth = 13L)
)

# ---------------------------------------------------------------------------
# Data contract
# ---------------------------------------------------------------------------

#' Read a Stan data contract written by analyses/49_write_stan_data.py.
#' Keys beginning with "_" are provenance (hashes, geometry) and are kept as
#' the "provenance" attribute so that fit_record() can carry them.
read_stan_data <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf(paste0("Stan data contract not found: %s\n",
                        "  Generate it with: .venv/bin/python analyses/49_write_stan_data.py"),
                 path), call. = FALSE)
  }
  d <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  prov <- d[startsWith(names(d), "_")]
  d <- d[!startsWith(names(d), "_")]
  d$y <- matrix(as.integer(unlist(d$y)), nrow = d$N, byrow = FALSE)
  d$D <- matrix(as.numeric(unlist(d$D)), nrow = d$N, byrow = FALSE)
  attr(d, "provenance") <- prov
  attr(d, "path") <- normalizePath(path)
  d
}

# ---------------------------------------------------------------------------
# Provenance: source hashes, layer hash, software versions
# ---------------------------------------------------------------------------

#' Every file a Stan model reads: the model and its #include files, resolved
#' relative to the model's own directory, recursively.
stan_include_files <- function(model_file) {
  model_file <- normalizePath(model_file, mustWork = TRUE)
  seen <- character()
  walk <- function(f) {
    if (f %in% seen) return(invisible())
    seen <<- c(seen, f)
    lines <- readLines(f, warn = FALSE)
    inc <- sub("^\\s*#include\\s+", "", grep("^\\s*#include\\s+", lines, value = TRUE))
    inc <- trimws(gsub("[\"<>]", "", inc))
    for (i in inc) {
      p <- if (startsWith(i, "/")) i else file.path(dirname(f), i)
      if (!file.exists(p)) stop(sprintf("include not found: %s (from %s)", i, f), call. = FALSE)
      walk(normalizePath(p))
    }
  }
  walk(model_file)
  seen
}

#' SHA-256 over the model source and every include, with file names, so that a
#' change to any of them changes the hash.
stan_source_hash <- function(model_file) {
  files <- stan_include_files(model_file)
  blob <- paste(vapply(files, function(f) {
    paste0(basename(f), "\n", paste(readLines(f, warn = FALSE), collapse = "\n"))
  }, character(1)), collapse = "\n----\n")
  list(files = files, sha256 = digest::digest(blob, algo = "sha256", serialize = FALSE))
}

#' Hash of this layer's own code (and the plotting harness when present).
layer_hash <- function() {
  files <- c(here::here("R", "fit_spatial.R"), here::here("R", "diagnostic_plots.R"))
  files <- files[file.exists(files)]
  blob <- paste(vapply(files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"),
                       character(1)), collapse = "\n----\n")
  digest::digest(blob, algo = "sha256", serialize = FALSE)
}

software_versions <- function() {
  v <- function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  list(R = R.version.string, platform = R.version$platform,
       rstan = v("rstan"), StanHeaders = v("StanHeaders"), posterior = v("posterior"),
       loo = v("loo"), bayesplot = v("bayesplot"), digest = v("digest"))
}

# ---------------------------------------------------------------------------
# Compilation, cached by source hash and toolchain, never beside the source
# ---------------------------------------------------------------------------

compile_model <- function(model_file, cache_dir = STAN_CACHE_DIR) {
  src <- stan_source_hash(model_file)
  ver <- software_versions()
  key <- digest::digest(list(src$sha256, ver$rstan, ver$StanHeaders, ver$platform), algo = "sha256")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  rds <- file.path(cache_dir, paste0(sub("\\.stan$", "", basename(model_file)), "_", substr(key, 1, 16), ".rds"))
  if (file.exists(rds)) {
    sm <- readRDS(rds)
    if (inherits(sm, "stanmodel")) return(sm)
  }
  sm <- rstan::stan_model(file = model_file, isystem = dirname(normalizePath(model_file)),
                          auto_write = FALSE)
  saveRDS(sm, rds)
  sm
}

# ---------------------------------------------------------------------------
# Diagnostics as pure functions
# ---------------------------------------------------------------------------

#' The maximum treedepth a fit was configured with. rstan stores the control
#' list it was called with; when no control list was passed the sampler's
#' default is 10.
configured_max_treedepth <- function(stan_args) {
  ctl <- stan_args[[1]]$control
  if (is.null(ctl) || is.null(ctl$max_treedepth)) return(10L)
  as.integer(ctl$max_treedepth)
}

#' Treedepth saturation: post-warmup iterations whose treedepth__ equals the
#' configured maximum. `sampler_params` is rstan::get_sampler_params(fit,
#' inc_warmup = FALSE): one matrix per chain with a treedepth__ column.
treedepth_saturation <- function(sampler_params, max_treedepth) {
  stopifnot(is.list(sampler_params), length(sampler_params) >= 1, is.numeric(max_treedepth))
  td <- unlist(lapply(sampler_params, function(x) {
    if (!"treedepth__" %in% colnames(x)) stop("sampler parameters lack a treedepth__ column", call. = FALSE)
    x[, "treedepth__"]
  }))
  if (any(td > max_treedepth)) {
    stop(sprintf("a treedepth of %d exceeds the configured maximum %d: the maximum passed is not this fit's",
                 max(td), max_treedepth), call. = FALSE)
  }
  n_at_max <- sum(td == max_treedepth)
  list(max_treedepth = as.integer(max_treedepth), n_at_max = as.integer(n_at_max),
       n_total = length(td), pct_at_max = 100 * n_at_max / length(td),
       max_observed = as.integer(max(td)))
}

divergences <- function(sampler_params) {
  dv <- unlist(lapply(sampler_params, function(x) x[, "divergent__"]))
  list(n_div = as.integer(sum(dv)), n_total = length(dv), pct_div = 100 * sum(dv) / length(dv))
}

#' E-BFMI per chain from the definition: mean squared energy change over the
#' energy variance.
ebfmi_per_chain <- function(sampler_params) {
  vapply(sampler_params, function(x) {
    e <- x[, "energy__"]
    mean(diff(e)^2) / stats::var(e)
  }, numeric(1))
}

#' Which draws variables belong to a base parameter name. "f" matches "f",
#' "f[1,1]", ...; it does not match "foo".
.vars_for <- function(all_vars, base_names) {
  if (length(base_names) == 0) return(character())
  pat <- paste0("^(", paste(gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", base_names), collapse = "|"), ")(\\[|$)")
  all_vars[grepl(pat, all_vars)]
}

.is_diagonal <- function(vars) {
  m <- regmatches(vars, regexec("\\[(\\d+),(\\d+)\\]$", vars))
  vapply(m, function(x) length(x) == 3 && x[2] == x[3], logical(1))
}

#' Rank-normalized R-hat, bulk and tail ESS for every monitored variable, with
#' each variable classed as gated, exempt (label-sensitive, reported only),
#' constant (a named structural constant such as the co-membership diagonal,
#' excluded) or stuck (zero variance where none is expected; a gate failure).
#' `draws` is a posterior draws_array; `gated`, `exempt`, `constants` are base
#' parameter names.
monitored_summary <- function(draws, gated, exempt = character(), constants = character()) {
  draws <- posterior::as_draws_array(draws)
  all_vars <- posterior::variables(draws)
  g <- .vars_for(all_vars, gated); e <- .vars_for(all_vars, exempt); k <- .vars_for(all_vars, constants)
  k_diag <- k[.is_diagonal(k)]          # only the diagonal of a matrix constant is structural
  k_off <- setdiff(k, k_diag)
  sel <- setdiff(unique(c(g, e, k_off)), k_diag)   # structural constants never enter the summary
  if (length(sel) == 0) stop("no monitored variable found in the draws", call. = FALSE)
  sub <- posterior::subset_draws(draws, variable = sel)
  s <- posterior::summarise_draws(sub, "mean", "sd", "rhat", "ess_bulk", "ess_tail")
  s <- as.data.frame(s)
  s$class <- ifelse(s$variable %in% e, "exempt", "gated")
  zero <- is.finite(s$sd) & s$sd == 0
  s$class[zero & s$class == "gated"] <- "stuck"
  if (length(k_diag)) {
    s <- rbind(s, data.frame(variable = k_diag, mean = NA_real_, sd = 0, rhat = NA_real_,
                             ess_bulk = NA_real_, ess_tail = NA_real_, class = "constant"))
  }
  rownames(s) <- NULL
  s
}

#' Per-chain posterior-mean matrices of a co-membership-like quantity and the
#' largest absolute off-diagonal difference between any two chains: detects
#' chains converged to different partitions (spec 4.2).
chain_split <- function(draws, base_name) {
  draws <- posterior::as_draws_array(draws)
  vars <- .vars_for(posterior::variables(draws), base_name)
  vars <- vars[!.is_diagonal(vars)]
  if (length(vars) == 0) return(list(max_abs_diff = NA_real_, n_entries = 0L))
  sub <- posterior::subset_draws(draws, variable = vars)
  a <- unclass(sub)                      # iterations x chains x variables
  cm <- apply(a, c(2, 3), mean)          # chains x variables
  if (nrow(cm) < 2) return(list(max_abs_diff = 0, n_entries = length(vars), per_chain = cm))
  d <- 0
  for (i in 1:(nrow(cm) - 1)) for (j in (i + 1):nrow(cm)) d <- max(d, max(abs(cm[i, ] - cm[j, ])))
  list(max_abs_diff = d, n_entries = length(vars), per_chain = cm)
}

#' The full diagnostic set for one stanfit (rule 16).
#' gated/exempt/constants/comembership are base parameter names (spec 4.2).
gp_diagnostics <- function(fit, gated, exempt = character(), constants = character(),
                           comembership = NULL) {
  draws <- posterior::as_draws_array(fit)
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  vars <- monitored_summary(draws, gated = gated, exempt = exempt, constants = constants)
  gt <- vars[vars$class == "gated", ]
  td <- treedepth_saturation(sp, configured_max_treedepth(fit@stan_args))
  dv <- divergences(sp)
  eb <- ebfmi_per_chain(sp)
  worst <- function(col, fun) if (nrow(gt)) gt$variable[fun(gt[[col]])] else NA_character_
  out <- list(
    variables = vars,
    n_gated = nrow(gt), n_exempt = sum(vars$class == "exempt"),
    n_constant = sum(vars$class == "constant"), n_stuck = sum(vars$class == "stuck"),
    stuck = vars$variable[vars$class == "stuck"],
    rhat_max = if (nrow(gt)) max(gt$rhat, na.rm = TRUE) else NA_real_,
    rhat_worst = worst("rhat", which.max),
    ess_bulk_min = if (nrow(gt)) min(gt$ess_bulk, na.rm = TRUE) else NA_real_,
    ess_bulk_worst = worst("ess_bulk", which.min),
    ess_tail_min = if (nrow(gt)) min(gt$ess_tail, na.rm = TRUE) else NA_real_,
    ess_tail_worst = worst("ess_tail", which.min),
    n_div = dv$n_div, n_total = dv$n_total, pct_div = dv$pct_div,
    treedepth = td,
    ebfmi = as.numeric(eb), ebfmi_min = min(eb),
    chains = length(sp), post_warmup_per_chain = nrow(sp[[1]])
  )
  if (!is.null(comembership)) out$comembership_split <- chain_split(draws, comembership)
  out
}

#' Pass or fail with reasons (spec 4.2). Treedepth is not a reason.
gate_fit <- function(d, thresholds = GATE) {
  reasons <- character()
  if (!is.finite(d$rhat_max) || d$rhat_max > thresholds$rhat_max) reasons <- c(reasons, "rhat")
  if (!is.finite(d$ess_bulk_min) || d$ess_bulk_min < thresholds$ess_min) reasons <- c(reasons, "ess_bulk")
  if (!is.finite(d$ess_tail_min) || d$ess_tail_min < thresholds$ess_min) reasons <- c(reasons, "ess_tail")
  if (d$n_div > thresholds$n_div_max) reasons <- c(reasons, "divergences")
  if (d$ebfmi_min < thresholds$ebfmi_min) reasons <- c(reasons, "ebfmi")
  if (d$n_stuck > 0) reasons <- c(reasons, "stuck")
  if (!is.null(d$comembership_split) && is.finite(d$comembership_split$max_abs_diff) &&
      d$comembership_split$max_abs_diff > thresholds$comembership_split_max) {
    reasons <- c(reasons, "comembership_split")
  }
  list(pass = length(reasons) == 0, reasons = reasons, thresholds = thresholds)
}

#' One-line diagnostic string, for reports.
diag_line <- function(d, label) {
  sprintf(paste0("- %s: R-hat max %.4f (%s) | bulk ESS min %.0f (%s) | tail ESS min %.0f (%s) | ",
                 "divergences %d/%d (%.2f%%) | treedepth at max %d: %d/%d (%.1f%%) | E-BFMI min %.3f | stuck %d"),
          label, d$rhat_max, d$rhat_worst, d$ess_bulk_min, d$ess_bulk_worst,
          d$ess_tail_min, d$ess_tail_worst, d$n_div, d$n_total, d$pct_div,
          d$treedepth$max_treedepth, d$treedepth$n_at_max, d$treedepth$n_total,
          d$treedepth$pct_at_max, d$ebfmi_min, d$n_stuck)
}

# ---------------------------------------------------------------------------
# Fitting: one rung, the ladder, the record, the cache
# ---------------------------------------------------------------------------

fit_once <- function(sm, dat, cfg, chains = 4L, seed = 20261016L, pars_save = NA, init = "random") {
  t0 <- Sys.time()
  fit <- rstan::sampling(
    sm, data = dat, chains = chains, iter = cfg$iter, warmup = cfg$warmup, seed = seed,
    control = list(adapt_delta = cfg$adapt_delta, max_treedepth = cfg$max_treedepth),
    pars = pars_save, init = init, refresh = 0)
  list(fit = fit, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

#' Fit with the ladder: the first rung whose gate passes is returned; otherwise
#' the last rung tried, with pass = FALSE and the reasons. Every rung's
#' diagnostics are kept.
fit_with_ladder <- function(sm, dat, gated, exempt = character(), constants = character(),
                            comembership = NULL, chains = 4L, seed = 20261016L,
                            ladder = GP_LADDER, pars_save = NA, init = "random",
                            thresholds = GATE, verbose = TRUE) {
  history <- list()
  for (i in seq_along(ladder)) {
    cfg <- ladder[[i]]
    r <- fit_once(sm, dat, cfg, chains = chains, seed = seed, pars_save = pars_save, init = init)
    d <- gp_diagnostics(r$fit, gated = gated, exempt = exempt, constants = constants,
                        comembership = comembership)
    g <- gate_fit(d, thresholds)
    history[[i]] <- list(rung = i, cfg = cfg, seconds = r$seconds, diag = d[setdiff(names(d), "variables")], gate = g)
    if (verbose) cat(sprintf("    rung %d (%.0f s): %s%s\n", i, r$seconds,
                             if (g$pass) "gate passed" else paste("gate failed:", paste(g$reasons, collapse = ", ")),
                             sprintf(" | treedepth at max %.1f%%", d$treedepth$pct_at_max)))
    res <- list(fit = r$fit, diag = d, gate = g, rung = i, rung_cfg = cfg, seconds = r$seconds,
                chains = chains, seed = seed, history = history)
    if (g$pass) return(res)
  }
  res
}

#' The JSON-able record that accompanies every fit (spec section 8).
fit_record <- function(res, model_file, dat, extra = list()) {
  src <- stan_source_hash(model_file)
  list(
    written = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    model_file = basename(model_file), stan_source_files = basename(src$files),
    stan_source_sha256 = src$sha256, fit_layer_sha256 = layer_hash(),
    software = software_versions(),
    data_path = attr(dat, "path"), data_provenance = attr(dat, "provenance"),
    data_shape = list(N = dat$N, K = dat$K, sherds = sum(dat$y)),
    sampler = list(chains = res$chains, seed = res$seed, rung = res$rung, rung_cfg = res$rung_cfg,
                   seconds = res$seconds),
    gate = res$gate[c("pass", "reasons", "thresholds")],
    diagnostics = res$diag[setdiff(names(res$diag), "variables")],
    variables = res$diag$variables,
    ladder_history = lapply(res$history, function(h) h[c("rung", "cfg", "seconds", "gate")]),
    extra = extra
  )
}

write_fit_record <- function(record, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(record, path, auto_unbox = TRUE, digits = NA, pretty = TRUE, null = "null", na = "null")
  invisible(path)
}

#' Cache key over everything that determines a fit (spec section 8).
cache_key <- function(model_file, dat, chains, seed, ladder, extra = list()) {
  digest::digest(list(stan_source_hash(model_file)$sha256, layer_hash(), software_versions(),
                      unclass(dat), attr(dat, "provenance"), chains, seed, ladder, extra),
                 algo = "sha256")
}

#' Fit through the cache: a hit returns the saved draws with the diagnostics and
#' gate recomputed under the policy of the present call; a miss fits, writes
#' both, and returns them. `label` names the files.
#'
#' Two kinds of option travel in `...` (Codex WB-3, 2026-10-19). Draw
#' generation: `pars_save` and `init` change the draws, so they enter the key
#' when supplied (calls that leave them at their defaults keep their keys).
#' Diagnostic policy: `gated`, `exempt`, `constants`, `comembership` and
#' `thresholds` change only how the draws are read, so a hit recomputes the
#' diagnostics and gate from the saved fit under the requested policy whenever
#' it differs from the one the fit was saved with (or the saved result predates
#' this bookkeeping). The returned result says so in `policy_recomputed`; the
#' record on disk keeps the gate of the original call.
DRAW_OPTIONS <- c("pars_save", "init")
POLICY_OPTIONS <- c("exempt", "constants", "comembership", "thresholds")
fit_policy <- function(gated, opts) {
  list(gated = gated, exempt = opts$exempt %||% character(), constants = opts$constants %||% character(),
       comembership = opts$comembership, thresholds = opts$thresholds %||% GATE)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
fit_cached <- function(label, model_file, dat, gated, ..., chains = 4L, seed = 20261016L,
                       ladder = GP_LADDER, extra = list(), cache_dir = FIT_CACHE_DIR, sm = NULL) {
  opts <- list(...)
  unknown <- setdiff(names(opts), c(DRAW_OPTIONS, POLICY_OPTIONS, "verbose"))
  if (length(unknown)) stop("fit_cached: unknown option(s) ", paste(unknown, collapse = ", "))
  draw_opts <- opts[intersect(names(opts), DRAW_OPTIONS)]
  key <- cache_key(model_file, dat, chains, seed, ladder, c(extra, if (length(draw_opts)) list(draw_options = draw_opts)))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(cache_dir, paste0(label, "_", substr(key, 1, 16)))
  policy <- fit_policy(gated, opts)
  if (file.exists(paste0(stem, ".rds")) && file.exists(paste0(stem, ".json"))) {
    res <- readRDS(paste0(stem, ".rds")); res$cached <- TRUE; res$record_path <- paste0(stem, ".json")
    res$policy_recomputed <- FALSE
    if (is.null(res$policy) || !identical(res$policy, policy)) {
      d <- gp_diagnostics(res$fit, gated = policy$gated, exempt = policy$exempt, constants = policy$constants,
                          comembership = policy$comembership)
      res$diag <- d; res$gate <- gate_fit(d, policy$thresholds); res$policy <- policy; res$policy_recomputed <- TRUE
    }
    return(res)
  }
  if (is.null(sm)) sm <- compile_model(model_file)
  res <- fit_with_ladder(sm, dat, gated = gated, ..., chains = chains, seed = seed, ladder = ladder)
  res$policy <- policy
  rec <- fit_record(res, model_file, dat, extra = c(extra, list(cache_key = key, label = label)))
  write_fit_record(rec, paste0(stem, ".json"))
  res$record <- rec; res$cached <- FALSE; res$record_path <- paste0(stem, ".json"); res$policy_recomputed <- FALSE
  saveRDS(res, paste0(stem, ".rds"))
  res
}
