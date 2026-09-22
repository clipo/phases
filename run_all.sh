#!/usr/bin/env bash
#
# run_all.sh — reproduce every analysis and figure for
#   "Are the phases real? Distinguishing bounded interaction groups from
#    spatially structured drift in lower Mississippi Valley decorated ceramics"
#
# Usage:
#   ./run_all.sh                 # uses `python3`
#   PYTHON=.venv/bin/python ./run_all.sh
#
# Prerequisites (see README.md):
#   pip install -r requirements.txt
#   pip install -e .             # makes the `mls_emergence` package importable
#
# The numbered scripts run in ascending order; outputs land in output/ and
# figures/. The run continues past any single failure and reports the count at
# the end.

set -u
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"

mkdir -p output figures
echo "Interpreter: $($PY --version 2>&1)"
echo

fail=0
run() {
  echo "=================================================================="
  echo ">>> $1"
  if ! "$PY" "$1"; then
    echo "!!! $1 exited non-zero (see message above)"
    fail=$((fail + 1))
  fi
}

# Numbered analysis + figure scripts, in order.
for s in analyses/[0-9][0-9]_*.py; do
  run "$s"
done

# Shared map helpers (make_map) and the house-style figure pipeline
# (make_figures: Figures 3, 6, 7, and S1).
run "analyses/make_map.py"
run "analyses/make_figures.py"

# Second pass for the figure scripts whose inputs are produced by a
# HIGHER-numbered script, so ascending order draws them from whatever the
# previous run left behind.
#
#   42_figS6_dynamic  reads output/tempo_posterior.npz (54),
#                     output/abc_smc_posterior.npz (38) and
#                     output/abc_smc_sbc_ranks.npz (39)
#   56_neutrality_ppc reads output/abc_pooled_posterior.npz (55)
#
# On a clean checkout the first pass either fails outright on the missing
# file or, worse, succeeds against a stale one. Re-running them here is the
# fix that keeps the numbering stable. If you add a figure script that reads
# another script's output, add it to this list or renumber it above its
# inputs.
echo
#   02_spatial/04_posterior_predictive.R reads output/gp_real_*.rds, which
#                     02_spatial/05_basin_fit.R writes (found 2026-09-22: the
#                     check had been comparing the previous run's fit with the
#                     current data on every rerun)
echo "=== second pass: figures whose inputs are produced later in the order ==="
run "analyses/42_figS6_dynamic.py"
run "analyses/56_neutrality_ppc.py"
echo "=== second pass: the spatial model's predictive check, on the fit it checks ==="
Rscript "analyses/02_spatial/04_posterior_predictive.R" || fail=$((fail + 1))

echo "=================================================================="
echo "Done. ${fail} script(s) reported a non-zero exit."
