#!/bin/bash
# Mechanical close-out after a pipeline rerun, safe to run unattended.
# Waits for the given runner PID, then refreshes generated tables, runs every
# check, and collects the values the manuscript text must be updated from into
# output/rerun/post_run_report.md. It does not edit the manuscript, commit,
# push, or sync: those need a person (or a session) to read the new values.
#   setsid nohup scripts/post_rerun_report.sh <runner-pid> > /dev/null 2>&1 < /dev/null & disown
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT" || exit 1
PY="$ROOT/.venv/bin/python"; OUT="$ROOT/output/rerun/post_run_report.md"
if [ -n "${1:-}" ]; then while kill -0 "$1" 2>/dev/null; do sleep 60; done; fi
{
echo "# Post-rerun report"; echo; echo "Written $(date -Is) by scripts/post_rerun_report.sh."; echo
L=$(grep -n "^Started" output/rerun/status.md | tail -1 | cut -d: -f1)
echo "## Pipeline"; echo; tail -n +"$L" output/rerun/status.md | grep -E "^Started|^Finished"
echo; echo "Failed steps:"; tail -n +"$L" output/rerun/status.md | awk -F'|' '/^\| [0-9]+ \|/ && $5+0 != 0'; echo
echo "## Checks"; echo '```'
"$PY" scripts/refresh_si_relaxation_tables.py 2>&1 | tail -1
"$PY" scripts/build_manifest.py 2>&1 | tail -1
"$PY" scripts/check_figure_claims.py 2>&1 | tail -8
"$PY" -m pytest -q -p no:cacheprovider tests 2>&1 | tail -1
echo '```'; echo
for f in output/revision_2026_09/report.md output/findings/scale_sweep.md output/findings/excess_locality.md \
         output/findings/excess_robustness.md output/findings/phases_as_groups.md output/findings/partition_posterior.md \
         output/findings/axis_geography.md output/findings/datum_conversion.md output/findings/minimum_sample_size.md \
         output/findings/other_departures.md output/findings/other_departures_c_highrep.md output/findings/unequal_populations.md \
         output/findings/gp_basin_fit.md output/findings/gp_posterior_predictive.md output/findings/closure_strength_posterior.md \
         output/signal_recovery.md output/findings/phase_partition_test.md output/findings/phase_recovery.md \
         output/findings/phases_under_drift.md output/findings/boundaries_in_settlement_gaps.md output/findings/settlement_and_rivers.md \
         output/findings/groupness_surface.md output/findings/mainfort_replication.md output/findings/source_effect.md \
         output/findings/partition_sensitivity.md output/findings/partition_ensemble.md output/findings/connectivity_mixing.md \
         output/findings/edge_effect.md output/findings/river_network_geometry.md; do
  echo "## $f"; echo; git diff --stat -- "$f" 2>/dev/null | tail -1; echo '```diff'
  git diff -U0 -- "$f" 2>/dev/null | grep '^[-+][^-+]' | grep -v "Produced\|2026-09-2[0-9]T" | head -60
  echo '```'; echo
done
} > "$OUT" 2>&1
