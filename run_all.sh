#!/usr/bin/env bash
#
# run_all.sh — reproduce the figures and analyses of
#   "What Are Phases? Geographic Structure, Neutral Copying, and the Search
#    for Social Boundaries in the Lower Mississippi Valley"
#
# Usage:
#   ./run_all.sh              # level 1: Figures 1 to 7 and S1 from the committed findings files (minutes)
#   ./run_all.sh --analyses   # level 2: also rerun the Python analyses that write those findings files (hours)
#   ./run_all.sh --stan       # level 3: also refit the Stan models (R, rstan; hours to days)
#   PYTHON=.venv/bin/python ./run_all.sh
#
# Prerequisites (see README.md):
#   pip install -r requirements.txt && pip install -e .
#   for --stan: R, a C++ toolchain, and Rscript analyses/00_setup/00_install_deps.R
#
# Each level reads what the one before it committed. The run continues past
# any single failure and reports the count at the end. MANIFEST.md names the
# script that writes every output.
#
# Not rerunnable here: analyses/rims/rim_r1.py (Figure 6's findings file,
# output/findings/rim_r1.md) reads an unpublished rim tabulation that is not
# in the release; see data/README.md. The figure is redrawn from the committed
# findings file in level 1.

set -u
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
LEVEL=1
for a in "$@"; do
  case "$a" in
    --analyses) LEVEL=2 ;;
    --stan)     LEVEL=3 ;;
    *) echo "unknown option: $a"; exit 2 ;;
  esac
done

mkdir -p output/findings output/probe_latent_groups figures
echo "Interpreter: $($PY --version 2>&1)"
echo "Level: $LEVEL"
echo

fail=0
run() {
  echo "=================================================================="
  echo ">>> $*"
  if ! "$@"; then
    echo "!!! $1 exited non-zero (see message above)"
    fail=$((fail + 1))
  fi
}
py() { run "$PY" "$@"; }

if [ "$LEVEL" -ge 2 ]; then
  echo "=== level 2: the data products and the Python analyses ==="
  py scripts/build_analysis_matrix.py
  py analyses/16_basin_membership.py
  py analyses/82_minimum_sample_size.py                       # Table S1
  py analyses/76_phase_recovery.py                            # Figure 3A, Table S2
  py analyses/79_phases_under_drift.py                        # Table S2
  py analyses/80_boundaries_in_settlement_gaps.py             # Figure 3B
  py analyses/81_settlement_and_rivers.py                     # Table S3
  py analyses/11_chronology_14c.py                            # Figure S1
  py scripts/build_mainfort2003_sites.py
  py analyses/83_mainfort_replication.py                      # Tables S7, S8
  py analyses/49_write_stan_data.py                           # data/stan/*.json
  # The neutral copying model: the screen, its second stage and the per-class
  # comparison (Table 1's neutral comparisons), then regional drift and the
  # planted markers (Figure 5, Tables 2 and S6), then composition by collection.
  py analyses/probe_latent_groups/resemblance_screen.py
  py analyses/probe_latent_groups/resemblance_stage2.py
  py analyses/probe_latent_groups/per_class_resemblance.py
  py analyses/probe_latent_groups/site_fst_decomposition.py
  py analyses/per_type/per_class_table.py --generator C3 --settings 0,0.01,0.03,0.1,0.3 --records 50 --points primary
  py analyses/per_type/per_class_table.py --generator C3 --settings 0,0.01,0.03,0.1,0.3 --records 50 --points secondary
  for cls in "Parkin_Punctated" "Barton/Kent/MPI" "Painted" "Walls_Engraved"; do
    for cl in small large; do
      py analyses/per_type/per_class_table.py --generator M --settings 0/1.5,0/2,0/3 --marker-class "$cls" --marker-cluster "$cl" --records 50 --points primary
      py analyses/per_type/per_class_table.py --generator M --settings 0/1.5,0/2,0/3 --marker-class "$cls" --marker-cluster "$cl" --records 50 --points secondary
    done
  done
  py analyses/per_type/drift_sensitivity_figure.py
  py analyses/per_type/composition_by_collection.py
fi

if [ "$LEVEL" -ge 3 ]; then
  echo "=== level 3: the Stan fits (every fit passes through R/fit_spatial.R) ==="
  run Rscript analyses/probe_latent_groups/refit_a0_basin.R           # A0 on the basin; Table S4
  run Rscript analyses/probe_latent_groups/a0_intercept_prior.R       # Figure 4
  run Rscript analyses/probe_latent_groups/validate_a1.R              # Figure 7A
  run Rscript analyses/probe_latent_groups/validate_a3.R              # Figure 7B
  run Rscript analyses/probe_latent_groups/a1_prior_check.R
  run Rscript analyses/probe_latent_groups/a3_pinned_range.R 6        # Table S5
  run Rscript analyses/probe_latent_groups/a3_pinned_range.R 12       # Table S5
  py analyses/probe_latent_groups/a0_predictive_summaries.py
  py analyses/probe_latent_groups/a3_prior_comparison.py
fi

echo "=== level 1: the figures, from the committed findings files ==="
py analyses/93_phillips1970_phase_map.py           # Figure 1
py analyses/102_fig2_two_histories.py              # Figure 2
py analyses/103_fig3_phase_recovery_map.py         # Figure 3
py analyses/104_fig4_spatial_field.py              # Figure 4
py analyses/105_fig5_per_class_drift.py            # Figure 5
py analyses/106_fig6_rim_attributes.py             # Figure 6
py analyses/107_fig7_latent_unit_recovery.py       # Figure 7
py analyses/11_chronology_14c.py                   # Figure S1 (figures/figS6_chronology)
# Figure S2 (figures/fig8_ranksize) is drawn by analyses/make_figures.py,
# which also redraws the earlier version's house-style figures:
#   $PY analyses/make_figures.py

echo "=================================================================="
echo "Done. ${fail} script(s) reported a non-zero exit."
