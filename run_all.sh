#!/usr/bin/env bash
#
# run_all.sh — reproduce the figures and analyses of
#   "What Are Phases? Geographic Structure, Neutral Copying, and the Search
#    for Social Boundaries in the Lower Mississippi Valley"
#
# Usage:
#   ./run_all.sh              # level 1: Figures 1 to 8 and S1 from the committed findings files (minutes)
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
# Level 2 first measures the simulator time scales
# (analyses/probe_latent_groups/time_scales.py), which the neutral copying
# model's scripts and the rim analysis read from output/probe_latent_groups/.
# That directory is not published, so level 2 must run before any of those
# scripts is run alone.

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
  py analyses/76_phase_recovery.py                            # Figure 4A, Table S2
  py analyses/79_phases_under_drift.py                        # Table S2
  py analyses/80_boundaries_in_settlement_gaps.py             # Figure 4B
  py analyses/81_settlement_and_rivers.py                     # Table S3
  py analyses/11_chronology_14c.py                            # Figure S1
  py scripts/build_mainfort2003_sites.py
  py analyses/83_mainfort_replication.py                      # Tables S9, S10
  py analyses/49_write_stan_data.py                           # data/stan/*.json
  # The neutral copying model: the screen, its second stage and the per-class
  # comparison (Table 1), then regional drift and the planted markers
  # (Figure 6, Table S6), then composition by collection and the painted class
  # by type.
  py analyses/probe_latent_groups/time_scales.py                    # simulator time scales, read by the scripts below
  py analyses/probe_latent_groups/resemblance_screen.py
  py analyses/probe_latent_groups/resemblance_stage2.py
  py analyses/probe_latent_groups/per_class_resemblance.py
  py analyses/probe_latent_groups/site_fst_decomposition.py
  py analyses/per_type/per_class_table.py --generator C3 --settings 0,0.01,0.03,0.1,0.3 --records 50 --points primary
  py analyses/per_type/per_class_table.py --generator C3 --settings 0,0.01,0.03,0.1,0.3 --records 50 --points secondary
  for cls in "Parkin_Punctated" "Barton/Kent/MPI" "Painted" "Walls_Engraved"; do
    for cl in small large; do
      # The settings are the ones the committed tables record in their own headers: the full amplitude-by-strength
      # grid at the primary point, and amplitude 0.1 only at the secondary point.
      py analyses/per_type/per_class_table.py --generator M --settings 0/1.5,0/2,0/3,0.01/1.5,0.01/2,0.01/3,0.03/1.5,0.03/2,0.03/3,0.1/1.5,0.1/2,0.1/3,0.3/1.5,0.3/2,0.3/3 --marker-class "$cls" --marker-cluster "$cl" --records 50 --points primary
      py analyses/per_type/per_class_table.py --generator M --settings 0.1/1.5,0.1/2,0.1/3 --marker-class "$cls" --marker-cluster "$cl" --records 50 --points secondary
    done
  done
  py analyses/per_type/drift_sensitivity_figure.py
  py analyses/per_type/composition_by_collection.py
  py analyses/per_type/painted_by_type.py                           # the painted class by type, 1940 collections
  # The two further sets of simulations (Tables S7 and S8).
  py analyses/per_type/drift_extension.py
  py analyses/per_type/drift_extension_joint.py
  py analyses/per_type/type_origination.py                          # Table S7
  py analyses/per_type/type_origination_joint.py
  py analyses/per_type/type_origination_calibration.py              # Table S8
  # The rim attributes: interior beveling, exterior notching and lip form
  # (Figure 7's findings file). Paste is recorded in the source and not analyzed.
  py scripts/build_mainfort_rim_counts.py
  py analyses/rims/rim_r1.py --points primary --kinds rim,decoration_current,decoration_pfg --records 50
  py analyses/rims/rim_r1.py --points secondary --kinds rim,decoration_current,decoration_pfg --records 50
  py analyses/rims/rim_r1.py --readings
fi

if [ "$LEVEL" -ge 3 ]; then
  echo "=== level 3: the Stan fits (every fit passes through R/fit_spatial.R) ==="
  run Rscript analyses/probe_latent_groups/refit_a0_basin.R           # A0 on the basin; Table S4
  run Rscript analyses/probe_latent_groups/a0_intercept_prior.R       # Figure 5
  # The validation records take their field scale and range from the basin fit above, so these
  # three steps must follow it and precede the validation fits.
  run Rscript analyses/probe_latent_groups/a0_pinned_nuisances.R      # medians of the basin fit
  py analyses/probe_latent_groups/a1_contracts.py                     # latent-unit validation records, data/stan/a1_*.json
  py analyses/probe_latent_groups/a3_contracts.py                     # barrier validation records, data/stan/a3_*.json
  run Rscript analyses/probe_latent_groups/validate_a1.R              # Figure 8A
  run Rscript analyses/probe_latent_groups/validate_a3.R              # Figure 8B
  run Rscript analyses/probe_latent_groups/a3_pinned_range.R 6        # Table S5
  run Rscript analyses/probe_latent_groups/a3_pinned_range.R 12       # Table S5
  # The same validation under the Normal(0, 2) intercept prior (the supplement's prior-sensitivity
  # entries and the Normal(0, 2) column of Table S5). Its records are simulated at the values of
  # the Normal(0, 2) basin fit, which a0_intercept_prior.R wrote above, and are kept in
  # data/stan/alphaN0-2/; its tables are the *_alphaN0-2.md files.
  run Rscript analyses/probe_latent_groups/a0_pinned_nuisances.R --alpha-sd=2
  py analyses/probe_latent_groups/a1_contracts.py --alpha-sd 2
  py analyses/probe_latent_groups/a3_contracts.py --alpha-sd 2
  run Rscript analyses/probe_latent_groups/validate_a1.R --alpha-sd=2
  run Rscript analyses/probe_latent_groups/validate_a3.R --alpha-sd=2
  run Rscript analyses/probe_latent_groups/a3_pinned_range.R 6 --alpha-sd=2
  run Rscript analyses/probe_latent_groups/a3_pinned_range.R 12 --alpha-sd=2
  py analyses/probe_latent_groups/a0_predictive_summaries.py
  py analyses/probe_latent_groups/a3_prior_comparison.py
  py analyses/probe_latent_groups/export_fit_diagnostics.py           # fit_diagnostics.md, from the fit records
fi

echo "=== level 1: the figures, from the committed findings files ==="
py analyses/93_phillips1970_phase_map.py           # Figure 1
py analyses/102_fig2_two_histories.py              # Figure 2
py analyses/103_fig3_model_schematic.py            # Figure 3
py analyses/104_fig4_phase_recovery_map.py         # Figure 4
py analyses/105_fig5_spatial_field.py              # Figure 5
py analyses/106_fig6_per_class_drift.py            # Figure 6
py analyses/107_fig7_rim_attributes.py             # Figure 7
py analyses/108_fig8_latent_unit_recovery.py       # Figure 8
py analyses/11_chronology_14c.py                   # Figure S1 (figures/figS1_chronology)
# Figure S2 (figures/figS2_ranksize) is drawn by analyses/make_figures.py,
# which also redraws the earlier version's house-style figures:
#   $PY analyses/make_figures.py

echo "=================================================================="
echo "Done. ${fail} script(s) reported a non-zero exit."
