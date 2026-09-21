# MANIFEST

Which script produces which artifact, and the order they run in.

**Generated. Do not edit by hand.** Written by `scripts/build_manifest.py`, which parses the scripts rather than trusting a list. Regenerate after adding or renaming any analysis:

```
    .venv/bin/python scripts/build_manifest.py
    .venv/bin/python scripts/build_manifest.py --check   # in CI
```

83 scripts produce 165 tracked artifacts.

Writes are detected from the call that performs them (`write_text`, `to_csv`, `savefig`, `save_all`, `savez`, `writeLines`, `saveRDS`), with paths resolved through each script's own module-level constants. Reads are deliberately NOT counted: an earlier regex version credited `21_signal_recovery.py` with producing `output/closure_posterior.npz`, which it only consumes and `53` writes.

A figure is listed once, under whichever extension the script names. `figstyle.save_all` writes four siblings from one call (`.png`, `.svg`, `.pdf`, `.tiff`), and only the `.svg` is tracked in git; `figures/*.png` is ignored, so the PNGs a manuscript build consumes are derived artifacts that must be regenerated rather than cloned.

## Order

Topologically sorted by sibling import: a script may be run at any point after everything above it. The order is derived from the `importlib.import_module` edges between scripts, and verified to have no violations apart from the mutual pair noted below.

**What this order is not.** It has not been executed as a single top-to-bottom pass, and several entries take hours (the ABC-SMC fits, the 500-realisation drift simulations, the Stan GP fits). Import edges also do not capture every dependency: a script that reads an artifact another wrote, without importing it, is ordered correctly here only because the numbering already reflects that. Treat this as the dependency record rule 15 asks for, not as a proven build script.

> **Mutual imports.** `21_signal_recovery.py` and `53_closure_strength_posterior.py`; `35_basin_pullout.py` and `36_canonical_phase_map.py` import each other at module level. This works because each uses the other only inside a function body, so neither needs the other's attributes at import time, and either may be run first. It does mean the order below is arbitrary between them.

| # | script | produces |
|---|---|---|
| 1 | `analyses/00_setup/02_record_environment.R` | `output/findings/environment_r.md` |
| 2 | `analyses/00_setup/02_record_environment.py` | `output/findings/environment_python.md` |
| 3 | `analyses/02_spatial/01_recovery.R` | `output/findings/spatial_recovery.md` |
| 4 | `analyses/02_spatial/02_prior_repair.R` | `output/findings/spatial_prior_repair.md` |
| 5 | `analyses/02_spatial/03_geometry_repair.R` | `output/findings/gp_geometry_repair.md` |
| 6 | `analyses/02_spatial/04_posterior_predictive.R` | `output/findings/gp_posterior_predictive.md` |
| 7 | `analyses/02_spatial/05_basin_fit.R` | `output/findings/gp_basin_fit.md` |
| 8 | `analyses/04_validation_report.py` | `output/validation_report.md` |
| 9 | `analyses/05_empirical_application.py` | `figures/05_ordinal_trajectory.png`<br>`figures/05_rank_size_area.png`<br>`figures/05_spatial_clusters.png`<br>`output/empirical_findings.md` |
| 10 | `analyses/06_empirical_two_level.py` | `figures/06_ca_trajectory.png`<br>`figures/06_ca_vs_14c.png`<br>`figures/06_mound_ditch_map.png`<br>`figures/06_mound_rank_size.png`<br>`output/empirical_findings_v2.md` |
| 11 | `analyses/07_refined_empirical.py` | `figures/07_early_late_contrast.png`<br>`figures/07_idss_bridge_rank.png`<br>`figures/07_signature_trajectory.png` |
| 12 | `analyses/08_coupling_robustness.py` | `figures/08_coupling_robustness.png`<br>`output/coupling_robustness.md` |
| 13 | `analyses/09_parkin_basin_restricted.py` | `figures/09_basin_rank_size.png`<br>`figures/09_basin_signature_trajectory.png`<br>`output/parkin_basin_restricted.md` |
| 14 | `analyses/10_neiman_power_diagnostics.py` | `output/neiman_power_diagnostics.md` |
| 15 | `analyses/11_chronology_14c.py` | `figures/figS6_chronology.svg`<br>`output/chronology_14c.md` |
| 16 | `analyses/12_sensitivity_grid.py` | `output/sensitivity_grid.md` |
| 17 | `analyses/13_neiman_distance_and_fit.py` | `figures/figS1_neiman.svg`<br>`output/neiman_distance_and_fit.md` |
| 18 | `analyses/15_continuum_test.py` | `output/continuum_test.md` |
| 19 | `analyses/23_phases_vs_spatial_drift.py` | `output/phases_vs_spatial_drift.md` |
| 20 | `analyses/26_cmv_phase_groupness.py` | `figures/figX_cmv_mds.png`<br>`output/cmv_phase_groupness.md` |
| 21 | `analyses/27_cmv_lmv_repertoire.py` | `figures/figX_cmv_lmv_repertoire.png`<br>`output/cmv_lmv_repertoire.md` |
| 22 | `analyses/30_regional_map.py` | `figures/figS7_regional.png` |
| 23 | `analyses/42_figS6_dynamic.py` | `figures/figS5_dynamic.svg` |
| 24 | `analyses/51_bayesian_rank_correlations.py` | `output/findings/bayesian_rank_correlations.md` |
| 25 | `analyses/52_bayesian_fit_increments.py` | `output/findings/bayesian_increment_test.md` |
| 26 | `analyses/58_dm_numerical_floor.py` | `output/findings/dm_numerical_floor.md` |
| 27 | `analyses/63_cmv_site_coords.py` | `data/processed/williams1954_cmv_coords.tsv` |
| 28 | `analyses/68_mound_lidar_heights.py` | `output/findings/mound_lidar_heights.md` |
| 29 | `analyses/69_mound_outline_measures.py` | `output/findings/mound_outline_measures.csv`<br>`output/findings/mound_outline_measures.md` |
| 30 | `analyses/70_historic_quad_tiles.py` | `output/lidar/manifest.json` |
| 31 | `analyses/make_figures.py` | `figures/fig2_validation.svg`<br>`figures/fig4_ca_ordination.svg`<br>`figures/fig5_idss_network.svg`<br>`figures/fig6_empirical_trajectory.svg`<br>`figures/fig7_idss_structure.svg`<br>`figures/fig8_ranksize.svg` |
| 32 | `analyses/make_map.py` | `figures/fig1_studyarea.svg` |
| 33 | `analyses/14_drainage_basin.py` | `output/drainage_basin.md` |
| 34 | `analyses/17_basin_results.py` | `output/basin_results.md` |
| 35 | `analyses/24_phases_drift_robustness.py` | `output/phases_drift_robustness.md` |
| 36 | `analyses/25_drift_vs_groups_demo.py` | `figures/figS2_drift_vs_groups.png`<br>`output/drift_vs_groups_demo.md` |
| 37 | `analyses/40_hierarchical_convergence.py` | `figures/hierarchical_convergence_slopes.svg`<br>`output/hierarchical_convergence.md` |
| 38 | `analyses/43_bayesian_fst.py` | `figures/bayesian_fst.svg`<br>`output/bayesian_fst.md` |
| 39 | `analyses/45_alvey_14c_robustness.py` | `output/alvey_14c_robustness.md` |
| 40 | `analyses/48_length_scale_recovery.py` | `figures/fig_length_scale_recovery.svg`<br>`output/findings/length_scale_recovery.md` |
| 41 | `analyses/48_revision_validation.py` | `output/revision_idealized_validation.json` |
| 42 | `analyses/49_write_stan_data.py` | `data/stan/basin_composition.json`<br>`data/stan/sim_truth.json` |
| 43 | `analyses/59_partition_sensitivity.py` | `figures/fig_partition_sensitivity.svg`<br>`output/findings/partition_sensitivity.md` |
| 44 | `analyses/61_river_network_geometry.py` | `output/findings/river_network_geometry.md` |
| 45 | `analyses/18_kandler_shennan_neutrality.py` | `output/kandler_shennan_neutrality.md` |
| 46 | `analyses/19_abc_transmission.py` | `output/abc_posterior.npz`<br>`output/abc_transmission.md` |
| 47 | `analyses/20_tempo_mode_ews.py` | `output/tempo_akaike.npz`<br>`output/tempo_mode_ews.md` |
| 48 | `analyses/28_macro_boundary.py` | `figures/figX_macro_transect.png`<br>`output/macro_boundary.md` |
| 49 | `analyses/31_within_region_structure.py` | `figures/figS8_within_region.png` |
| 50 | `analyses/33_time_aware_emergence.py` | `figures/figS3_emergent_phases.png`<br>`output/time_aware_emergence.md`<br>`output/time_aware_parkin_prob.csv`<br>`output/time_aware_runs.csv` |
| 51 | `analyses/41_hierarchical_convergence_validation.py` | `figures/hierarchical_convergence_calibration.svg`<br>`output/hierarchical_convergence_validation.md` |
| 52 | `analyses/44_bayesian_fst_validation.py` | `figures/bayesian_fst_coverage.svg`<br>`output/bayesian_fst_validation.md` |
| 53 | `analyses/46_radiocarbon_table.py` | `output/radiocarbon_dates_used.md` |
| 54 | `analyses/47_basin_scope_check.py` | `output/findings/basin_scope_check.md` |
| 55 | `analyses/50_perbin_bayesian_fst.py` | `output/findings/perbin_bayesian_fst.md` |
| 56 | `analyses/54_tempo_mode_posterior.py` | `output/findings/tempo_mode_posterior.md`<br>`output/tempo_posterior.npz` |
| 57 | `analyses/56_neutrality_ppc.py` | `figures/fig_neutrality_ppc.svg`<br>`output/findings/neutrality_ppc.md` |
| 58 | `analyses/57_fst_prior_predictive.py` | `figures/fig_fst_prior_predictive.svg`<br>`output/findings/fst_prior_justification.md` |
| 59 | `analyses/60_partition_ensemble.py` | `figures/fig_partition_ensemble.svg`<br>`output/findings/partition_ensemble.md` |
| 60 | `analyses/34_emergence_robustness.py` | `figures/figS4_emergence_robustness.png`<br>`output/emergence_robustness.csv`<br>`output/emergence_robustness.md` |
| 61 | `analyses/38_abc_smc_transmission.py` | `output/abc_smc_posterior.npz`<br>`output/abc_smc_transmission.md` |
| 62 | `analyses/39_abc_smc_validation.py` | `figures/abc_smc_crosscheck.svg`<br>`figures/abc_smc_sbc.svg`<br>`output/abc_smc_sbc_ranks.npz`<br>`output/abc_smc_validation.md` |
| 63 | `analyses/55_abc_stability.py` | `output/abc_pooled_posterior.npz`<br>`output/findings/abc_stability.md` |
| 64 | `analyses/16_basin_membership.py` | `data/processed/basin_members_broad.txt`<br>`data/processed/basin_members_curated.txt`<br>`data/processed/basin_members_watershed.txt` |
| 65 | `analyses/21_signal_recovery.py` | `figures/fig5_recovery.svg`<br>`figures/fig6_empirical_trajectory.svg`<br>`output/revision_recovery.json`<br>`output/signal_recovery.md` |
| 66 | `analyses/22_generator_diagnostic.py` | `output/generator_diagnostic.md` |
| 67 | `analyses/53_closure_strength_posterior.py` | `figures/fig_closure_strength_posterior.svg`<br>`output/closure_posterior.npz`<br>`output/findings/closure_strength_posterior.md` |
| 68 | `analyses/62_gini_simpson_bias.py` | `output/findings/gini_simpson_bias.md` |
| 69 | `analyses/29_concept_figure.py` | `figures/fig3_concept.png` |
| 70 | `analyses/35_basin_pullout.py` | `figures/fig9_parkin_pullout.png`<br>`output/basin_pullout.md`<br>`output/basin_pullout_prob.csv`<br>`output/basin_pullout_runs.csv` |
| 71 | `analyses/36_canonical_phase_map.py` | `figures/fig1_phases.png` |
| 72 | `analyses/37_lmv_drift_groups.py` | `figures/fig8_lmv_drift_groups.png`<br>`output/lmv_drift_groups.md`<br>`output/lmv_drift_groups_example.csv`<br>`output/lmv_drift_groups_runs.csv` |
| 73 | `analyses/47_revision_analysis.py` | `output/revision_2026_09/baseline.csv`<br>`output/revision_2026_09/boundary_grid.csv`<br>`output/revision_2026_09/calibrated_rates.csv`<br>`output/revision_2026_09/calibration.csv`<br>`output/revision_2026_09/co_membership.csv`<br>`output/revision_2026_09/initialization_sensitivity.csv`<br>`output/revision_2026_09/leave_one_out.csv`<br>`output/revision_2026_09/manifest.json`<br>`output/revision_2026_09/report.md`<br>`output/revision_2026_09/summary.json` |
| 74 | `analyses/50_revision_figures.py` | `figures/fig8_lmv_drift_groups.svg`<br>`figures/fig9_parkin_pullout.svg`<br>`figures/figS2_drift_vs_groups.svg`<br>`figures/figS3_emergent_phases.svg`<br>`figures/figS4_emergence_robustness.svg`<br>`figures/figS8_within_region.svg` |
| 75 | `analyses/64_unequal_populations.py` | `output/findings/unequal_populations.md`<br>`output/unequal_populations.json` |
| 76 | `analyses/65_other_departures.py` | `output/findings/other_departures.md`<br>`output/other_departures.json` |
| 77 | `analyses/66_cmv_class_composition.py` | `output/findings/cmv_class_composition.md` |
| 78 | `analyses/71_scale_sweep.py` | `output/findings/scale_sweep.csv`<br>`output/findings/scale_sweep.md` |
| 79 | `analyses/72_excess_locality.py` | `output/findings/excess_locality.csv`<br>`output/findings/excess_locality.md` |
| 80 | `analyses/74_phase_partition_test.py` | `figures/fig12_phase_partition.svg`<br>`output/findings/phase_partition_ensembles.csv`<br>`output/findings/phase_partition_test.md` |
| 81 | `analyses/73_connectivity_mixing.py` | `output/findings/connectivity_mixing.md` |
| 82 | `analyses/75_groupness_surface.py` | `figures/fig13_groupness_surface.svg`<br>`output/findings/groupness_surface.md` |
| 83 | `analyses/76_phase_recovery.py` | `figures/fig14_phase_recovery.svg`<br>`output/findings/phase_recovery.csv`<br>`output/findings/phase_recovery.md` |

## Artifacts, by path

| artifact | produced by |
|---|---|
| `data/processed/basin_members_broad.txt` | `16_basin_membership.py` |
| `data/processed/basin_members_curated.txt` | `16_basin_membership.py` |
| `data/processed/basin_members_watershed.txt` | `16_basin_membership.py` |
| `data/processed/williams1954_cmv_coords.tsv` | `63_cmv_site_coords.py` |
| `data/stan/basin_composition.json` | `49_write_stan_data.py` |
| `data/stan/sim_truth.json` | `49_write_stan_data.py` |
| `figures/05_ordinal_trajectory.png` | `05_empirical_application.py` |
| `figures/05_rank_size_area.png` | `05_empirical_application.py` |
| `figures/05_spatial_clusters.png` | `05_empirical_application.py` |
| `figures/06_ca_trajectory.png` | `06_empirical_two_level.py` |
| `figures/06_ca_vs_14c.png` | `06_empirical_two_level.py` |
| `figures/06_mound_ditch_map.png` | `06_empirical_two_level.py` |
| `figures/06_mound_rank_size.png` | `06_empirical_two_level.py` |
| `figures/07_early_late_contrast.png` | `07_refined_empirical.py` |
| `figures/07_idss_bridge_rank.png` | `07_refined_empirical.py` |
| `figures/07_signature_trajectory.png` | `07_refined_empirical.py` |
| `figures/08_coupling_robustness.png` | `08_coupling_robustness.py` |
| `figures/09_basin_rank_size.png` | `09_parkin_basin_restricted.py` |
| `figures/09_basin_signature_trajectory.png` | `09_parkin_basin_restricted.py` |
| `figures/abc_smc_crosscheck.svg` | `39_abc_smc_validation.py` |
| `figures/abc_smc_sbc.svg` | `39_abc_smc_validation.py` |
| `figures/bayesian_fst.svg` | `43_bayesian_fst.py` |
| `figures/bayesian_fst_coverage.svg` | `44_bayesian_fst_validation.py` |
| `figures/fig12_phase_partition.svg` | `74_phase_partition_test.py` |
| `figures/fig13_groupness_surface.svg` | `75_groupness_surface.py` |
| `figures/fig14_phase_recovery.svg` | `76_phase_recovery.py` |
| `figures/fig1_phases.png` | `36_canonical_phase_map.py` |
| `figures/fig1_studyarea.svg` | `make_map.py` |
| `figures/fig2_validation.svg` | `make_figures.py` |
| `figures/fig3_concept.png` | `29_concept_figure.py` |
| `figures/fig4_ca_ordination.svg` | `make_figures.py` |
| `figures/fig5_idss_network.svg` | `make_figures.py` |
| `figures/fig5_recovery.svg` | `21_signal_recovery.py` |
| `figures/fig6_empirical_trajectory.svg` | `21_signal_recovery.py`, `make_figures.py` |
| `figures/fig7_idss_structure.svg` | `make_figures.py` |
| `figures/fig8_lmv_drift_groups.png` | `37_lmv_drift_groups.py` |
| `figures/fig8_lmv_drift_groups.svg` | `50_revision_figures.py` |
| `figures/fig8_ranksize.svg` | `make_figures.py` |
| `figures/fig9_parkin_pullout.png` | `35_basin_pullout.py` |
| `figures/fig9_parkin_pullout.svg` | `50_revision_figures.py` |
| `figures/figS1_neiman.svg` | `13_neiman_distance_and_fit.py` |
| `figures/figS2_drift_vs_groups.png` | `25_drift_vs_groups_demo.py` |
| `figures/figS2_drift_vs_groups.svg` | `50_revision_figures.py` |
| `figures/figS3_emergent_phases.png` | `33_time_aware_emergence.py` |
| `figures/figS3_emergent_phases.svg` | `50_revision_figures.py` |
| `figures/figS4_emergence_robustness.png` | `34_emergence_robustness.py` |
| `figures/figS4_emergence_robustness.svg` | `50_revision_figures.py` |
| `figures/figS5_dynamic.svg` | `42_figS6_dynamic.py` |
| `figures/figS6_chronology.svg` | `11_chronology_14c.py` |
| `figures/figS7_regional.png` | `30_regional_map.py` |
| `figures/figS8_within_region.png` | `31_within_region_structure.py` |
| `figures/figS8_within_region.svg` | `50_revision_figures.py` |
| `figures/figX_cmv_lmv_repertoire.png` | `27_cmv_lmv_repertoire.py` |
| `figures/figX_cmv_mds.png` | `26_cmv_phase_groupness.py` |
| `figures/figX_macro_transect.png` | `28_macro_boundary.py` |
| `figures/fig_closure_strength_posterior.svg` | `53_closure_strength_posterior.py` |
| `figures/fig_fst_prior_predictive.svg` | `57_fst_prior_predictive.py` |
| `figures/fig_length_scale_recovery.svg` | `48_length_scale_recovery.py` |
| `figures/fig_neutrality_ppc.svg` | `56_neutrality_ppc.py` |
| `figures/fig_partition_ensemble.svg` | `60_partition_ensemble.py` |
| `figures/fig_partition_sensitivity.svg` | `59_partition_sensitivity.py` |
| `figures/hierarchical_convergence_calibration.svg` | `41_hierarchical_convergence_validation.py` |
| `figures/hierarchical_convergence_slopes.svg` | `40_hierarchical_convergence.py` |
| `output/abc_pooled_posterior.npz` | `55_abc_stability.py` |
| `output/abc_posterior.npz` | `19_abc_transmission.py` |
| `output/abc_smc_posterior.npz` | `38_abc_smc_transmission.py` |
| `output/abc_smc_sbc_ranks.npz` | `39_abc_smc_validation.py` |
| `output/abc_smc_transmission.md` | `38_abc_smc_transmission.py` |
| `output/abc_smc_validation.md` | `39_abc_smc_validation.py` |
| `output/abc_transmission.md` | `19_abc_transmission.py` |
| `output/alvey_14c_robustness.md` | `45_alvey_14c_robustness.py` |
| `output/basin_pullout.md` | `35_basin_pullout.py` |
| `output/basin_pullout_prob.csv` | `35_basin_pullout.py` |
| `output/basin_pullout_runs.csv` | `35_basin_pullout.py` |
| `output/basin_results.md` | `17_basin_results.py` |
| `output/bayesian_fst.md` | `43_bayesian_fst.py` |
| `output/bayesian_fst_validation.md` | `44_bayesian_fst_validation.py` |
| `output/chronology_14c.md` | `11_chronology_14c.py` |
| `output/closure_posterior.npz` | `53_closure_strength_posterior.py` |
| `output/cmv_lmv_repertoire.md` | `27_cmv_lmv_repertoire.py` |
| `output/cmv_phase_groupness.md` | `26_cmv_phase_groupness.py` |
| `output/continuum_test.md` | `15_continuum_test.py` |
| `output/coupling_robustness.md` | `08_coupling_robustness.py` |
| `output/drainage_basin.md` | `14_drainage_basin.py` |
| `output/drift_vs_groups_demo.md` | `25_drift_vs_groups_demo.py` |
| `output/emergence_robustness.csv` | `34_emergence_robustness.py` |
| `output/emergence_robustness.md` | `34_emergence_robustness.py` |
| `output/empirical_findings.md` | `05_empirical_application.py` |
| `output/empirical_findings_v2.md` | `06_empirical_two_level.py` |
| `output/findings/abc_stability.md` | `55_abc_stability.py` |
| `output/findings/basin_scope_check.md` | `47_basin_scope_check.py` |
| `output/findings/bayesian_increment_test.md` | `52_bayesian_fit_increments.py` |
| `output/findings/bayesian_rank_correlations.md` | `51_bayesian_rank_correlations.py` |
| `output/findings/closure_strength_posterior.md` | `53_closure_strength_posterior.py` |
| `output/findings/cmv_class_composition.md` | `66_cmv_class_composition.py` |
| `output/findings/connectivity_mixing.md` | `73_connectivity_mixing.py` |
| `output/findings/dm_numerical_floor.md` | `58_dm_numerical_floor.py` |
| `output/findings/environment_python.md` | `00_setup/02_record_environment.py` |
| `output/findings/environment_r.md` | `00_setup/02_record_environment.R` |
| `output/findings/excess_locality.csv` | `72_excess_locality.py` |
| `output/findings/excess_locality.md` | `72_excess_locality.py` |
| `output/findings/fst_prior_justification.md` | `57_fst_prior_predictive.py` |
| `output/findings/gini_simpson_bias.md` | `62_gini_simpson_bias.py` |
| `output/findings/gp_basin_fit.md` | `02_spatial/05_basin_fit.R` |
| `output/findings/gp_geometry_repair.md` | `02_spatial/03_geometry_repair.R` |
| `output/findings/gp_posterior_predictive.md` | `02_spatial/04_posterior_predictive.R` |
| `output/findings/groupness_surface.md` | `75_groupness_surface.py` |
| `output/findings/length_scale_recovery.md` | `48_length_scale_recovery.py` |
| `output/findings/mound_lidar_heights.md` | `68_mound_lidar_heights.py` |
| `output/findings/mound_outline_measures.csv` | `69_mound_outline_measures.py` |
| `output/findings/mound_outline_measures.md` | `69_mound_outline_measures.py` |
| `output/findings/neutrality_ppc.md` | `56_neutrality_ppc.py` |
| `output/findings/other_departures.md` | `65_other_departures.py` |
| `output/findings/partition_ensemble.md` | `60_partition_ensemble.py` |
| `output/findings/partition_sensitivity.md` | `59_partition_sensitivity.py` |
| `output/findings/perbin_bayesian_fst.md` | `50_perbin_bayesian_fst.py` |
| `output/findings/phase_partition_ensembles.csv` | `74_phase_partition_test.py` |
| `output/findings/phase_partition_test.md` | `74_phase_partition_test.py` |
| `output/findings/phase_recovery.csv` | `76_phase_recovery.py` |
| `output/findings/phase_recovery.md` | `76_phase_recovery.py` |
| `output/findings/river_network_geometry.md` | `61_river_network_geometry.py` |
| `output/findings/scale_sweep.csv` | `71_scale_sweep.py` |
| `output/findings/scale_sweep.md` | `71_scale_sweep.py` |
| `output/findings/spatial_prior_repair.md` | `02_spatial/02_prior_repair.R` |
| `output/findings/spatial_recovery.md` | `02_spatial/01_recovery.R` |
| `output/findings/tempo_mode_posterior.md` | `54_tempo_mode_posterior.py` |
| `output/findings/unequal_populations.md` | `64_unequal_populations.py` |
| `output/generator_diagnostic.md` | `22_generator_diagnostic.py` |
| `output/hierarchical_convergence.md` | `40_hierarchical_convergence.py` |
| `output/hierarchical_convergence_validation.md` | `41_hierarchical_convergence_validation.py` |
| `output/kandler_shennan_neutrality.md` | `18_kandler_shennan_neutrality.py` |
| `output/lidar/manifest.json` | `70_historic_quad_tiles.py` |
| `output/lmv_drift_groups.md` | `37_lmv_drift_groups.py` |
| `output/lmv_drift_groups_example.csv` | `37_lmv_drift_groups.py` |
| `output/lmv_drift_groups_runs.csv` | `37_lmv_drift_groups.py` |
| `output/macro_boundary.md` | `28_macro_boundary.py` |
| `output/neiman_distance_and_fit.md` | `13_neiman_distance_and_fit.py` |
| `output/neiman_power_diagnostics.md` | `10_neiman_power_diagnostics.py` |
| `output/other_departures.json` | `65_other_departures.py` |
| `output/parkin_basin_restricted.md` | `09_parkin_basin_restricted.py` |
| `output/phases_drift_robustness.md` | `24_phases_drift_robustness.py` |
| `output/phases_vs_spatial_drift.md` | `23_phases_vs_spatial_drift.py` |
| `output/radiocarbon_dates_used.md` | `46_radiocarbon_table.py` |
| `output/revision_2026_09/baseline.csv` | `47_revision_analysis.py` |
| `output/revision_2026_09/boundary_grid.csv` | `47_revision_analysis.py` |
| `output/revision_2026_09/calibrated_rates.csv` | `47_revision_analysis.py` |
| `output/revision_2026_09/calibration.csv` | `47_revision_analysis.py` |
| `output/revision_2026_09/co_membership.csv` | `47_revision_analysis.py` |
| `output/revision_2026_09/initialization_sensitivity.csv` | `47_revision_analysis.py` |
| `output/revision_2026_09/leave_one_out.csv` | `47_revision_analysis.py` |
| `output/revision_2026_09/manifest.json` | `47_revision_analysis.py` |
| `output/revision_2026_09/report.md` | `47_revision_analysis.py` |
| `output/revision_2026_09/summary.json` | `47_revision_analysis.py` |
| `output/revision_idealized_validation.json` | `48_revision_validation.py` |
| `output/revision_recovery.json` | `21_signal_recovery.py` |
| `output/sensitivity_grid.md` | `12_sensitivity_grid.py` |
| `output/signal_recovery.md` | `21_signal_recovery.py` |
| `output/tempo_akaike.npz` | `20_tempo_mode_ews.py` |
| `output/tempo_mode_ews.md` | `20_tempo_mode_ews.py` |
| `output/tempo_posterior.npz` | `54_tempo_mode_posterior.py` |
| `output/time_aware_emergence.md` | `33_time_aware_emergence.py` |
| `output/time_aware_parkin_prob.csv` | `33_time_aware_emergence.py` |
| `output/time_aware_runs.csv` | `33_time_aware_emergence.py` |
| `output/unequal_populations.json` | `64_unequal_populations.py` |
| `output/validation_report.md` | `04_validation_report.py` |

## Written by more than one script

Not necessarily wrong (a figure may be produced by the analysis that owns it and by a figure driver), but each is a place where running the pipeline in a different order gives a different file, so the owner should be the later entry in the order above.

- `figures/fig6_empirical_trajectory.svg`: `21_signal_recovery.py`, `make_figures.py`

## Not covered here

- `analyses/figstyle.py`, `make_map.py` are helpers imported by other scripts rather than entry points.
- `analyses/00_setup/*.R` set up and record the R toolchain; they produce no analysis artifact.
- Scripts that only print are omitted, since they leave nothing to regenerate.
- Inputs under `data/raw/` and `data/Shapefiles/` are tracked sources, not derived artifacts; see `data/DATA_INVENTORY.md`.
