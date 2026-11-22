# Corrections

Corrections to this repository after its first public release, newest first.
Each entry states what was wrong, what changed, and whether any number in the
manuscript is affected.

---

## 2026-11-21 — paste removed from the rim analysis; reproduction instructions corrected

**Affects:** `analyses/rims/rim_r1.py`, `analyses/rims/rim_r1_readings.py`,
`output/findings/rim_r1.md`, Figure 7, `README.md`, `run_all.sh`,
`data/stan/basin_composition.json`.
**Manuscript affected:** **Yes, by removal only.** The paper no longer reports
paste (Bell against Neeley's Ferry). The values for interior beveling, exterior
notching, lip form and the decorated classes are unchanged.

### What was wrong

1. The rim analysis reported paste beside the three rim attributes. The paper
   excludes plain sherds because paste and temper can reflect firing conditions
   and vessel function, and the Bell and Neeley's Ferry pastes are type-level
   distinctions of the 1951 typology. Reporting paste was inconsistent with
   that exclusion.
2. `README.md` said that reruns reproduce the committed findings files. That
   holds only under the NumPy version that produced each file. The simulator
   time scales and the per-class, regional-drift and rim tables were produced
   under NumPy 2.5.4, and the repository pins NumPy below 2.5.
3. `README.md` gave `pytest` as the test command, and neither
   `requirements.txt` nor `environment.yml` installs pytest.
4. `data/stan/basin_composition.json` was an older write of
   `analyses/49_write_stan_data.py`. It lacked the `prior_alpha_mean` and
   `prior_alpha_sd` fields the spatial field model's Stan program declares,
   and the provenance block.

### What changed

1. Paste is removed from the rim simulations, the readings and Figure 7. The
   findings file lists the paste counts and reports no statistic for them.
2. `README.md` states the version each seeded result needs and gives the
   commands for a NumPy 2.5.4 environment. Checked on 2026-11-21: under 2.5.4
   the time scales and the rim analysis are reproduced; the per-class and
   regional-drift tables were not rerun.
3. `README.md` installs the `dev` extras before running `pytest`.
4. `data/stan/basin_composition.json` is the current output of
   `analyses/49_write_stan_data.py`. The counts are identical to the earlier
   file. Distances involving three sites differ in the twelfth decimal place.
5. New: `analyses/per_type/painted_by_type.py` and
   `output/findings/painted_by_type.md`, the painted class by type in the 1940
   collections.

---

## 2026-11-01 — the repository now carries the restructured paper

**Affects:** the whole release.
**Manuscript affected:** the paper is a new manuscript, "What Are Phases?
Geographic Structure, Neutral Copying, and the Search for Social Boundaries in
the Lower Mississippi Valley", replacing "Are the phases real?". Its analyses
are the spatial field, latent-unit and barrier Stan models under
`stan/spatial/` with `R/fit_spatial.R`, the neutral copying model and its
diversity-matched comparisons under `analyses/probe_latent_groups/` and
`analyses/per_type/`, the rim analysis under `analyses/rims/`, and the phase
recovery and settlement analyses 76 to 83. The figures are drawn by
`analyses/93_phillips1970_phase_map.py` and `102_` to `108_` from the committed
findings files. The earlier version's scripts, figures and result reports
remain in the repository; `README.md` says which files the current paper
depends on. The two rim data files that were withheld are included from 2026-11-14 (`data/README.md`).

---

## 2026-09-21 to 2026-09-22 — the analysis matrix double-counted sherds; rebuilt from two sources

**Affects:** every count-based result. `data/processed/analysis_matrix.csv`
replaces the `pfg-cpl-mainfort` sheet of `data/raw/mainfort-pfg-cpl.xlsx` as
the decorated-class matrix.
**Manuscript affected:** **Yes.** Every quantitative result of the paper was
rerun from the new matrix and membership.

### What was wrong

The sheet every earlier result rested on summed Lipo's (2001) compilation with
a `Mainfort` sheet. The two are not independent: the `Mainfort` sheet's counts
are transcriptions of the same Lower Mississippi Survey files, so Phillips,
Ford and Griffin's sherds were counted twice. The `Mainfort` sheet's Parkin
Punctated column was also doubled relative to Mainfort's published Table 1
(`scripts/check_matrix_against_mainfort2003.py`;
`output/findings/matrix_integrity.md`: twice the published count at 34 of 38
shared sites, ratio 2.00 to 2.15 at all 38).

### What changed

`scripts/build_analysis_matrix.py` builds the matrix from two sources with
one source per assemblage and never a sum: the `cpl-pfg` tab of Lipo's
compilation where an assemblage appears there, otherwise `data/raw/PFGData.xlsx`.
Mainfort's (2003a) published table is kept as a separate comparison set
(`data/processed/mainfort2003_matrix.csv`). Each row's source and its earlier
and corrected totals are in `data/processed/analysis_matrix_provenance.csv`
(Parkin, for example: 865 Parkin Punctated from the definitive counts against
1,298 in Mainfort's table). At the same time the Parchman phase was excluded on
geography and a 75-decorated-sherd minimum adopted
(`output/findings/minimum_sample_size.md`), giving the 28-assemblage study
set with 14,101 decorated sherds, and five assemblage coordinates and five
settlement coordinates were corrected from PFG's legal descriptions
(`data/raw/coordinate_corrections.csv`,
`data/raw/settlement_coordinate_corrections.csv`). The datum of PFG UTMs read
as NAD83 was corrected to NAD27 where demonstrable, moving 13 of 28 members by
209 m with no change to membership or results
(`output/findings/datum_conversion.md`).

---

## 2026-08-31 — the Bayesian cultural F_ST was computed on 55 assemblages, not the 29-assemblage basin

**Affects:** `analyses/43_bayesian_fst.py`, `analyses/44_bayesian_fst_validation.py`,
and the file they write, `output/bayesian_fst.md`.
**Manuscript affected:** **No.** See "Scope" below.

### What was wrong

`basin_group_counts()` grouped the output of
`07_refined_empirical.prepare_inputs()` by `inp.cluster_of` over
`inp.have_coords_ids` and never restricted to the canonical drainage-basin
membership in `data/processed/basin_members_curated.txt`. It therefore fitted
the whole curated set with coordinates, **55 assemblages in 5 spatial
clusters**, while being named for the basin and writing a report headed
"Observed St. Francis basin". The 26 extra assemblages include `Walls`,
`Wall`, `Chuccalissa` and `Parchman`, which are the Mississippi-proximal sites
that `analyses/16_basin_membership.py` names in its own docstring as the reason
the earlier latitude cut was replaced by a hydrological rule.

### What the corrected numbers are

Measured by the new `analyses/47_basin_scope_check.py`, which fits both scopes
side by side:

| | as published (55 assemblages, k=5) | corrected (29 assemblages, k=3) |
|---|---|---|
| Balding-Nichols F_ST median | 0.0674 | **0.0734** |
| BN 95% credible interval | [0.0346, 0.1160] | [0.0260, 0.1583] |
| Gini-Simpson F_ST median | 0.0327 | **0.0179** |
| Gini-Simpson 95% interval | [0.0301, 0.0352] | [0.0157, 0.0201] |
| frequentist plug-in | 0.0330 | 0.0179 |

Both fits are numerically healthy (R-hat 1.001, minimum ESS above 4500, zero
divergent transitions), so this is a scope error and not a sampling problem.
The three spatial clusters in the corrected fit are what the silhouette
criterion selects on the basin's own coordinates, which is the number the
manuscript describes.

### Scope: no manuscript number changes

The manuscript's figures and reported statistics reach the data through
`make_figures._load_curated()`, which *does* apply the basin membership. The
defect was confined to analyses 43 and 44, whose output is a standalone
repository analysis and is not cited in the manuscript. It is corrected here
because the published code produced a mislabelled quantity, which anyone
reproducing this repository would have obtained.

### A note on the regenerated outputs and library versions

`output/bayesian_fst.md` and `figures/bayesian_fst.*` in this commit were
regenerated with the corrected code, but under **newer library versions than
this repository pins**: pymc 6.3.1 / pytensor 3.3.0 / arviz 1.3.0 / numpy 2.4.6,
against the pinned pymc 6.1.0 / pytensor 3.1.3 in `requirements-lock.txt`. The
README notes that PyMC output shifts across versions, and it does here: the
uncorrected fit gave a Balding-Nichols median of 0.0665 under the pinned stack
and 0.0674 under the newer one.

That drift, about 0.001, is an order of magnitude smaller than the correction it
sits inside (Gini-Simpson 0.0327 to 0.0179), so the shipped numbers are the
right ones to a far better tolerance than the error being fixed. They are
nonetheless not bit-identical to what the locked environment will produce.
Regenerating under the lock file is the authoritative route:

    pip install -r requirements-lock.txt && python analyses/43_bayesian_fst.py

Shipping the stale, demonstrably mislabelled output was judged worse than
shipping a corrected output with a stated version caveat.

### What changed in the code

- `basin_group_counts(inp, scope="basin")` now restricts to the canonical
  membership and re-selects the cluster count on the basin's own coordinates.
  `scope="region"` reproduces the previous 55-assemblage behaviour, retained
  deliberately as a legitimate measurement at a wider spatial scale rather than
  removed, so the scale dependence is visible instead of hidden. An unknown
  scope raises.
- `tests/inference/test_basin_scope.py` asserts that the fitted set equals
  `data/processed/basin_members_curated.txt` element by element, with the
  `scope="region"` probe as the point where the rival account ("the sets are
  effectively the same") predicts no difference and is wrong.
- `analyses/43_bayesian_fst.py` now reports the full diagnostic set (R-hat,
  bulk and tail ESS, divergences as count and percentage, treedepth saturation,
  E-BFMI) for **both** the primary fit and the Beta(1,3) prior-sensitivity
  refit. Previously the refit, which is the stated defence of the prior choice,
  reported an interval with no diagnostics behind it.

---

## 2026-08-31 — `h5py` was not declared

**Affects:** `pyproject.toml`, `requirements.txt`.
**Manuscript affected:** No.

`h5netcdf` dispatches to `h5py` for the actual HDF5 backend. Without it,
`idata.to_netcdf()` raises `ImportError` — but only *after* the markdown report
has been written, so a run appeared to succeed while silently saving no
posterior. Added with a version floor of 3.10.
