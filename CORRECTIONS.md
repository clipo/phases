# Corrections

Corrections to this repository after its first public release, newest first.
Each entry states what was wrong, what changed, and whether any number in the
manuscript is affected.

---

## 2026-11-23 — the validation fits under the Normal(0, 2) intercept prior rerun

**Affects:** `run_all.sh`, `README.md`, `R/intercept_arm.R` (new),
`analyses/probe_latent_groups/` (`a0_pinned_nuisances.R`, `a1_contracts.py`,
`a3_contracts.py`, `validate_a1.R`, `validate_a3.R`, `a3_pinned_range.R`,
`a3_prior_comparison.py`, `export_fit_diagnostics.py`),
`tests/probe/test_intercept_arm.py` (new), `output/findings/` (the four
validation tables named `*_alphaN0-2.md`, `a3_prior_comparison.md`,
`fit_diagnostics.md` and `.csv`), `docs/manuscript/SUPPLEMENTAL_TEXT.*`,
`MANIFEST.md`.
**Manuscript affected:** **Yes, in values reported in the supplement. No
conclusion changed. The main text is unchanged.**

### What was wrong

The entry of 2026-11-22 left the validation tables under the Normal(0, 2)
class-intercept prior as the fits of October 2026. No committed script option
could rerun them, so the Normal(0, 2) column of Table S5 and the supplement's
latent-unit entry under that prior could not be reproduced from this
repository.

### What changed

- The record writers and the validation scripts take the argument
  `--alpha-sd=2` (`--alpha-sd 2` for the two Python scripts). With it they
  simulate the records at the field scale and range of the spatial field
  model's basin fit under Normal(0, 2), keep them in `data/stan/alphaN0-2/`,
  fit them under that prior and write the `*_alphaN0-2` tables. Each
  validation script stops if a record's stated prior is not the prior of the
  run. `./run_all.sh --stan` now runs both priors.
- The four validation tables under Normal(0, 2) were fitted on 2026-11-23
  (R 4.6.2, rstan 2.32.7, StanHeaders 2.39.1.1) and replace the October
  tables. Their records are simulated at a range of 202.6 km, where the
  October tables print 204 km.
- Latent-unit model under Normal(0, 2): 13 of the 14 recovery fits fail the
  diagnostic gate (12 in October), all on R-hat, and so does the fit to the
  record with no units. On that record the posterior median of the
  between-unit scale is 0.44 (0.04 to 1.20), where October gave 0.52 (0.05 to
  1.39). The medians on the five records at a true value of 0.5 are 0.27 to
  0.48, as before.
- Barrier model under Normal(0, 2): all 39 fits pass the gate at each of the
  three range settings, as before. Every Normal(0, 2) cell of Table S5 moved.
  The largest change is the upper end of the free-range row at a true value
  of 0.25, 0.84 where October gave 1.14. The two statements the table
  supports hold on the new values: with the range pinned, every range of
  medians lies above the prior median of 1.35, and the ranges of medians at
  adjacent true values overlap in every setting.
- `fit_diagnostics.md` now lists 283 fits: 145 under Normal(0, 10), 137 under
  Normal(0, 2) and one under the prior centered on the pooled log-odds. The
  October file remains as `fit_diagnostics_2026-10.md`.
- Supplement: the latent-unit entry under Normal(0, 2), the recovery-fits
  entry of the sampler diagnostics, and the Normal(0, 2) column of Table S5.

---

## 2026-11-22 — the Stan validation fits and four per-class tables rerun in pipeline order

**Affects:** `run_all.sh`, `README.md`,
`analyses/probe_latent_groups/` (`refit_a0_basin.R`, `validate_a1.R`,
`validate_a3.R`, `a3_pinned_range.R`, `a3_prior_comparison.py`,
`a1_prior_check.R`), `analyses/108_fig8_latent_unit_recovery.py`,
`output/findings/` (the `a0_`, `a1_`, `a3_` and `fit_diagnostics` files and
four `per_type_generators_*_secondary*` tables), `figures/`,
`docs/manuscript/`.
**Manuscript affected:** **Yes, in reported values. No conclusion changed.**

### What was wrong

1. `./run_all.sh --stan` could not run from a fresh clone. It omitted
   `a0_pinned_nuisances.R`, `a1_contracts.py` and `a3_contracts.py`, which
   write the simulated validation records, so the first validation fit stopped
   for lack of its data. `a3_pinned_range.R`, given no records, wrote an empty
   table over the committed one and exited normally.
2. The committed validation tables were fitted to records simulated with the
   field scale and range of the spatial field model's basin fit under the
   earlier Normal(0, 2) intercept prior (range 204 km). Run in order, the
   pipeline simulates them from the basin fit under the paper's Normal(0, 10)
   prior (range 99 km). The latent-unit validation of Figure 8A had been
   fitted under Normal(0, 2), with three fits repeated under Normal(0, 10).
3. `./run_all.sh --analyses` ran the planted-marker tables with three settings
   at both operating points. The committed tables record fifteen at the
   primary point and three at the secondary (corrected in pull request 12).
4. `refit_a0_basin.R` embedded a git diff against a commit this repository
   does not contain, which printed an error into its findings file.

### What changed

- Every Stan fit of level 3 was rerun on 2026-11-22 in the corrected order
  (R 4.6.2, rstan 2.32.7, StanHeaders 2.39.1.1), and the findings files,
  Figures 5 and 8 and the diagnostic figures are the results of that run. All
  14 latent-unit recovery fits now use Normal(0, 10). The October table under
  Normal(0, 2) is kept as `a1_validation_alphaN0-2.md`, and the October
  per-fit diagnostics as `fit_diagnostics_2026-10.md`.
- Conclusions that hold: all 39 free-range barrier fits pass the diagnostic
  gate and the barrier model does not separate no barrier from a barrier of
  1.0. Twelve of the 14 latent-unit recovery fits fail the gate. The spatial
  share is 0.93 under the main prior.
- Values that changed in the paper: the spatial-share intervals by 0.01. In
  the supplement, the prior-predictive and sensitivity values, Table S4, the
  recovery paragraphs, the diagnostics entries and the Normal(0, 10) column of
  Table S5. Barrier medians with the range free are 0.45 to 0.78 at no barrier
  and 0.82 to 2.31 at a barrier of 1.0. These ranges of medians do not
  overlap, so the supplement now states the overlap on the intervals: three of
  the six medians at 1.0 lie inside the 90 percent interval of every
  no-barrier record. With the range pinned at 12 km, 31 of 39 fits pass (37
  before).
- One spatial-field recovery fit (share 0.9, range 5 km) ended its third
  sampler setting with 12 divergent transitions in 24,000 draws under the
  current R and StanHeaders versions. A fourth setting was added and passes.
- Four secondary-point per-class tables differ from the earlier committed
  files in five band endpoints (four in the last printed digit, one lower
  endpoint 0.0050 against 0.0033). The difference persists with the original
  code, identical inputs and two sets of scipy and pandas versions. Its cause
  was not found. The tables are regenerated to the values a rerun under NumPy
  2.5.4 gives. The paper quotes none of the five values. The ten
  primary-point tables and the regional-drift summary reproduce exactly.
- The three validation scripts now stop with an error when their records are
  absent.

---

## 2026-11-22 — the painted class and the count of assemblages with later collections were misdescribed

**Affects:** `analyses/per_type/painted_by_type.py`,
`analyses/per_type/composition_by_collection.py`,
`output/findings/painted_by_type.md`,
`output/findings/composition_by_collection.md`,
`analyses/rims/rim_r1_readings.py`, `output/findings/rim_r1.md`,
`docs/manuscript/`.
**Manuscript affected:** **Yes, in wording and two counts. No reported
statistic changed, and the analysis matrix is unchanged.**

### What was wrong

1. The paper and `painted_by_type.md` said the painted class pools three types
   (Old Town Red, Nodena Red and White, Avenue Polychrome). Lipo (2001: Table
   7.4, note a) defines his painted category as five types, adding Carson
   Red-and-Buff and Hollywood White, and the matrix rows taken from his table
   follow it. Rows built from the survey spreadsheet by name use three. In the
   1940 collections of the 28 assemblages the two further types total 36
   sherds: 24 at Hollywood, counted in the painted class, and 12 at Williamson,
   not counted.
2. The paper and `composition_by_collection.md` said nine assemblages combine
   a 1940 collection with one Lipo made in 1996 and 1997. Six do. The three
   small differences (Hollywood +25, Barton Ranch +9, Neeleys Ferry +4) equal,
   in all but one sherd, 1940 counts of types his table pools: Carson
   Red-on-Buff and Hollywood White Slipped at Hollywood, and Manly Punctated at
   the other two. That Manly Punctated was pooled with Parkin Punctated is
   inferred from the equal counts.
3. The rim findings file and the paper said the rims were recorded by one
   observer. The source does not state that.
4. A sentence cited House (1993) for the Hollywood radiocarbon ranges, which
   that paper does not report, and a sentence on time-averaging stated its
   dependence on variant lifetime the wrong way round.

### What changed

Both scripts now compute and report the account in items 1 and 2, and a test
reads the expected counts from the spreadsheet directly. The manuscript states
the class definition, the six assemblages, and the size of the inconsistency.
The other sentences were corrected. The paper now spells site 12-O-5 Cranor
Place, following Lipo (2001); the data and findings files keep the survey's
spelling in the row name `Cramor_Place`.

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
