# Are the phases real?

Code and data for **"Are the phases real? Distinguishing bounded interaction
groups from spatially structured drift in lower Mississippi Valley decorated
ceramics"** by Carl P. Lipo, Robert J. DiNapoli, and Mark E. Madsen.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/clipo/phases/main)

Everything in the paper is reproducible from this repository: install an
environment, run one script, and the figures and result tables regenerate. Launch
a zero-install environment in the browser with Binder (badge above), or build
locally with conda, the conda-lock lockfile, or Docker (see [Installation](#installation)).

## The question

The *phase* is the foundational unit of Mississippian culture history. It rests
on a silent premise: that assemblages more similar to one another than to their
neighbors record a bounded social group whose members interacted preferentially
among themselves. That premise is rarely tested, and it has a strong
alternative. Under neutral cultural transmission, potters copy one another's
decorated styles with no one favoring particular designs; style frequencies
drift, and when potters interact more with near neighbors than distant ones, the
drift concentrates similar styles nearby. The result is the same
within-greater-than-between similarity, produced by geography rather than by a
social boundary.

This repository tests, for the late pre-contact lower Mississippi Valley,
whether decorated-ceramic similarity records bounded interaction groups or
spatially structured drift. It compares two generative processes (neutral
spatial drift versus bounded groups) on the real geography of the assemblages
and calibrates the comparison to the record's small sample sizes and
time-averaging with a signal-recovery experiment. The headline result: within
regions the structure is reproduced by neutral spatial drift and does not
require bounded groups, so the phase, at the scale at which it is drawn, is
largely indistinguishable from drift on geography.

## Quick start

```bash
git clone https://github.com/clipo/phases.git
cd phases
conda env create -f environment.yml      # or see the pip path below
conda activate phases
pip install -e .                          # makes the mls_emergence package importable
./run_all.sh                              # regenerate every figure and result table
```

The full run takes roughly an hour (Monte Carlo simulations); it caches its
heavy results, so re-runs are fast. To regenerate just one figure, run its
script directly (see the [Figures](#figures) table). The manuscript itself is in
`docs/manuscript/` and rebuilds with pandoc (see [Manuscript](#manuscript)).

## Installation

Three ways to get a working environment, in order of reproducibility. All use
Python 3.11 or newer.

### Conda / mamba (recommended)

Conda-forge ships the compiled GEOS/PROJ/GDAL libraries that the geospatial
steps need, so no system packages are required:

```bash
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate phases
pip install -e .                          # makes the mls_emergence package importable
```

For an exactly pinned, per-platform build (linux-64, osx-64, osx-arm64, win-64),
install from the lockfile instead of `environment.yml`:

```bash
pip install conda-lock                              # once
conda-lock install --name phases conda-lock.yml
conda activate phases
pip install -e .
```

### pip + venv

This path needs the system GEOS and PROJ libraries first:

```bash
# macOS:         brew install geos proj
# Debian/Ubuntu: sudo apt-get install libgeos-dev libproj-dev proj-data proj-bin
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .                          # makes the mls_emergence package importable
```

### Docker

A self-contained container (conda environment baked in):

```bash
docker build -t phases .
docker run --rm -it \
  -v "$PWD/output:/repo/output" -v "$PWD/figures:/repo/figures" \
  phases ./run_all.sh
```

### Checking the install

```bash
python3 -c "import mls_emergence, geopandas, cartopy; print('environment OK')"
```

If that prints `environment OK`, every analysis script will run. The repository
is self-contained and has no dependency on any sibling repository.

## Reproducing the analysis

Run everything (outputs go to `output/` and `figures/`):

```bash
./run_all.sh
# or, to use a specific interpreter:
PYTHON=.venv/bin/python ./run_all.sh
```

`run_all.sh` runs the numbered scripts in `analyses/` in ascending order, then
the two shared figure pipelines (`make_map.py`, `make_figures.py`). It continues
past any single failure and reports the count at the end. The full run takes
roughly an hour: several scripts run Monte Carlo simulations with many seeds. The
heaviest are `33_time_aware_emergence.py` and `35_basin_pullout.py` (500 drift
realizations each), `34_emergence_robustness.py` (a 432-run factor grid),
`37_lmv_drift_groups.py` (500 realizations on the wider valley), and
`19_abc_transmission.py` (40,000 ABC simulations). Each caches its per-run
results to `output/`, so re-runs are fast.

The Bayesian analyses (`38`–`41`) need PyMC and arviz, installed by
`requirements.txt`. `39_abc_smc_validation.py` is among the heaviest (a few
hundred short MCMC fits, parallelized across cores). PyMC/PyTensor numerical
output can shift across library versions, so the pinned versions in
`requirements.txt` reproduce the deposited MCMC result tables exactly; a
different PyMC version may change the trailing digits.

To regenerate a single figure or result, run its script directly. For example:

```bash
python3 analyses/36_canonical_phase_map.py    # Figure 1 (phase-territory map)
python3 analyses/21_signal_recovery.py        # Figures 4 and 5 (signal recovery)
python3 analyses/30_regional_map.py           # Figure 10 (regional map)
```

## Figures

Every manuscript figure is produced by one script and written to `figures/` in
four formats (`.png`, `.pdf`, `.svg`, and a local-only `.tif`). Run the script in
the table to regenerate that figure.

| Figure | Script | Output (in `figures/`) |
|---|---|---|
| Fig 1  | `36_canonical_phase_map.py`   | `fig1_phases` |
| Fig 2  | `29_concept_figure.py`        | `fig2_concept` |
| Fig 3  | `make_figures.py`             | `fig3_ca_ordination` |
| Fig 4  | `21_signal_recovery.py`       | `fig4_recovery` |
| Fig 5  | `21_signal_recovery.py`       | `fig5_empirical_trajectory` |
| Fig 6  | `make_figures.py`             | `fig6_idss_structure` |
| Fig 7  | `make_figures.py`             | `fig7_ranksize` |
| Fig 8  | `37_lmv_drift_groups.py`      | `fig8_lmv_drift_groups` |
| Fig 9  | `35_basin_pullout.py`         | `fig9_parkin_pullout` |
| Fig 10 | `30_regional_map.py`          | `fig10_regional` |
| Fig 11 | `31_within_region_structure.py` | `fig11_within_region` |
| Fig S1 | `make_figures.py`             | `figS1_validation` |
| Fig S2 | `13_neiman_distance_and_fit.py` | `figS2_neiman` |
| Fig S3 | `25_drift_vs_groups_demo.py`  | `figS3_drift_vs_groups` |
| Fig S4 | `33_time_aware_emergence.py`  | `figS4_emergent_phases` |
| Fig S5 | `34_emergence_robustness.py`  | `figS5_emergence_robustness` |
| Fig S6 | `20_tempo_mode_ews.py`        | `figS6_dynamic` |
| Fig S7 | `11_chronology_14c.py`        | `figS7_chronology` |

`make_figures.py` produces Figures 3, 6, 7, and S1 (the house-style validation
and basin diagnostics). The remaining figures are produced by their own numbered
scripts. `figures/` also holds earlier-stage and diagnostic images (`05_*`–`09_*`,
`figX_*`, `fig1_studyarea`, `fig5_idss_network`, `figS9_canonical_phases`) kept
for provenance; no main-text result depends on them.

## Generated output

Running the analyses populates two directories. Both are regenerated from source
and can be deleted and rebuilt at any time.

- **`figures/`** — the manuscript figures (above) plus the retained
  earlier-stage images.
- **`output/`** — text and data results:
  - **Result reports** (`*.md`), one per analysis, holding the numbers cited in
    the paper. Examples: `validation_report.md` (criterion GO/NO-GO),
    `basin_results.md` (consolidated basin numbers), `signal_recovery.md`
    (size-controlled recovery), `chronology_14c.md` (radiocarbon),
    `lmv_drift_groups.md` and `emergence_robustness.md` (drift simulations).
    These reports are version-controlled as the canonical reproducibility
    baseline. Every Monte Carlo analysis is seeded, so rerunning the pipeline
    reproduces them exactly and `git diff output/*.md` stays empty.
  - **Monte Carlo caches** (`*.csv`, `*.npz`) so repeat runs are fast:
    `emergence_robustness.csv`, `basin_pullout_runs.csv`, `basin_pullout_prob.csv`,
    `lmv_drift_groups_runs.csv`, `time_aware_runs.csv`, `abc_posterior.npz`. These
    are regenerable and are not tracked. Delete them to force a full from-scratch
    recompute (e.g. `rm output/*.csv output/*.npz`, then `./run_all.sh`); the
    result reports above must come back identical.
  - `diag_macro_transect.png` — the between-region (central-vs-lower valley)
    transect diagnostic, retained for provenance; it is not a manuscript figure
    (the paper reports the within-region drift result for each phase scheme).

## Repository layout

```
analyses/         numbered analysis and figure scripts (each has a header docstring)
src/              the mls_emergence package (transmission models, signatures, I/O)
data/             source ceramic tables, radiocarbon, and geospatial layers (provenance in data/README.md)
docs/manuscript/  the manuscript (MAIN_TEXT.md, SUPPLEMENTAL_TEXT.md, references.bib, american-antiquity.csl, built docx/PDF)
figures/          manuscript figures (committed); run_all.sh regenerates them
output/           generated result tables and Monte Carlo caches (created by run_all.sh)
run_all.sh        reproduce every analysis and figure
environment.yml / conda-lock.yml / requirements.txt / Dockerfile / postBuild   reproducible environments
pyproject.toml    package metadata for `pip install -e .`
```

## Script index

Every script begins with a docstring describing what it does and what it writes.
The numbered scripts run in order under `run_all.sh`.

**Basin definition and chronology**

| Script | Purpose |
|---|---|
| `14_drainage_basin.py`, `16_basin_membership.py`, `09_parkin_basin_restricted.py` | Define the St. Francis basin assemblage set (drainage / watershed rules) |
| `11_chronology_14c.py` | Calibrate the basin radiocarbon and test the seriation-vs-calendar order (Fig S7) |

**Basin transmission analyses**

| Script | Purpose |
|---|---|
| `10_neiman_size_diagnostic.py` | Neiman (1995:17) sample-size / time-averaging diagnostic of the neutrality estimators |
| `13_neiman_distance_and_fit.py` | Neiman diversity-distance diagnostic (Fig S2) |
| `15_continuum_test.py` | Clustering-tendency tests on the CA ordination |
| `17_basin_results.py` | Consolidated basin Results numbers |
| `18_kandler_shennan_neutrality.py` | Non-equilibrium neutrality test |
| `19_abc_transmission.py` | Time-averaging-aware ABC inference of the copying bias |
| `20_tempo_mode_ews.py` | Tempo-and-mode model comparison and early-warning probe (Fig S6) |
| `12_sensitivity_grid.py` | Forking-paths sensitivity of the basin result |

**Bayesian robustness analyses**

| Script | Purpose |
|---|---|
| `38_abc_smc_transmission.py` | ABC-SMC (sequential Monte Carlo) posterior of the copying bias, the calibrated upgrade to `19` |
| `39_abc_smc_validation.py` | Simulation-based calibration, coverage, and cross-check of the ABC-SMC posterior against `19` |
| `40_hierarchical_convergence.py` | Hierarchical Bayesian posterior of the four-signature convergence: `P(all four slopes rise together)` |
| `41_hierarchical_convergence_validation.py` | Calibration, prior sensitivity, and posterior-predictive checks for the convergence model |
| `43_bayesian_fst.py` | Balding-Nichols Bayesian cultural F_ST for the basin: the BN F_ST parameter, a conjugate Gini-Simpson readout, and a structure-vs-panmixia Bayes factor |
| `44_bayesian_fst_validation.py` | Coverage, simulation-based calibration, small-sample-bias, MCMC-health, and Bayes-factor-sanity checks for the BN F_ST |

These use PyMC and arviz (in `requirements.txt`). They restate the transmission-bias inference and the convergence criterion as calibrated Bayesian posteriors, and reach the same conclusions as the primary analyses: a copying bias tightly constrained near neutral, and no convergence of the four signatures (`P = 0`).

`43`/`44` fit the between-cluster cultural F_ST with the Balding-Nichols model (`F ~ Uniform(0,1)`, `alpha = (1-F)/F`, per-cluster frequencies marginalized inside a Dirichlet-multinomial), the same model the hyperlocality project uses, so the two analyses share one Bayesian approach. The flat prior is on F_ST itself, so the estimate is not pulled toward the drift level, and the marginalized likelihood is well conditioned (no funnel-taming `target_accept` needed). One fit yields two readouts: the BN F_ST parameter (basin median 0.07, 95% CI [0.03, 0.11]) and a conjugate reconstruction of the exact Gini-Simpson estimator the manuscript reports (`p_g | x_g ~ Dirichlet(alpha*pi + x_g)`; basin median 0.033, 95% CI [0.030, 0.035], matching the frequentist plug-in). This is estimation uncertainty conditional on the clustering, distinct from the size-controlled per-bin F_ST that carries the main-text argument. The reported Bayes factor tests structure against *exact* panmixia, a far weaker null than drift; with large samples panmixia is rejected trivially, so it is a model-adequacy and parity check, not the structure test. The operative structure-vs-drift test is the separate stochastic-drift null.

**Criterion validation and signal recovery**

| Script | Purpose |
|---|---|
| `04_validation_report.py` | Blind validation of the convergence criterion on synthetic mechanisms |
| `21_signal_recovery.py` | Record-matched, size-controlled signal-recovery calibration (Figs 4, 5) |
| `22_generator_diagnostic.py` | Which signatures are recoverable at the record's resolution |

**Drift versus bounded groups (the central test)**

| Script | Purpose |
|---|---|
| `23_phases_vs_spatial_drift.py`, `24_phases_drift_robustness.py`, `25_drift_vs_groups_demo.py` | Drift-versus-bounded-groups generative comparison on the real coordinates (Fig S3) |
| `33_time_aware_emergence.py` | Basin positive control: phase-like structure from neutral drift, sampled time-transgressively (Fig S4) |
| `34_emergence_robustness.py` | Robustness of the drift result across interaction length, mixing, and innovation (Fig S5) |
| `35_basin_pullout.py` | Parkin pull-out: does the Parkin phase separate from the other named phases by more than drift (Fig 9) |
| `37_lmv_drift_groups.py` | Phase-like group counts from neutral drift across the wider lower Mississippi Valley (Fig 8) |

**Regional extension (southeast-Missouri phases)**

| Script | Purpose |
|---|---|
| `26_cmv_phase_groupness.py`, `27_cmv_lmv_repertoire.py` | Williams (1954) southeast-Missouri repertoire and within-region groupness |
| `28_macro_boundary.py` | Between-region transect diagnostic (retained for provenance, not a manuscript figure) |
| `31_within_region_structure.py` | Within-region structure for each phase scheme, lower valley vs central valley (Fig 11) |

**Maps and figure pipelines**

| Script | Purpose |
|---|---|
| `29_concept_figure.py` | Conceptual figure: what the test distinguishes (Fig 2) |
| `30_regional_map.py` | Regional two-scale map (Fig 10) |
| `36_canonical_phase_map.py` | Canonical phase-territory map (Fig 1) |
| `make_figures.py` | House-style pipeline for Figures 3, 6, 7, and S1, plus shared plotting helpers |
| `make_map.py` | Shared map helpers, the river basemap, and river-network distance |
| `figstyle.py` | Figure house style (fonts, palette, multi-format export) |

Scripts `05`–`07` are earlier-stage analyses retained for provenance.

## Data

Source data live in `data/`, with full provenance, licensing, and per-file
descriptions in [`data/README.md`](data/README.md). In brief: the decorated
sherd counts come from Phillips, Ford & Griffin (1951) as compiled by Lipo
(2001) and from Williams (1954) for southeast Missouri; the radiocarbon is the
Mainfort (2001) compilation; the geospatial layers are the project's Lower
Mississippi Valley Survey shapefiles plus a clip of HydroRIVERS for the regional
maps. State boundaries and the locator insets are drawn from Natural Earth, which
cartopy downloads on demand.

## Manuscript

The paper is in `docs/manuscript/`: `MAIN_TEXT.md` and `SUPPLEMENTAL_TEXT.md`
(sources), `references.bib` (bibliography), and the built `.docx` and `.pdf` for
each. Rebuild with pandoc:

```bash
pandoc docs/manuscript/MAIN_TEXT.md -o docs/manuscript/MAIN_TEXT.pdf \
  --citeproc --bibliography=docs/manuscript/references.bib \
  --csl docs/manuscript/american-antiquity.csl \
  --pdf-engine=xelatex -H docs/manuscript/_pdf_header.tex \
  --resource-path=docs/manuscript:.
```

Build the Supplemental Text the same way, substituting `SUPPLEMENTAL_TEXT` for
`MAIN_TEXT`, and produce the `.docx` by changing the output extension. The
American Antiquity citation style (`american-antiquity.csl`) is vendored in
`docs/manuscript/`, so the build is self-contained.

## Licensing

- **Code** (`src/`, `analyses/`): MIT License, see [`LICENSE`](LICENSE).
- **Data** (`data/`): Creative Commons Attribution 4.0 International (CC BY 4.0),
  see [`LICENSE-data`](LICENSE-data), which also lists data provenance.

## Citing

If you use this code or data, please cite the paper (in preparation) and this
repository:

> Lipo, Carl P., Robert J. DiNapoli, and Mark E. Madsen. Are the phases real?
> Distinguishing bounded interaction groups from spatially structured drift in
> lower Mississippi Valley decorated ceramics.

The archived version of the code and data is deposited at Zenodo: https://doi.org/10.5281/zenodo.20796240 (private until the paper is published, so it will not resolve before then).

## Authors

Carl P. Lipo, Robert J. DiNapoli, and Mark E. Madsen.
