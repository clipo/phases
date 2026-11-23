# What are phases?

Code and data for **"What Are Phases? Geographic Structure, Neutral Copying,
and the Search for Social Boundaries in the Lower Mississippi Valley"** by
Carl P. Lipo, Beau DiNapoli, and Mark E. Madsen.

The repository holds the decorated-ceramic counts, site coordinates, rim
attribute counts and radiocarbon data the paper uses; every analysis script
with the findings file it writes; the Stan programs for the spatial field,
latent-unit and barrier models; the manuscript figures and the scripts that
draw them; and the manuscript sources with their built Word and PDF files.
`MANIFEST.md` names the script that produces each result and figure. The
numbers in the paper come from the findings files under `output/findings/`.

**Contents:** [The questions](#the-questions) ·
[Three levels of reproduction](#three-levels-of-reproduction) ·
[Quick start](#quick-start) · [Figures and tables](#figures-and-tables) ·
[Reproducing the analyses](#reproducing-the-analyses) ·
[Use the tests on your own data](#use-the-tests-on-your-own-data) ·
[Installation](#installation) · [Data](#data) · [Manuscript](#manuscript) ·
[Repository layout](#repository-layout) · [Licensing](#licensing) ·
[Citing](#citing)

## The questions

Culture historians built phases to describe continuous spatial and temporal
variation in artifact classes, yet researchers still often treat phases as
social units and communities. The paper asks three questions about the late
pre-contact phases of the St. Francis basin in the Lower Mississippi Valley:

1. Why did culture historians see phases in these data?
2. Is there evidence of social boundaries between the phases beyond the
   geographic distribution of settlement?
3. What would it take to detect such boundaries in a record like this one?

The answers, on 28 time-averaged decorated-ceramic assemblages:

1. **Geography.** Clustering on site coordinates alone recovers Phillips's
   phases for all but two assemblages. A spatial field model assigns 93
   percent of the between-site variance in composition to a smooth geographic
   surface, and the phase boundaries run through the least-occupied ground
   between settlement clusters.
2. **No evidence of social boundaries.** The difference between the basin's
   two spatial clusters is one that neutral copying at matched diversity can
   produce, although most settings produce less. Three common decorated types,
   and rim attributes tabulated separately on 19 sites, vary between sites
   more than the tested neutral settings produce. The two share no pattern and
   support no phase boundary.
3. **Data this record lacks.** Neither of the two models built to detect
   bounded units reliably recovered known units on records simulated at these
   sites, plausibly because of the gaps between settlement clusters. Detecting
   a boundary would take assemblages from those gaps, dated or stratified
   contexts, and trait systems recorded separately from the decorated types.

The phases are real as a geographic pattern, not resolved as bounded units by
the two models tested, and not demonstrated as social units.

## Three levels of reproduction

The pipeline has three levels. Each level reads what the one before it
committed, so a reader can stop at any of them.

1. **Figures from the committed findings (minutes, Python only).** The eight
   main-text figures and Figure S1 are drawn from the findings files under
   `output/findings/` and from the data. No model is refitted and nothing is
   simulated. Each figure script asserts the values the text quotes before it
   draws, so a figure cannot silently drift from a changed findings file. This
   is what `./run_all.sh` does by default.
2. **The Python analyses that write those findings files (hours).** The phase
   recovery and settlement analyses, the minimum sample size, the radiocarbon
   chronology, the Mainfort (2003a) comparison, the neutral copying
   simulations (the resemblance screen, the per-class tables with regional
   drift and planted markers, the composition-by-collection analysis, the
   breakdown of the painted class by type, the two
   further sets of simulations of Tables S7 and S8), the rim analysis, and the
   Stan data contracts. `./run_all.sh --analyses` runs them.
3. **The Stan fits (R, rstan, a C++ toolchain; hours to days).** The spatial
   field model on the basin with its prior checks and recovery runs, the
   latent-unit model and the barrier model on records from their own
   generators, and the prior and range sensitivity refits. File names and
   findings files call these three models A0, A1 and A3. `./run_all.sh --stan` runs them
   after the Python analyses. The fits are cached by content in
   `output/gp_cache/`, which is not published; the findings files carry every
   posterior summary together with its sampler diagnostics (R-hat, bulk and
   tail effective sample size, divergences, treedepth saturation, E-BFMI).

The rim-attribute analysis (`analyses/rims/rim_r1.py`, Figure 7) runs at level
2. It covers three attributes: interior beveling, exterior notching and lip
form. Its lip-form counts come from an unpublished tabulation by Robert
Mainfort (`data/raw/mainfort_miss_rims.xls`), the counts for beveling and
notching come from his published Table 1
(`data/raw/mainfort2003_rim_table1_counts.csv`), and
`scripts/build_mainfort_rim_counts.py` combines the two into
`data/processed/mainfort_rim_counts.csv`. See `data/README.md`. The tabulation
also records paste (Bell or Neeley's Ferry). The paper treats paste as a
technological attribute and does not analyze it, so the findings file lists
the paste counts and reports no statistic for them.

The neutral copying model's scripts and the rim analysis read simulator time
scales from `output/probe_latent_groups/`, which is not published. Level 2
writes them first (`analyses/probe_latent_groups/time_scales.py`, under a
minute), so run that script before running any of the others alone.

## Quick start

```bash
git clone https://github.com/clipo/phases.git
cd phases
conda env create -f environment.yml      # or the pip path under Installation
conda activate phases
pip install -e .                          # makes the mls_emergence package importable
./run_all.sh                              # level 1: Figures 1 to 8 and S1 from the committed findings
```

To redraw one figure, run its script (see [Figures and tables](#figures-and-tables)).
The manuscript is in `docs/manuscript/` and rebuilds with pandoc
(see [Manuscript](#manuscript)).

## Figures and tables

Each figure script writes `.png`, `.pdf`, `.svg` and a local-only `.tif` with
one stem under `figures/`. The `.png` is what the manuscript embeds.

| Figure | Script | Output stem | Reads |
|---|---|---|---|
| 1 | `analyses/93_phillips1970_phase_map.py` | `fig1_phillips1970_phases` | `data/raw/phillips1970_phase_outlines.csv`, the analysis matrix |
| 2 | `analyses/102_fig2_two_histories.py` | `fig2_two_histories` | nothing; an illustrative schematic from a fixed seed |
| 3 | `analyses/103_fig3_model_schematic.py` | `fig3_model_schematic` | the analysis matrix (site coordinates); a schematic from a fixed seed |
| 4 | `analyses/104_fig4_phase_recovery_map.py` | `fig4_phase_recovery_map` | `output/findings/phase_recovery.md`, `phase_recovery.csv`, `boundaries_in_settlement_gaps.md` |
| 5 | `analyses/105_fig5_spatial_field.py` | `fig5_spatial_field` | `output/findings/a0_intercept_prior.md` |
| 6 | `analyses/106_fig6_per_class_drift.py` | `fig6_per_class_drift` | `output/findings/per_type_drift_sensitivity.md` |
| 7 | `analyses/107_fig7_rim_attributes.py` | `fig7_rim_attributes` | `output/findings/rim_r1.md` |
| 8 | `analyses/108_fig8_latent_unit_recovery.py` | `fig8_latent_unit_recovery` | `output/findings/a1_validation.md`, `a3_validation.md` |
| S1 | `analyses/11_chronology_14c.py` | `figS1_chronology` | `data/raw/14CDatesFromMainfort2001.csv`, `data/raw/intcal20.14c` |
| S2 | `analyses/make_figures.py` | `figS2_ranksize` | the settlement tables under `data/` |

The tables of the main text and supplement come from these
findings files:

| Table | Findings file | Produced by |
|---|---|---|
| 2 | `output/findings/per_type_drift_sensitivity.md`, `per_type_record_row.md` | `analyses/per_type/drift_sensitivity_figure.py` |
| S1 | `output/findings/minimum_sample_size.md` | `analyses/82_minimum_sample_size.py` |
| S2 | `output/findings/phase_recovery.md`, `phases_under_drift.md` | `analyses/76_phase_recovery.py`, `analyses/79_phases_under_drift.py` |
| S3 | `output/findings/settlement_and_rivers.md` | `analyses/81_settlement_and_rivers.py` |
| S4 | `output/findings/a0_basin_refit.md` | `analyses/probe_latent_groups/refit_a0_basin.R` |
| S5 | `output/findings/a3_validation.md`, `a3_pinned_rho6.md`, `a3_pinned_rho12.md` and their `_alphaN0-2` siblings | `validate_a3.R`, `a3_pinned_range.R` |
| S6 | `output/findings/per_type_drift_sensitivity.md` | `analyses/per_type/drift_sensitivity_figure.py` |
| S7 | `output/findings/type_origination.md` | `analyses/per_type/type_origination.py` |
| S8 | `output/findings/type_origination_calibration.md` | `analyses/per_type/type_origination_calibration.py` |
| S9, S10 | `output/findings/mainfort_replication.md` | `analyses/83_mainfort_replication.py` |

Tables 1 and 3 of the main text are prose tables with no computed cells.

`figures/` also carries the figures of the earlier version of this paper
(`fig1_phases`, `fig9_parkin_pullout`, `fig12_` to `fig18_`, `figS1_` to
`figS8_` and others) and `analyses/` the numbered scripts that drew them. They
are retained because several of those scripts still produce findings files the
current paper reads (76, 79, 80, 81, 82, 83), because `CORRECTIONS.md` and the
git history refer to them, and because the earlier version is archived at
Zenodo under this repository's name. No result of the current paper depends on
a figure outside the two tables above.

## Reproducing the analyses

### Level 1: figures

```bash
./run_all.sh
# or, to use a specific interpreter:
PYTHON=.venv/bin/python ./run_all.sh
```

### Level 2: the Python analyses

```bash
./run_all.sh --analyses
```

This runs, in order: `scripts/build_analysis_matrix.py` (the decorated-class
matrix and its provenance file); `analyses/16_basin_membership.py` (the
28-assemblage basin set); `analyses/82_minimum_sample_size.py`;
`analyses/76_phase_recovery.py`, `79_phases_under_drift.py`,
`80_boundaries_in_settlement_gaps.py` and `81_settlement_and_rivers.py`;
`analyses/11_chronology_14c.py`; `scripts/build_mainfort2003_sites.py` and
`analyses/83_mainfort_replication.py`; `analyses/49_write_stan_data.py` (the
Stan data contracts under `data/stan/`); then
`analyses/probe_latent_groups/time_scales.py` (the simulator time scales the
later scripts read); then the neutral copying simulations:
`analyses/probe_latent_groups/resemblance_screen.py` (the
297-cell screen at six seeds for each copying-network family, 2,970 cells in
all), `resemblance_stage2.py`, `per_class_resemblance.py`,
`site_fst_decomposition.py`, the per-class tables with regional drift and
planted markers (`analyses/per_type/per_class_table.py`, generators C3 and M
at the two operating points, then `drift_sensitivity_figure.py`),
`analyses/per_type/composition_by_collection.py`,
`analyses/per_type/painted_by_type.py` (the painted class by type in the 1940
collections), the two further sets of
simulations (`drift_extension.py` and `type_origination.py` with their
`_joint` and `_calibration` scripts, Tables S7 and S8), and the rim analysis
(`scripts/build_mainfort_rim_counts.py`, then `analyses/rims/rim_r1.py` at
both operating points and with `--readings`). The simulations are the
slow part; they use `multiprocessing` and every record is seeded from a
`SeedSequence` named in the script's header. The operating points are registered in
`analyses/probe_latent_groups/operating_points.json`.

**A seeded rerun reproduces a committed findings file only under the NumPy
version that produced it.** NumPy 2.5 changed the draws
`Generator.multinomial` returns from a given state. The simulator time
scales (`analyses/probe_latent_groups/time_scales.py`,
`output/findings/simulator_time_scales.md`) and the per-class,
regional-drift and rim tables (`analyses/per_type/` and `analyses/rims/`)
were produced under NumPy 2.5.4. The later scripts read the time scales, so
a run under another version differs in its draws, and the scripts that read
the time scales stop with an error because the regenerated values do not
match the registered operating points. **Level 2 therefore needs the NumPy
2.5.4 environment below.** Level 1 and the test suite run under the pinned
environment. The simulation results carried from the
earlier version of the study were produced under NumPy 2.4.6, which is the
version `requirements-lock.txt` pins. The two sets have not been regenerated
under one version. To reproduce the first set value for value, run those
scripts in a second environment that differs only in NumPy:

```bash
python -m venv .venv254
.venv254/bin/pip install -r requirements.txt
.venv254/bin/pip install --no-deps "numpy==2.5.4"     # after the requirements, because they pin numpy below 2.5
.venv254/bin/pip install --no-deps -e .
PYTHON=.venv254/bin/python ./run_all.sh --analyses
```

Checked on 2026-11-21 for the time scales and the rim analysis: under NumPy
2.5.4, `time_scales.py` reproduces `simulator_time_scales.md`, and
`analyses/rims/rim_r1.py` at both operating points and with `--readings`
reproduces every value of `output/findings/rim_r1.md`. Under NumPy 2.4.6,
`time_scales.py` gives a local time scale of 171 steps at the primary
operating point where the committed file has 182, and the rim commands give
different neutral bands. The per-class and regional-drift
tables were not rerun on that date. Scripts that draw no random numbers, such
as the data products and `painted_by_type.py`, do not depend on the version.

### Level 3: the Stan fits

The models are Stan programs under `stan/spatial/`: `composition_gp_marginal.stan`
(the spatial field model, A0), `composition_gp_groups.stan` (the latent-unit
model, A1) and `composition_gp_barrier.stan` (the barrier model, A3), sharing
`kernels.stanfunctions`.
They are fitted through `R/fit_spatial.R`, which attaches the diagnostics to
every fit and applies the diagnostic gate the paper states (R-hat at most
1.01, bulk and tail effective sample size at least 400, no divergences,
E-BFMI at least 0.3).

```bash
Rscript analyses/00_setup/00_install_deps.R        # rstan, posterior, bayesplot, jsonlite, here, digest, ggplot2, systemfonts, ragg
Rscript analyses/00_setup/01_verify_toolchain.R    # compiles and samples a model with a known posterior
./run_all.sh --stan
```

`--stan` runs `analyses/probe_latent_groups/refit_a0_basin.R` (the spatial
field model on the basin and the prior checks of Table S4), `a0_intercept_prior.R` (the three
intercept priors of Figure 5), `a0_pinned_nuisances.R`, `a1_contracts.py` and
`a3_contracts.py` (the simulated validation records, which take their field
scale and range from the basin fit), `validate_a1.R` and `validate_a3.R` (the
recovery tests of Figure 8), `a3_pinned_range.R` at 6 and 12 km (Table S5),
the same record writers and validation scripts again with `--alpha-sd=2`
(the Normal(0, 2) intercept prior: the `*_alphaN0-2.md` tables and the
Normal(0, 2) column of Table S5, on records kept in `data/stan/alphaN0-2/`),
and then the Python summaries that read the fit records,
`a0_predictive_summaries.py`, `a3_prior_comparison.py` and
`export_fit_diagnostics.py`. Each script sets `options(mc.cores = 4L)`. The
committed Stan findings were produced on 2026-11-22, and the Normal(0, 2)
validation tables on 2026-11-23, with R 4.6.2, rstan 2.32.7 and StanHeaders
2.39.1.1 on macOS (arm64), as the header of each findings file records. Sampler output differs in the second or third digit
between R and StanHeaders versions at the same seed, so a rerun reproduces
the conclusions and the gate outcomes of most fits, not the files.
`a1_prior_check.md`, `a0_basin_refit_alphaN0-2.md`,
`a0_predictive_summaries_alphaN0-2.md` and `fit_diagnostics_2026-10.md` are
the fits of October 2026 and are not rerun by `--stan`.

### Tests

```bash
pip install -e ".[dev]"                                       # pytest is in the dev extras, not in requirements.txt or environment.yml
pytest                                                        # Python: the package, the data loaders, the probe, the per-class and rim code, phasecheck
Rscript -e 'testthat::test_file("tests/testthat/test-fit-layer.R")'   # R: the diagnostics layer on synthetic draws with analytic answers
Rscript -e 'testthat::test_file("tests/testthat/test-kernels.R")'
```

Some tests are skipped in a fresh clone, each with its reason printed by
`pytest -rs`: two that read the simulator time scales (generated at level 2
and not committed), four that compare against files of the authors' audit
workspace (not published), two that need the parent model package of the
earlier version of the study (`monument-mls`, see `pyproject.toml`), and three
that compare against tables produced under NumPy 2.5.4. Checked in a fresh
clone on 2026-11-21 under the pinned environment: 441 passed, 11 skipped, 1
failed, the failure being the known issue below.

Known issue: `tests/phasecheck/test_copying_model.py::test_model_comparison_tells_a_restricted_record_from_an_unrestricted_one`
fails on 3 of its 80 records under numpy 2.4.6, scipy 1.18.1 and pandas 3.0.6.
The code has not changed since the test passed; the difference is in the
numerical libraries, and it is under investigation. Under numpy 2.5, two more
tests in that file fail because seeded results change with the numpy version
(checked 2026-11-14 with numpy 2.5.4). `pyproject.toml` therefore pins numpy
below 2.5.

## Use the tests on your own data

`phasecheck` applies the phase tests of the earlier version of this paper to
any table of assemblages. It needs one row per assemblage: a name, latitude and
longitude, a phase, and one column of sherd counts per class.

```bash
pip install -e .
phasecheck my_assemblages.xlsx --out my_report            # the phase lines against other lines
phasecheck my_assemblages.xlsx --model --jobs 8 --out my_report   # and the copying model with no groups
```

It answers a series of questions in order and writes a report: does location
alone recover the phases; do the phase boundaries separate the assemblages
better than other lines; does each assemblage fit its own phase; and, with
`--model`, does the record differ between phases by more than copying across
distance produces, and could the comparison detect a boundary on a record like
yours. The tool implements those tests, not the spatial field, latent-unit and
barrier models or the diversity-matched neutral comparison of the current main
text. The paper's own record is the worked example in `examples/phasecheck/`.
The guide is [docs/PHASECHECK.md](docs/PHASECHECK.md).

## Installation

The Python pipeline (levels 1 and 2) runs on Windows, macOS (Intel and Apple
Silicon) and Linux with Python 3.11 or newer. The geospatial stack (geopandas,
cartopy, pyproj, shapely) is what makes the installation path matter:
conda-forge ships the compiled GEOS, PROJ and GDAL libraries on all three
operating systems, so the conda path needs no system packages. The pip path
needs those libraries installed separately, which is straightforward on macOS
and Linux and awkward on native Windows. On Windows, use conda or WSL2.

### Conda (recommended)

```bash
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate phases
pip install -e .
```

For an exactly pinned, per-platform build, install from the lockfile instead:

```bash
pip install conda-lock
conda-lock install --name phases conda-lock.yml
conda activate phases
pip install -e .
```

### pip and venv

```bash
# macOS:         brew install geos proj
# Debian/Ubuntu: sudo apt-get install libgeos-dev libproj-dev proj-data proj-bin
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt           # or requirements-lock.txt for the full pinned closure
pip install -e .
```

### R and Stan (level 3)

Level 3 needs R (4.3 or newer), a C++ toolchain (`build-essential` on Ubuntu;
Xcode command line tools on macOS; Rtools on Windows) and the packages that
`analyses/00_setup/00_install_deps.R` installs. rstan compiles the models on
first use and caches them under `output/stan_cache/`.

### Container and Binder

`Dockerfile` builds an image with the conda environment and pandoc; it covers
levels 1 and 2 and the manuscript build, not R. Mount `output/` and `figures/`
to keep what a run writes:

```bash
docker build -t phases .
docker run --rm -it -v "$PWD/output:/repo/output" -v "$PWD/figures:/repo/figures" phases ./run_all.sh
```

A prebuilt image of the earlier release is at `ghcr.io/clipo/phases:1.0.0`; it
predates the current scripts and should not be used for this version.
[Binder](https://mybinder.org/v2/gh/clipo/phases/main) builds the conda
environment in the browser for a quick look at level 1.

### Checking the install

```bash
python -c "import mls_emergence, geopandas, cartopy; print('environment OK')"
```

## Data

Source data live in `data/`, with per-file descriptions, sources, the matrix
construction rule and licensing in [`data/README.md`](data/README.md). In
brief: the decorated sherd counts are the Phillips, Ford and Griffin (1951)
survey counts (`data/raw/PFGData.xlsx`) and Lipo's (2001) compilation of those
counts with his 1996 to 1997 collections (`data/raw/mainfort-pfg-cpl.xlsx`),
one source per assemblage and never summed; Mainfort's (2003a) published table
is a separate comparison set; the rim attribute counts for beveling and
notching are Mainfort's (2003b) Table 1, and those for lip form are
from his unpublished tabulation (`data/raw/mainfort_miss_rims.xls`); the radiocarbon is the Mainfort
(2001) compilation calibrated against IntCal20; the settlement tables are the
regional compilation of Lipo and Dunnell (2007) with the published corrections
in `data/raw/`; the geospatial layers are the project's Lower Mississippi
Valley shapefiles and clips of HydroRIVERS. State boundaries and locator insets
come from Natural Earth, which cartopy downloads on demand.

## Manuscript

The paper is in `docs/manuscript/`: `MAIN_TEXT.md` and `SUPPLEMENTAL_TEXT.md`
(sources), `references.bib`, the American Antiquity citation style
(`american-antiquity.csl`), the Word styling (`reference.docx`), and the built
`.docx` and `.pdf` of each. Rebuild from `docs/manuscript/` with pandoc:

```bash
pandoc MAIN_TEXT.md --citeproc --bibliography=references.bib \
  --csl=american-antiquity.csl --reference-doc=reference.docx -o MAIN_TEXT.docx
pandoc MAIN_TEXT.md --citeproc --bibliography=references.bib \
  --csl=american-antiquity.csl --include-in-header=_pdf_header.tex -o MAIN_TEXT.pdf
```

Build the Supplemental Text the same way. Two things differ between a fresh
build and the committed files. The committed `.docx` files were adjusted by
hand in Word after the pandoc build (table widths and similar styling, with
the text edits carried back into the Markdown), so a rebuilt `.docx` has the
same content in pandoc's default layout. The committed `.pdf` files were made
from those Word files with LibreOffice
(`soffice --headless --convert-to pdf MAIN_TEXT.docx`), not with the LaTeX
command above. `docs/manuscript/BUILD.md` has the details, and `scripts/check_figure_claims.py` and `scripts/check_signs.py`
check the figures' printed values and the text's signed values against the
findings files.

## Repository layout

```
analyses/                 numbered analysis and figure scripts; per_type/, probe_latent_groups/ and rims/ hold the current paper's simulations and fits
analyses/00_setup/        R dependency install, toolchain check, environment recording
src/mls_emergence/        the package: transmission (the copying model), probe (generators, readouts, rivals), per_type, rims, signatures, dataio
src/phasecheck/           the phase tests as a tool for other data (docs/PHASECHECK.md)
stan/spatial/             the Stan programs for the spatial field (A0), latent-unit (A1) and barrier (A3) models, and shared kernel functions
R/                        fit_spatial.R (fitting, diagnostics, gate, cache) and diagnostic_plots.R
scripts/                  matrix and table builders, the manifest builder, the figure-claim and sign checks, the release sync
tests/                    pytest suites; tests/testthat/ for R
data/                     source tables and geospatial layers (data/README.md)
output/findings/          the findings files the manuscript cites, one per analysis
output/*.md, *.csv        result reports of the earlier version's pipeline, kept as its baseline
figures/                  manuscript figures; run_all.sh redraws them
docs/manuscript/          the manuscript and its build files
examples/phasecheck/      the worked example for phasecheck
MANIFEST.md               generated: which script writes which output
CORRECTIONS.md            corrections since the first public release
run_all.sh                reproduce by level (figures; --analyses; --stan)
```

## Licensing

- **Code** (`src/`, `analyses/`, `scripts/`, `R/`, `stan/`): MIT License, see [`LICENSE`](LICENSE).
- **Data** (`data/`): Creative Commons Attribution 4.0 International (CC BY 4.0),
  see [`LICENSE-data`](LICENSE-data) and the provenance in `data/README.md`.

## Citing

If you use this code or data, please cite the paper (in preparation) and this
repository:

> Lipo, Carl P., Beau DiNapoli, and Mark E. Madsen. What are phases? Geographic
> structure, neutral copying, and the search for social boundaries in the Lower
> Mississippi Valley.

`CITATION.cff` carries the metadata. The Zenodo deposit at
https://doi.org/10.5281/zenodo.20796240 is version 1.0.0 of this repository
(22 June 2026), the code and data of the earlier version of the paper; a new
version will be deposited when the current paper is published.

## Authors

Carl P. Lipo, Beau DiNapoli, and Mark E. Madsen.
