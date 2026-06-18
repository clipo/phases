# Are the phases real?

Code and data for **"Are the phases real? Distinguishing bounded interaction
groups from spatially structured drift in central Mississippi Valley decorated
ceramics"** by Carl P. Lipo, Robert J. DiNapoli, and Mark E. Madsen.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/clipo/phases/main)

Launch a zero-install environment in the browser with Binder (badge above), or
build locally with conda, the conda-lock lockfile, or Docker (see Installation).

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

This repository tests, for the late pre-contact central Mississippi Valley,
whether decorated-ceramic similarity records bounded interaction groups or
spatially structured drift. It compares two generative processes (neutral
spatial drift versus bounded groups) on the real geography of the assemblages
and calibrates the comparison to the record's small sample sizes and
time-averaging with a signal-recovery experiment. The headline result: within
regions the structure is reproduced by neutral spatial drift and does not
require bounded groups, so the phase, at the scale at which it is drawn, is
largely indistinguishable from drift on geography.

## Manuscript

The paper is in `docs/manuscript/`: `MAIN_TEXT.md` (source), `references.bib`
(bibliography), and the built `MAIN_TEXT.pdf` and `MAIN_TEXT.docx`, with the Supplemental
Text (`SUPPLEMENTAL_TEXT.md` and its built docx/PDF) alongside. It is a
draft in progress. Rebuild it with pandoc:

```bash
pandoc docs/manuscript/MAIN_TEXT.md -o docs/manuscript/MAIN_TEXT.pdf \
  --citeproc --bibliography=docs/manuscript/references.bib \
  --csl docs/manuscript/american-antiquity.csl \
  --pdf-engine=xelatex -H docs/manuscript/_pdf_header.tex \
  --resource-path=docs/manuscript:.
```

Build the Supplemental Text the same way, substituting `SUPPLEMENTAL_TEXT` for
`MAIN_TEXT`. The American Antiquity citation style (`american-antiquity.csl`) is
vendored in `docs/manuscript/` so the build is self-contained.

## Repository layout

```
analyses/         numbered analysis and figure scripts (each has a header docstring)
src/              the mls_emergence package (transmission models, signatures, I/O)
data/             source ceramic tables, radiocarbon, and geospatial reference layers (provenance in data/README.md)
docs/manuscript/  the manuscript (MAIN_TEXT.md, references.bib, american-antiquity.csl, built docx/PDF)
figures/          manuscript figures (committed); run_all.sh regenerates them
output/           generated text/markdown results       (created by run_all.sh)
```

## Installation

Three ways to get a working environment, in order of reproducibility. All use
Python 3.11 or newer.

### Conda / mamba (recommended)

Conda-forge ships the compiled GEOS/PROJ/GDAL libraries, so no system packages
are needed:

```bash
conda env create -f environment.yml     # or: mamba env create -f environment.yml
conda activate phases
pip install -e .                         # makes the mls_emergence package importable
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
pip install -e .          # makes the mls_emergence package importable
```

### Docker

A self-contained container (conda environment baked in):

```bash
docker build -t phases .
docker run --rm -it \
  -v "$PWD/output:/repo/output" -v "$PWD/figures:/repo/figures" \
  phases ./run_all.sh
```

## Reproducing the analysis

Run everything (outputs go to `output/` and `figures/`):

```bash
./run_all.sh
# or, to use a specific interpreter:
PYTHON=.venv/bin/python ./run_all.sh
```

Or run any single step directly, for example:

```bash
python3 analyses/21_signal_recovery.py        # record-matched signal recovery
python3 analyses/23_phases_vs_spatial_drift.py # drift-vs-groups comparison
python3 analyses/29_concept_figure.py          # the conceptual figure
```

The full run takes roughly an hour: several scripts run Monte Carlo simulations
with many seeds. The heaviest are `33_time_aware_emergence.py` (500 drift
realizations) and `34_emergence_robustness.py` (a 432-run factor grid); each
caches its per-run results to `output/`, so re-runs are fast. The run continues
past any single failure and reports the count at the end.

## Script index

Every script begins with a docstring describing what it does and what it writes.
The ones that produce the paper's results and figures:

| Script | Purpose |
|---|---|
| `09`, `14`, `16` | Define the St. Francis basin assemblage set (drainage / watershed rules) |
| `11_chronology_14c.py` | Calibrate the basin radiocarbon and test the seriation-vs-calendar order |
| `13_neiman_distance_and_fit.py` | Neiman diversity-distance diagnostic |
| `15_continuum_test.py` | Clustering-tendency tests on the CA ordination |
| `17_basin_results.py` | Consolidated basin Results numbers |
| `18_kandler_shennan_neutrality.py` | Non-equilibrium neutrality test |
| `19_abc_transmission.py` | Time-averaging-aware ABC inference of the copying bias |
| `20_tempo_mode_ews.py` | Tempo-and-mode model comparison and early-warning probe |
| `21_signal_recovery.py` | Record-matched, size-controlled signal-recovery calibration |
| `22_generator_diagnostic.py` | Which signatures are recoverable at the record's resolution |
| `23`–`25` | Drift-versus-bounded-groups generative comparison on real coordinates |
| `26`–`28` | Central-valley (Williams 1954) vs lower-valley repertoire and macro transect |
| `29_concept_figure.py` | Conceptual figure (what the test distinguishes) |
| `30_regional_map.py` | Regional two-scale map |
| `31_within_region_structure.py` | Within-region structure, central valley vs lower valley |
| `33_time_aware_emergence.py` | Basin positive control: phase-like structure from neutral drift on the river network, sampled time-transgressively, and the probability of sharing Parkin's drift group (Figure S4) |
| `34_emergence_robustness.py` | Robustness of the drift result across interaction length, mixing, and innovation (Figure S5) |
| `35_basin_pullout.py` | Parkin pull-out: does the Parkin phase separate from the other named phases by more than drift (Figure 9) |
| `37_lmv_drift_groups.py` | Phase-like group counts from neutral drift across the wider lower Mississippi Valley (Figure 8) |
| `make_figures.py` | House-style pipeline for most main-text and supplemental figures |
| `make_map.py` | Shared map helpers, the river basemap, and river-network distance |

Scripts `04`–`08` and `10` are earlier-stage analyses retained for provenance.
Scripts `08` and `10` are coupled to the parent monument-signaling model and
require the optional `monument-mls` package (below); no main-text figure depends
on them.

## Optional parent model (monument-mls)

Two supplemental analyses couple the ceramic simulator to a parent
costly-signaling model of monumental construction. They need the `signaling`
package from the separate `monument-mls` repository:

```bash
pip install -e ../monument-mls
```

Without it, those two scripts exit with a clear message and everything else runs
normally.

## Licensing

- **Code** (`src/`, `analyses/`): MIT License, see [`LICENSE`](LICENSE).
- **Data** (`data/`): Creative Commons Attribution 4.0 International (CC BY 4.0),
  see [`LICENSE-data`](LICENSE-data), which also lists data provenance.

## Citing

If you use this code or data, please cite the paper (in preparation) and this
repository:

> Lipo, Carl P., Robert J. DiNapoli, and Mark E. Madsen. Are the phases real?
> Distinguishing bounded interaction groups from spatially structured drift in
> central Mississippi Valley decorated ceramics.

## Authors

Carl P. Lipo, Robert J. DiNapoli, and Mark E. Madsen.
