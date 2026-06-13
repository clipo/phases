# Are the phases real?

Code and data for **"Are the phases real? Distinguishing bounded interaction
groups from spatially structured drift in central Mississippi Valley decorated
ceramics"** by Carl P. Lipo and Robert J. DiNapoli.

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
(bibliography), and the built `MAIN_TEXT.pdf` and `MAIN_TEXT.docx`. It is a
draft in progress. Rebuild it with pandoc:

```bash
pandoc docs/manuscript/MAIN_TEXT.md -o docs/manuscript/MAIN_TEXT.pdf \
  --citeproc --bibliography=docs/manuscript/references.bib \
  --pdf-engine=xelatex -H docs/manuscript/_pdf_header.tex \
  --resource-path=docs/manuscript:.
```

## Repository layout

```
analyses/         numbered analysis and figure scripts (each has a header docstring)
src/              the mls_emergence package (transmission models, signatures, I/O)
data/             source ceramic tables, radiocarbon, and geospatial reference layers
docs/manuscript/  the manuscript (MAIN_TEXT.md, references.bib, built docx/PDF)
figures/          manuscript figures (committed); run_all.sh regenerates them
output/           generated text/markdown results       (created by run_all.sh)
```

## Installation

Python 3.11 or newer. The geospatial stack needs system GEOS and PROJ:

```bash
# macOS
brew install geos proj
# Debian/Ubuntu
sudo apt-get install libgeos-dev libproj-dev proj-data proj-bin
```

Then, from the repository root:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
pip install -e .          # makes the `mls_emergence` package importable
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

The full run takes a while: several scripts run Monte Carlo simulations with
many seeds. The run continues past any single failure and reports the count at
the end.

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
| `make_figures.py` | House-style pipeline for most main-text and supplemental figures |
| `make_map.py` | Shared map helpers |

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

> Lipo, Carl P., and Robert J. DiNapoli. Are the phases real? Distinguishing
> bounded interaction groups from spatially structured drift in central
> Mississippi Valley decorated ceramics.

## Authors

Carl P. Lipo and Robert J. DiNapoli.
