# Data provenance

This directory holds the data files the analyses read. Files in `raw/` are
primary data transcribed or compiled from published sources, or corrections to
them with the source of each correction in the file. Files in `processed/` are
derived from the raw data by scripts in this repository and can be
regenerated. `stan/` holds the data contracts the Stan models read. The
shapefiles and the top-level settlement tables are base geographic and
settlement data for the Lower Mississippi Valley (LMV). Full bibliographic
entries for every source cited below are in `docs/manuscript/references.bib`.
All data are released under CC BY 4.0 (`LICENSE-data`) except where a file is
noted as withheld.

## `raw/` — primary data and documented corrections

| File | Contents | Source |
|---|---|---|
| `PFGData.xlsx` (sheet `SherdData`; sheet `Sheet2` type attributes) | The definitive Phillips, Ford and Griffin survey counts: 393 assemblages by 27 ceramic types. The analysis matrix takes its counts from here for every assemblage not in Lipo's compilation. | Phillips, Ford and Griffin (1951), transcribed by Lipo (2001). |
| `PFGData_sherds.csv`, `PFGData_types.csv` | The same counts and type attributes as CSV. | As above. |
| `mainfort-pfg-cpl.xlsx` (tab `cpl-pfg`), `mainfort-pfg-cpl.csv` | Lipo's (2001) compilation: the PFG counts plus his 1996 to 1997 controlled surface collections at seven assemblages (Beck, Belle Meade, Cramor Place, Castile Landing, Nickel, Rose Mound, Holden Lake). The analysis matrix takes this tab's row for every assemblage it holds. The workbook's other sheets, including the `Mainfort` sheet and the summed `pfg-cpl-mainfort` sheet the earlier version of this study used, are kept so the comparison in `output/findings/matrix_integrity.md` can be checked; no current result reads them. | Lipo (2001); Phillips, Ford and Griffin (1951). |
| `mainfort-pfg-cplXY.txt` | Latitude and longitude for the compiled assemblages. | Georeferenced from the Lower Mississippi Survey site files. |
| `mainfort-spatial-data.xlsx` (sheet `From_Lipo`) | The PFG site table's UTM coordinates, used to check which file coordinates are PFG UTMs read in the wrong datum (`scripts/measure_datum.py`, `output/findings/datum_conversion.md`). | Phillips, Ford and Griffin (1951) site table. |
| `coordinate_corrections.csv` | Five assemblage coordinates replaced, each with its source (author-supplied coordinates checked against PFG's legal descriptions and the BLM PLSS cadastral service, or PFG's own UTM). | Stated per row. |
| `settlement_coordinate_corrections.csv` | Five settlement-table coordinates replaced by PFG's UTM where the compilation's coordinate fell outside the recorded section. | Stated per row. |
| `mound_height_corrections.csv` | Parkin's mound height and count from the published descriptions, superseding the compilation's entry. | Morse (1981, 1990), as keyed in the file. |
| `phillips1970_phase_outlines.csv` | The Nodena, Parkin, Walls and Kent phase areas of Phillips (1970) as redrawn by Lipo (2001: Figure 2.6), read from the vector artwork and georeferenced from the figure's frame by `scripts/build_phillips1970_outlines.py`, which records the check in the file header. Figure 1 and the phase assignments. | Phillips (1970); Lipo (2001). |
| `mainfort2003_table1.csv` | Decorated and plain type counts for the 39 sites of Mainfort's Table 1, transcribed from the page image and verified cell for cell by a second reader (header). The comparison set of the supplement's second case; never merged into the analysis matrix. | Mainfort (2003a: Table 1). |
| `mainfort2003_rim_table1_counts.csv` | Rim attribute counts (beaded applique strip, noded, widely spaced nodes, notched top, notched exterior, beveled, N) for the sites of Mainfort's rim study, transcribed from the page image. Supplies the beveling and notching counts of the rim analysis. | Mainfort (2003b: Table 1). |
| `14CDatesFromMainfort2001.csv`, `.xls` | Radiocarbon determinations for the basin (sample, provenience, uncorrected years BP, published calibrated range). | Mainfort (2001). |
| `intcal20.14c` | IntCal20 atmospheric calibration curve, used to recalibrate the determinations. | Reimer et al. (2020). |
| `pfg-cpl-frequency.xlsx` | Seriation-ordered decorated-class frequency matrix of the earlier version's analyses. | Lipo (2001). |
| `williams1954_cmv_counts.tsv` | Southeast-Missouri type counts by site, used by the earlier version's regional extension. | Williams (1954). |

**Withheld.** `raw/mainfort_miss_rims.xls`, Robert Mainfort's unpublished
tabulation of rim attributes by site, and `processed/mainfort_rim_counts.csv`,
the counts `scripts/build_mainfort_rim_counts.py` derives from it together
with his published Table 1, are not in this release. The current rim analysis
(`analyses/rims/rim_r1.py`) reads lip form and Bell paste from the derived
file, so it cannot be rerun here until he agrees to their release. Its
findings file, `output/findings/rim_r1.md`, and the figure drawn from it are
committed.

## How the analysis matrix is built

`scripts/build_analysis_matrix.py` writes `processed/analysis_matrix.csv`
(38 assemblages by the ten decorated classes) and
`processed/analysis_matrix_provenance.csv`, which names the source of every
row. The rule, by author ruling of 2026-09-21: **one source per assemblage,
never a sum.** An assemblage in the `cpl-pfg` tab of Lipo's compilation takes
that row; otherwise it takes its `PFGData.xlsx` row, matched on the normalized
site name. Decorated classes only; plain wares are excluded as non-stylistic.
`Barton/Kent/MPI` pools Barton Incised, Kent Incised and Mound Place Incised.

The earlier version of this study used the workbook's `pfg-cpl-mainfort`
sheet, which summed the compilation with a `Mainfort` sheet. The two are not
independent: the `Mainfort` sheet's counts are transcriptions of the same
survey files, so PFG sherds were counted twice, and that sheet's Parkin
Punctated column is doubled relative to Mainfort's published table
(`output/findings/matrix_integrity.md`, `scripts/check_matrix_against_mainfort2003.py`).
Mainfort's published counts exceed PFG's at ten sites and may include
collections of his own, which is why his table is kept as a separate
comparison (`processed/mainfort2003_matrix.csv`) rather than corrected and
merged. `CORRECTIONS.md` records the change.

The study set is the 28 assemblages of `processed/basin_members_curated.txt`
(14,101 decorated sherds; phases Parkin 11, Kent 8, Walls 9), produced by
`analyses/16_basin_membership.py`: the St. Francis basin assemblages with at
least 75 decorated sherds, the minimum justified in
`output/findings/minimum_sample_size.md` (Table S1), with the Parchman phase
excluded on geography (its three assemblages lie across the Mississippi, 30 to
46 km from the nearest other assemblage) and Nodena's single surviving
assemblage treated as unassigned.

## `processed/` — derived (regenerable)

| File | Contents | Produced by |
|---|---|---|
| `analysis_matrix.csv` | The decorated-class matrix, 38 assemblages. | `scripts/build_analysis_matrix.py` |
| `analysis_matrix_provenance.csv` | Source used for each assemblage, with the earlier summed total beside the corrected one. | `scripts/build_analysis_matrix.py` |
| `assemblage_sources.csv` | Which assemblages have a compilation share and whether each is in the basin set. | `scripts/build_analysis_matrix.py` |
| `basin_members_curated.txt`, `basin_members_broad.txt`, `basin_members_watershed.txt` | Basin membership under the curated drainage rule (the study set), the broader corridor rule and the watershed rule. | `analyses/16_basin_membership.py` |
| `mainfort2003_matrix.csv` | Mainfort's (2003a) Table 1 as a decorated-type matrix with his phase assignments, for the supplement's second case. | `scripts/build_analysis_matrix.py` |
| `mainfort2003_sites.csv` | Coordinates for the 39 sites of that table, with the matching rule in the file header. | `scripts/build_mainfort2003_sites.py` |
| `williams1954_cmv_coords.tsv` | Coordinates for the Williams (1954) assemblages (earlier version). | Georeferenced to the LMV Survey grid. |

## `stan/` — data contracts for the Stan models

`basin_composition.json` is the 28-assemblage count matrix with site
coordinates and the distance matrix as the spatial field, latent-unit and
barrier models read it. `a0_recovery_share*_rho*.json` are the four records
simulated from the spatial field model at known spatial share and range for
its recovery checks (Table S4), with their planted values in
`a0_recovery_truth.json`. Both are written by `analyses/49_write_stan_data.py`,
which records the sha256 of the counts and coordinates in each file.
`sim_rho*.json` and `sim_truth.json` are the earlier version's recovery
records, read by `analyses/02_spatial/01_recovery.R`.

## Top-level settlement tables

| File | Contents | Source |
|---|---|---|
| `LMVData.xlsx` | The regional settlement compilation: site number and name, type, period, UTM coordinates in two zone sheets, mound count, height and area, fortification. The settlement analyses read this file with the corrections in `raw/`. | Lipo and Dunnell (2007) compilation. |
| `LMVData_locations.csv`, `LMVData-22March2006.csv`, `LMVData-22March2006.xls`, `ArcGIS Files/LMVData-22March2006.xls` | Earlier exports of the same compilation read by the earlier version's scripts and by the map helpers. | As above. |
| `Alvey LMV 14C Database_220714.xlsx` (also under `14C_data/`) | A wider radiocarbon compilation for the valley, used to ask whether adding its dates moves the basin chronology (`analyses/45_alvey_14c_robustness.py`). | Alvey's LMV radiocarbon database, file dated 14 July 2022. |

## `Shapefiles/` and `HydroRIVERS_v10_na.gdb.zip` — geographic layers

Project-compiled layers for the Lower Mississippi Valley, each a standard ESRI
shapefile set: `LMVHydrology` (stream centerlines), `LMVMajorRivers`,
`LMVcounties4`, `LMVgeology3`, `LMVstates`. `hydrorivers_lmv_cmv.gpkg` and
`hydrorivers_basin_all_orders.gpkg` are clips of HydroRIVERS v1.0 (Lehner and
Grill 2013, North America layer) made by `scripts/clip_hydrorivers.py` from
`HydroRIVERS_v10_na.gdb.zip`; the second carries every stream order for the
basin and is what the settlement-and-waterways analysis (Table S3) reads.
HydroRIVERS is distributed under the HydroSHEDS license (free for
non-commercial and commercial use with attribution);
https://www.hydrosheds.org/products/hydrorivers. The `LMV*` layers stop near
36.6 degrees N; the figure scripts extend state boundaries with Natural Earth
data, which cartopy downloads on demand.

## `CMV/` — central Mississippi Valley point layers

Point shapefiles (`GEOFILE`, `PFGOUT`, `PHILLIPS`) carrying assemblage
locations on the Lower Mississippi Survey grid, used in georeferencing and by
the earlier version's regional maps.

## Notes

- The decorated-ceramic record is time-averaged surface material, and the
  inference is object-mediated (counts of sherds, not potters); the manuscript
  states both limits.
- Coordinates are site locations georeferenced from published maps and survey
  tables, not survey-grade fixes. Where a file coordinate was demonstrably a
  PFG UTM read in the wrong datum, it is read as NAD27
  (`output/findings/datum_conversion.md`); the measured effect on every
  result is null at basin scale. Please treat site locations with the care
  archaeological site data warrants.
