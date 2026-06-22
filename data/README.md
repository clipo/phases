# Data provenance

This directory holds every data file the analyses read. Files in `raw/` are
primary data transcribed or compiled from published sources; files in
`processed/` are derived from the raw data by scripts in this repository and can
be regenerated. Shapefiles and the top-level settlement tables are base
geographic and settlement data for the Lower Mississippi Valley (LMV). Full
bibliographic entries for every source cited below are in
`docs/manuscript/references.bib`. All data are released under CC BY 4.0
(`LICENSE-data`).

## `raw/` — primary data

| File | Contents | Source |
|---|---|---|
| `PFGData_sherds.csv` | Assemblage-by-type sherd counts (site name, site number, 27 ceramic types) for the Lower Mississippi Survey collections, the 266-assemblage source matrix. | Phillips, Ford & Griffin (1951), transcribed into machine-readable form by Lipo (2001). |
| `PFGData_types.csv` | Attribute table for the culture-historical types (type name, temper, surface treatment, decoration). | Phillips, Ford & Griffin (1951). |
| `mainfort-pfg-cpl.csv` | Curated decorated-class count matrix (assemblage x 10 collapsed decorated classes) used for the basin and valley transmission analyses; assemblages keyed by LMV Survey trinomial/grid designation. | Decorated counts from Phillips, Ford & Griffin (1951) / Lipo (2001); assemblage set follows the phase assignments of Mainfort (1996b). |
| `mainfort-pfg-cplXY.txt` | Latitude/longitude for the curated assemblages. | Georeferenced from the Lower Mississippi Survey site files. |
| `pfg-cpl-frequency.csv` | Seriation-ordered decorated-class frequency matrix (seriation number, assemblage, classes). | Derived from the PFG counts; ordering from frequency seriation (Lipo 2001). |
| `williams1954_cmv_counts.tsv` | Southeast-Missouri (central-valley) decorated and plain type counts by site, with physiographic region, printed page, type, count, percent, and stated total. | Transcribed from the tables in Williams (1954). |
| `14CDatesFromMainfort2001.csv` | Radiocarbon determinations for the central valley (sample ID, provenience, uncorrected years BP, published 1-sigma calibrated range). | Compiled by Mainfort (2001). |
| `intcal20.14c` | IntCal20 atmospheric radiocarbon calibration curve. | Reimer et al. (2020); distributed with the curve, used to recalibrate the Mainfort dates. |

## `processed/` — derived (regenerable)

| File | Contents | Produced by |
|---|---|---|
| `basin_members_curated.txt` | St. Francis basin assemblage membership under the curated drainage rule (within 20 km of the St. Francis / Tyronza / L'Anguille system). | `analyses/16_basin_membership.py` |
| `basin_members_broad.txt` | Basin membership under the broader corridor rule. | `analyses/16_basin_membership.py` |
| `williams1954_cmv_coords.tsv` | Coordinates (lat/lon) for the Williams (1954) assemblages, with the match method (grid vs. name). | Georeferenced to the LMV Survey grid from `williams1954_cmv_counts.tsv`. |

## Top-level settlement tables

| File | Contents | Source |
|---|---|---|
| `LMVData_locations.csv` | Site gazetteer and settlement attributes (number, name, area class, site type, period, UTM northing/easting, mound count/height/area, fortification). | Regional geographic/settlement database of Lipo & Dunnell (2007). |
| `LMVData-22March2006.csv` | Binary settlement features (mound, defensive ditch, St-Francis-type fortification, platform mound, mound counts and heights, UTM coordinates); a 22 March 2006 snapshot of the same database. | Regional geographic/settlement database of Lipo & Dunnell (2007). |

These provide the rank-size mound distribution and fortification proportions for
the settlement cross-check (main-text Figure 7).

## `Shapefiles/` — LMV base GIS layers

Project-compiled GIS layers for the Lower Mississippi Valley, used for the
study-area maps (Figures 1, 9, and S4) and, for the hydrography, the
river-network distance kernel in the drift simulations. Each layer is a standard
ESRI shapefile set (`.shp/.shx/.dbf/.prj` plus spatial-index siblings).

| Layer | Contents |
|---|---|
| `LMVHydrology.*` | River and stream centerlines (the basin hydrography graph is built from this layer). |
| `LMVMajorRivers.*` | Major river polylines/polygons (St. Francis, Tyronza, Mississippi). |
| `LMVcounties4.*` | County boundaries. |
| `LMVgeology3.*` | Surface geology polygons. |
| `LMVstates.*` | State boundaries. |
| `hydrorivers_lmv_cmv.gpkg` | Detailed river network across the full two-region extent, for the regional map (Figure 10). |

The `LMV*` local layers cover the western (St. Francis) side and stop near 36.6
degrees N; the basin-scale figure scripts extend the eastern and northern state
boundaries with Natural Earth data, which the scripts download on demand.

`hydrorivers_lmv_cmv.gpkg` is a clip of HydroRIVERS v1.0 (Lehner and Grill 2013,
North America layer) to the Figure 10 extent (longitude -92.4 to -87.8, latitude
33.7 to 37.9), reprojected on the fly to UTM 15N. It gives uniform river detail
across Arkansas, Missouri, Tennessee, and Mississippi for the regional map, with
the `ORD_FLOW` field driving the line width. HydroRIVERS is distributed under the
HydroSHEDS license (free for non-commercial and commercial use with attribution);
source: https://www.hydrosheds.org/products/hydrorivers.

## `CMV/` — central Mississippi Valley point layers

Point shapefiles used in georeferencing and mapping the central-valley and PFG
assemblages: `GEOFILE.*`, `PFGOUT.*`, and `PHILLIPS.*`. These carry assemblage
locations on the Lower Mississippi Survey grid.

## Notes

- The decorated-ceramic record is **time-averaged** surface material, and the
  inference is **object-mediated** (counts of sherds, not potters); the
  manuscript Methods bound both constraints.
- The chronology is **seriation-derived**; the radiocarbon files exist to break
  that circularity with an independent (if coarse) calendar anchor.
- Coordinates are approximate site locations georeferenced from published maps
  and survey tables (USGS topographic sheets and the Lower Mississippi Survey
  grid), not precise survey-grade fixes. They are included because the analysis
  is geographic; please treat site locations with the care archaeological site
  data warrants.
