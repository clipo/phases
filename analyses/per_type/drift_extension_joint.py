"""Joint tabulations of the drift-extension cell table (exploratory; run 2026-11-14 without the Codex build audit).

Reads output/findings/drift_extension_cells.csv, written by drift_extension.py,
and counts, per operating point, the cells at which sets of the record's
summaries each lie inside their own 10th to 90th percentile band:

  level        the three common-class F values (the set L of the reading rules)
  rare         how many of the seven rare classes have their record occupancy inside its band
  space3       D_dist, A_lat of class 0 and A_lat of class 2 (the two classes with a latitude trend in the record)

These are counts of marginal coverage over dependent grid cells (cells that
share a window share simulated histories). They are descriptive (rule 18) and
are not the registered reading rules R1 to R5, which use all four spatial
readouts; the tabulations here were chosen after the run was seen and are
labeled exploratory for that reason.

Usage: .venv-np254/bin/python analyses/per_type/drift_extension_joint.py
Output: output/findings/drift_extension_joint.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analyses" / "per_type"))
import drift_extension_readings as R

CELLS = ROOT / "output" / "findings" / "drift_extension_cells.csv"
OUT = ROOT / "output" / "findings" / "drift_extension_joint.md"
SPACE3 = ("D_dist", "A_lat_0", "A_lat_2")
STRICT, HALF, NONE = ("lat100N", "lat100S"), ("lat50N", "lat50S"), ("c0",)


def joint(table: pd.DataFrame) -> dict:
    """Joint counts for one operating point's status table (504 cells; checked by the reading rules)."""
    out = R.read(table)                                         # validates the grid and gives L
    w = table[table["stat"].isin(R.REQUIRED)].pivot(index=R.CELL, columns="stat", values="status")
    level = (w[list(R.F3)] == "inside").all(axis=1)
    if int(level.sum()) != len(out["L"]):
        raise RuntimeError("the level set here differs from the reading rules' L")
    n_f = (w[list(R.F3)] == "inside").sum(axis=1)
    rare = (w[list(R.RARE)] == "inside").sum(axis=1)
    rare_above = (w[list(R.RARE)] == "above").sum(axis=1)
    space3 = (w[list(SPACE3)] == "inside").all(axis=1)
    coupling = w.index.get_level_values("coupling"); fw = w.index.get_level_values("F_W")
    by = lambda names: coupling.isin(names)
    return {"n_cells": len(w), "n_evaluable": len(out["E_L"]), "n_level": int(level.sum()),
            "level_amplitudes": sorted(set(fw[level])), "level_rare_max": int(rare[level].max()) if level.any() else None,
            "level_rare_counts": rare[level].value_counts().sort_index().to_dict(),
            "level_rare_above_min": int(rare_above[level].min()) if level.any() else None,
            "n_rare5": int((rare >= 5).sum()), "rare5_max_f_inside": int(n_f[rare >= 5].max()) if (rare >= 5).any() else None,
            "any_f_max_rare": int(rare[n_f >= 1].max()) if (n_f >= 1).any() else None,
            "level_by_span": level.groupby("span").sum().reindex(list(R.SPANS)).astype(int).to_dict(),
            "level_by_window": level.groupby("window").sum().reindex(list(R.WINDOWS)).astype(int).to_dict(),
            "level_and_rare5_by_span_window": int((level & (rare >= 5)).sum()),
            "a_lat2": {k: (int(((w["A_lat_2"] == "inside") & by(v)).sum()), int(by(v).sum())) for k, v in
                       (("none", NONE), ("half", HALF), ("strict", STRICT))},
            "a_lat2_above_c0": int(((w["A_lat_2"] == "above") & by(NONE)).sum()),
            "space3": {c: (int((space3 & (coupling == c)).sum()), int((coupling == c).sum())) for c in R.COUPLINGS},
            "space3_strict_fw": {f"{a:g}": int((space3 & by(STRICT) & (fw == a)).sum()) for a in sorted(set(fw))},
            "level_space3": int((level & space3).sum()),
            "level_space3_couplings": sorted(set(coupling[level & space3])),
            "level_space3_rare_max": int(rare[level & space3].max()) if (level & space3).any() else None,
            "level_space3_a_lat1": w.loc[level & space3, "A_lat_1"].value_counts().to_dict(),
            "crosstab": pd.crosstab(n_f.rename("common-class F inside"), rare.rename("rare classes inside"))}


def main() -> int:
    t = pd.read_csv(CELLS)
    L = ["# Drift extension: joint tabulations (exploratory)", "",
         ("Produced by `.venv-np254/bin/python analyses/per_type/drift_extension_joint.py` from "
          "`output/findings/drift_extension_cells.csv`. Exploratory: the grid was run on 2026-11-14 without the Codex build audit, "
          "on one seed set of 50 records per cell, and these tabulations were chosen after the run was seen. They are not the "
          "registered reading rules, whose statements are in `output/findings/drift_extension.md`. Every count is a count of grid "
          "cells at which each named summary of the record lies inside its own 10th to 90th percentile band (marginal coverage). "
          "Cells that share a window share simulated histories, so the counts are not counts of independent confirmations."), ""]
    for point in ("primary", "secondary"):
        j = joint(t[t["point"] == point])
        s2 = j["a_lat2"]
        L += [f"## {point.capitalize()} point", "",
              f"- Level: all three common-class F values inside their bands at {j['n_level']} of {j['n_evaluable']} evaluable cells "
              f"({j['n_cells']} tested), all at F_W in {{{', '.join(f'{a:g}' for a in j['level_amplitudes'])}}}. By span limit: "
              + ", ".join(f"{k} {v}" for k, v in j["level_by_span"].items()) + ". By window: "
              + ", ".join(f"{k} {v}" for k, v in j["level_by_window"].items()) + ".",
              (f"- Rare classes in those cells: number of the seven with record occupancy inside its band, as cells per count: "
              f"{j['level_rare_counts']} (maximum {j['level_rare_max']}). In every one of those cells at least "
              f"{j['level_rare_above_min']} of the seven have the record's occupancy above the band (the record has the class at more sites)."),
              (f"- Converse: {j['n_rare5']} cells have at least 5 of the 7 rare classes inside; the largest number of common-class F "
              f"values inside among them is {j['rare5_max_f_inside']}. Among cells with at least one common-class F inside, the largest "
              f"number of rare classes inside is {j['any_f_max_rare']}. Cells with level and at least 5 rare classes inside: "
              f"{j['level_and_rare5_by_span_window']}."),
              (f"- Latitude trend of the painted types (A_lat, class 2) inside its band: dates independent of location (c0) "
              f"{s2['none'][0]} of {s2['none'][1]} cells (record above the band at {j['a_lat2_above_c0']}); half coupling (lat50N, lat50S) "
              f"{s2['half'][0]} of {s2['half'][1]}; strict ordering (lat100N, lat100S) {s2['strict'][0]} of {s2['strict'][1]}."),
              "- D_dist and the latitude trends of classes 0 and 2 all inside, by coupling: "
              + ", ".join(f"{c} {a} of {b}" for c, (a, b) in j["space3"].items()) + ". Under strict ordering, by F_W: "
              + ", ".join(f"{a}: {n}" for a, n in j["space3_strict_fw"].items()) + " (of 18 each).",
              f"- Level and those three spatial readouts together: {j['level_space3']} cells"
              + (f", all under {', '.join(j['level_space3_couplings'])}; rare classes inside at most {j['level_space3_rare_max']}; "
                 f"status of the latitude trend of class 1 in them: {j['level_space3_a_lat1']}." if j["level_space3"] else "."), "",
              "Cells by number of common-class F values inside (rows) and number of rare classes inside (columns):", "",
              "| F inside \\ rare inside | " + " | ".join(str(c) for c in j["crosstab"].columns) + " |",
              "|---|" + "---|" * len(j["crosstab"].columns)]
        L += [f"| {i} | " + " | ".join(str(int(x)) for x in row) + " |" for i, row in zip(j["crosstab"].index, j["crosstab"].to_numpy())]
        L += [""]
    OUT.write_text("\n".join(L)); print("\n".join(L)); print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
