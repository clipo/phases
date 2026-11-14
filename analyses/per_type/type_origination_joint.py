"""Joint tabulations of the type-origination cell table (exploratory).

Reads output/findings/type_origination_cells.csv, written by
type_origination.py, and reports per operating point where the level of
variation (F of ranks 1 to 3) and the occupancy of at least five of the
seven rare ranks are both inside their bands, what richness those cells
have, and how many of their commonest types originated in the stored steps.
Counts of marginal coverage over dependent grid cells from one seed set;
descriptive (rule 18). The tabulations were chosen after the run was seen.

Usage: .venv-np254/bin/python analyses/per_type/type_origination_joint.py
Output: output/findings/type_origination_joint.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analyses" / "per_type"))
import type_origination as to

CELLS = ROOT / "output" / "findings" / "type_origination_cells.csv"
OUT = ROOT / "output" / "findings" / "type_origination_joint.md"


def joint(table: pd.DataFrame) -> dict:
    s = to.sets(table)
    med = table.pivot(index=to.CELL, columns="stat", values="median")
    status = table.pivot(index=to.CELL, columns="stat", values="status")
    lr5 = s["L"] & s["R5"]; k_in = status["K_obs"] == "inside"
    by = lambda level: {f"{k:g}" if not isinstance(k, str) else k: int(v) for k, v in lr5.groupby(level).sum().items()}
    close = lr5 & s["S3"]
    return {"n_cells": len(s), "L": int(s["L"].sum()), "LR5": int(lr5.sum()), "LR": int((s["L"] & s["R"]).sum()),
            "LR5_theta0": int(lr5[s.index.get_level_values("theta") == 0].sum()),
            "by_theta": by("theta"), "by_F_W": by("F_W"), "by_span": by("span"), "by_coupling": by("coupling"),
            "rich_range": (float(med.loc[lr5, "K_obs"].min()), float(med.loc[lr5, "K_obs"].max())),
            "beyond_range": (float(med.loc[lr5, "beyond10"].min()), float(med.loc[lr5, "beyond10"].max())),
            "k_inside_by_theta": {f"{k:g}": int(v) for k, v in k_in.groupby("theta").sum().items()},
            "k_inside_max_rare": int(s.loc[k_in, "n_rare"].max()), "k_inside_and_LR5": int((k_in & lr5).sum()),
            "all_four": int((lr5 & s["K"] & s["S3"]).sum()),
            "close": [(c, float(med.loc[c, "K_obs"]), float(med.loc[c, "beyond10"])) for c in s.index[close]],
            "new_top3": (float(med.loc[lr5, "new_in_top3"].min()), float(med.loc[lr5, "new_in_top3"].max())),
            "new_top10": (float(med.loc[lr5, "new_in_top10"].min()), float(med.loc[lr5, "new_in_top10"].max()))}


def main() -> int:
    t = pd.read_csv(CELLS)
    L = ["# Type origination: joint tabulations (exploratory)", "",
         ("Produced by `.venv-np254/bin/python analyses/per_type/type_origination_joint.py` from "
          "`output/findings/type_origination_cells.csv`. Exploratory: unaudited, one seed set of 50 records per cell, tabulations "
          "chosen after the run was seen. Every count is a count of grid cells at which each named summary of the record lies inside "
          "its own 10th to 90th percentile band. Cells that share a window, amplitude and theta share simulated histories."), ""]
    for point in ("primary", "secondary"):
        j = joint(t[t["point"] == point])
        L += [f"## {point.capitalize()} point", "",
              f"- Level (F of ranks 1 to 3 inside): {j['L']} of {j['n_cells']} cells. Level and at least five of the seven rare ranks inside: "
              f"{j['LR5']} of {j['n_cells']} (all seven: {j['LR']}); at theta = 0: {j['LR5_theta0']}.",
              f"- Those {j['LR5']} cells by theta: {j['by_theta']}; by F_W: {j['by_F_W']}; by span limit: {j['by_span']}; by date arrangement: {j['by_coupling']}.",
              f"- Richness in those cells (record 10): cell medians from {j['rich_range'][0]:g} to {j['rich_range'][1]:g} types; share of sherds beyond "
              f"rank 10 (record 0): cell medians from {j['beyond_range'][0]:.3f} to {j['beyond_range'][1]:.3f}.",
              f"- Cells with richness inside its band, by theta: {j['k_inside_by_theta']}. Largest number of rare ranks inside among them: "
              f"{j['k_inside_max_rare']}. Cells with richness inside and level and five rare ranks inside: {j['k_inside_and_LR5']}.",
              f"- Types originated in the stored steps, in the cells with level and five rare ranks: among the three commonest, cell medians "
              f"from {j['new_top3'][0]:g} to {j['new_top3'][1]:g}; among the ten commonest, from {j['new_top10'][0]:g} to {j['new_top10'][1]:g}.",
              f"- Level, five rare ranks, richness with the share beyond rank 10, and the three spatial readouts all inside: {j['all_four']} cells.",
              "- Level, five rare ranks and the three spatial readouts inside (richness aside): "
              + ("none." if not j["close"] else "; ".join(
                  f"span {c[0]}, window {c[1]}, dates {c[2]}, F_W {c[3]:g}, theta {c[4]:g} (median richness {k:g}, share beyond rank 10 {b:.3f})"
                  for c, k, b in j["close"]) + "."), ""]
    OUT.write_text("\n".join(L)); print("\n".join(L)); print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
