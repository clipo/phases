"""export_phasecheck_example.py - write the paper's basin record in phasecheck's input format.

The worked example for `phasecheck`, and the record its validation test runs on:
the 28 St. Francis basin assemblages of the primary analysis, with the phase
each takes from Phillips's (1970) areas, as one table, and the river-distance
matrix the paper's boundary excess uses, as a second.

Outputs: examples/phasecheck/st_francis_basin.csv
         examples/phasecheck/st_francis_river_km.csv
Usage: PYTHONPATH=src .venv/bin/python scripts/export_phasecheck_example.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))
OUT_TABLE = ROOT / "examples" / "phasecheck" / "st_francis_basin.csv"
OUT_DIST = ROOT / "examples" / "phasecheck" / "st_francis_river_km.csv"


def main() -> int:
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    rev = importlib.import_module("47_revision_analysis")
    counts, coords = mf._load_curated()
    names = [str(i) for i in counts.index]
    xy = coords.to_numpy(float)
    labels, _ = ph.assign_primary_phases(names, xy)
    data = rev.load_sets()["basin"]
    if [str(n) for n in data["names"]] != names:
        raise RuntimeError("assemblage order differs between loaders")
    table = pd.DataFrame({"name": names, "latitude": xy[:, 0], "longitude": xy[:, 1], "phase": labels})
    table = pd.concat([table, counts.reset_index(drop=True).astype(int)], axis=1)
    # Position along the seriation order the paper's model samples at (0 earliest, 1 latest).
    table.insert(4, "seriation_order", data["ranks"])
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT_TABLE, index=False)
    pd.DataFrame(data["d"], index=names, columns=names).to_csv(OUT_DIST, index_label="name", float_format="%.17g")
    print(f"{len(names)} assemblages, {counts.shape[1]} classes -> {OUT_TABLE.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
