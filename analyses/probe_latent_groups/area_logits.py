"""Between-area logit differences and per-class logit SDs across sites (correction of the A0 adequacy sentence, 2026-10-19).

Replaces the uncommitted figures in the A0 adequacy note of 2026-10-17
(validate_a1.R; discussion log, step 6): the "-0.26 logits" there is the
Parkin-minus-Kent difference in the site mean of log(Barton/Kent/MPI over
Parkin_Punctated), not a property of Parkin_Punctated, and "1 to 4 logits"
and "1.0 to 2.9" had no procedure. Descriptive only.
Usage: .venv/bin/python analyses/probe_latent_groups/area_logits.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "analyses"))
from mls_emergence.dataio.basin import load_curated

OUT = ROOT / "output" / "findings" / "area_logits.md"


def main() -> int:
    m36 = importlib.import_module("36_canonical_phase_map")
    counts, coords = load_curated(); y = counts.to_numpy(float); cls = list(counts.columns)
    lab = np.asarray(m36.assign_primary_phases([str(i) for i in counts.index], coords.to_numpy(float))[0])
    areas = ("Parkin", "Kent", "Walls")
    pool = {a: y[lab == a].sum(0) for a in areas}
    lg = lambda t: np.log(t / (t.sum() - t))
    diff = lg(pool["Parkin"]) - lg(pool["Kent"])
    sd_logit = np.log((y + .5) / (y.sum(1, keepdims=True) - y + .5)).std(0, ddof=1)
    lr = np.log(y[:, 1] / y[:, 0])
    old = lr[lab == "Parkin"].mean() - lr[lab == "Kent"].mean()
    L = ["# Between-area logit differences and per-class logit SDs across sites", "",
         ("Produced by `analyses/probe_latent_groups/area_logits.py` from the curated 28-assemblage matrix with Phillips's "
          f"areas (`36_canonical_phase_map.assign_primary_phases`; {', '.join(f'{a} {int((lab == a).sum())}' for a in areas)}). "
          "Pooled logit share per area, Parkin minus Kent; per-class SD across sites of the logit share with a half-count "
          "correction. Descriptive; replaces the uncommitted figures of the 2026-10-17 A0 adequacy note."), "",
         "| class | Parkin-area share | Kent | Walls | pooled logit, Parkin minus Kent | SD of logit share across sites |", "|---|---|---|---|---|---|"]
    for c in range(len(cls)):
        d = "absent from the Parkin area" if not np.isfinite(diff[c]) else f"{diff[c]:+.2f}"
        L.append(f"| {cls[c]} | {pool['Parkin'][c] / pool['Parkin'].sum():.4f} | {pool['Kent'][c] / pool['Kent'].sum():.4f} | "
                 f"{pool['Walls'][c] / pool['Walls'].sum():.4f} | {d} | {sd_logit[c]:.2f} |")
    L += ["", (f"Site mean of log(Barton/Kent/MPI over Parkin_Punctated), Parkin area minus Kent area: {old:+.2f} "
               "(the figure the 2026-10-17 note labeled as the reference class's own difference)."), ""]
    OUT.write_text("\n".join(L)); print("\n".join(L)); print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
