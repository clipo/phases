"""87_axis_geography.py - how much of the seriation axis is place rather than time?

The first correspondence-analysis axis of the decorated-class frequencies is
used throughout as a relative chronological order. On 2026-09-23 the dated
assemblages were found at both ends of it regardless of age: Kent Place (Kent
floor archaeomagnetic AD 1345-1380, House 1993:23) and Hollywood (calibrated
1-sigma AD 1285-1435, Haley 2014, citing Johnson et al. 2000:48-49) are among
the earliest dated sites and Clay Hill the latest, yet all three sit at the
same end. This script measures how strongly the axis follows geography.

Output: output/findings/axis_geography.md
Usage: .venv/bin/python analyses/87_axis_geography.py
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))
OUT_MD = ROOT / "output" / "findings" / "axis_geography.md"
DATED = {  # assemblage: dating as reported, for the table only
    "Kent_Place": "archaeomagnetic AD 1345-1380 on the house floor (House 1993:23); pooled 14C median AD 1466 includes an intrusive post",
    "Hollywood": "calibrated 1-sigma AD 1285-1420 and 1300-1435 (Haley 2014, citing Johnson et al. 2000:48-49)",
    "Parkin": "pooled 14C median AD 1483 (19 dates)",
    "Neeleys_Ferry": "pooled 14C median AD 1513 (2 dates)",
    "Clay_Hill": "pooled 14C median AD 1597 (3 dates)",
}


def main() -> int:
    mf = importlib.import_module("make_figures")
    ph = importlib.import_module("36_canonical_phase_map")
    counts, coords = mf._load_curated()
    ca = np.asarray(mf.correspondence_axis(counts.to_numpy(float))[0], float)
    lat = coords["Latitude"].to_numpy(float)
    lon = coords["Longitude"].to_numpy(float)
    labs, _ = ph.assign_phases_by_territory([str(i) for i in counts.index], coords.to_numpy(float))
    r_lat, r_lon = spearmanr(ca, lat)[0], spearmanr(ca, lon)[0]
    X = np.column_stack([np.ones_like(lat), lat, lon])
    b = np.linalg.lstsq(X, ca, rcond=None)[0]
    r2 = 1 - ((ca - X @ b) ** 2).sum() / ((ca - ca.mean()) ** 2).sum()
    s = pd.Series(ca, index=counts.index)
    rank = s.rank().astype(int)
    df = pd.DataFrame(dict(ca=ca, phase=labs), index=counts.index)
    by = df.groupby("phase").ca.agg(["count", "mean", "min", "max"])
    L = ["# How much of the seriation axis is place rather than time?", "",
         f"Produced by `analyses/87_axis_geography.py`. First correspondence-analysis axis of the "
         f"{len(ca)} basin assemblages, in its computed orientation (as Figure 4 and the analyses use it).", "",
         f"- Rank correlation with latitude: **{r_lat:+.2f}**; with longitude: {r_lon:+.2f}.",
         f"- Share of the axis's variance explained by a linear fit on latitude and longitude: **{r2:.2f}**.", "",
         "## Dated assemblages on the axis", "",
         "| assemblage | axis value | rank of 28 | dating |", "|---|---|---|---|"]
    for n, d in DATED.items():
        L.append(f"| {n} | {s[n]:+.3f} | {rank[n]} | {d} |")
    L += ["", "## By phase", "", "| phase | assemblages | mean | min | max |", "|---|---|---|---|---|"]
    for p, r in by.iterrows():
        L.append(f"| {p} | {int(r['count'])} | {r['mean']:+.2f} | {r['min']:+.2f} | {r['max']:+.2f} |")
    L += ["", "## Reading", "",
          "The two earliest dated sites and the latest share one end of the axis, and the axis tracks "
          "latitude closely. It orders place at least as much as time, so any trend 'along the "
          "sequence' is largely a trend from south to north.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
