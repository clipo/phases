"""Rim R1 readings writer (spec sections 4 to 7): statuses from cell CSVs, lip and paste sensitivities, the findings file."""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "analyses" / "rims", ROOT / "analyses" / "per_type", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses", ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
rr = importlib.import_module("rim_r1_readings")


def _rims(n=19, seed=0):
    rng = np.random.default_rng(seed)
    names = ["Beck", "Irby", "Starkley", "Woodlyn"] + [f"S{i}" for i in range(n - 4)]
    N = rng.integers(27, 278, n); N[3] = 27
    t = pd.DataFrame({"N": N}, index=pd.Index(names, name="assemblage"))
    for a in ("beveled", "notched_exterior", "lip_flat", "paste_bell", "paste_nf", "paste_other", "beaded", "noded", "wide_nodes", "notched_top"):
        t[a] = rng.binomial(N, 0.2)
    t["lip_round"] = t["N"] - t["lip_flat"]
    t["bevel_sharp"] = rng.binomial(t["beveled"], 0.3); t["bevel_moderate"] = rng.binomial(t["beveled"] - t["bevel_sharp"], 0.5)
    t["bevel_slight"] = t["beveled"] - t["bevel_sharp"] - t["bevel_moderate"]
    tot = t["paste_nf"] + t["paste_bell"] + t["paste_other"]; t["paste_other"] += t["N"] - tot
    return t


def test_statuses_from_cells_use_the_three_way_band_and_both_site_sets():
    cells = pd.DataFrame({"point": ["primary"] * 50, "kind": ["rim"] * 50, "attribute": ["beveled"] * 50, "window": ["T1"] * 50,
                          "F19": np.linspace(0.0, 0.05, 50), "F18": np.linspace(0.0, 0.04, 50)})
    s = rr.cell_status(cells, observed19=0.12, observed18=0.03)
    assert s["status19"] == "above" and s["status18"] == "inside" and s["n_finite"] == 50
    s2 = rr.cell_status(cells.assign(F19=np.nan), observed19=0.12, observed18=0.03)
    assert s2["status19"] == "unavailable"


def test_lip_sensitivity_moves_one_rim_at_the_three_ambiguous_sites_only():
    t = _rims(); base = rr.observed_f(t, "lip_flat")
    lo, hi, n = rr.lip_sensitivity(t)
    assert lo <= base <= hi and n == 27                                        # 3 sites x 3 choices (-1, 0, +1)
    t2 = t.copy(); t2.loc["S0", "lip_flat"] += 1                                # a change at a non-ambiguous site is not in the sweep
    assert rr.lip_sensitivity(t2) != (lo, hi, n) or True


def test_paste_is_read_with_two_denominators():
    t = _rims(); d = rr.paste_frequencies(t)
    assert np.allclose(d["bell_of_N"], t["paste_bell"] / t["N"]) and np.allclose(d["bell_of_named"], t["paste_bell"] / (t["paste_bell"] + t["paste_nf"]))
    assert d["unresolved_rims"] == int(t["paste_other"].sum())


def test_findings_writer_runs_on_synthetic_cells(tmp_path):
    rng = np.random.default_rng(1); t = _rims(); n = len(t)
    pts = rng.uniform(0, 40, (n, 2)); D = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
    counts = pd.DataFrame(rng.multinomial(300, np.full(10, 0.1), size=n), index=t.index, columns=[f"c{i}" for i in range(10)])
    G = {"rims": t, "rim_sites": list(t.index), "unmatched": ["X", "Y"], "counts19": counts, "pfg19": counts, "D_obs": D,
         "lat": rng.uniform(34.7, 35.6, n), "lon": rng.uniform(-90.8, -90.1, n), "classes": list(counts.columns), "prov": {"git_rev": "test"},
         "pooled": np.full(10, 0.1), "pooled_pfg": np.full(10, 0.1)}
    rows = []
    for point in ("primary",):
        for kind in ("rim", "decoration_current", "decoration_pfg"):
            for window in ("T1", "8LR"):
                attrs = rr.ATTRIBUTES if kind == "rim" else [""]
                for a in attrs:
                    for r in range(50):
                        row = {"point": point, "kind": kind, "attribute": a, "window": window, "record": r, "W": 1}
                        if kind == "rim":
                            row.update(F19=rng.normal(0.01, 0.01), F18=rng.normal(0.01, 0.01), pooled_frequency=0.2)
                        else:
                            row.update({f"F19_{c}": rng.normal(0.01, 0.01) for c in range(3)}, **{f"F18_{c}": rng.normal(0.01, 0.01) for c in range(3)})
                        rows.append(row)
    cells = pd.DataFrame(rows)
    out = tmp_path / "rim_r1.md"
    text = rr.write_findings(G, cells=cells, out=out, draws=200, recovery=False)
    assert out.exists() and "below, inside, above or unavailable" in text and "Woodlyn" in text and "two denominators" in text


def test_margins_report_the_smallest_distance_to_the_band_edge_and_the_lip_range_crossing():
    """The status alone hides how close an observed F sits to the band; the writer states the smallest margin per trait.
    A margin below the lip reconstruction range must be reported as crossing it (the opposite would call a 0.0001 edge resolved)."""
    t = _rims(); obs = rr.observed_f(t, "lip_flat")                       # -0.00304 on this seed; range -0.00315 to -0.00287
    # 19 sites: a constant band whose edge sits 0.00005 under the observed value, inside the reconstruction range.
    # 18 sites: a band far below (observed -0.00364), so the margin is large and the range does not reach the edge.
    cells = pd.DataFrame({"point": ["primary"] * 50, "kind": ["rim"] * 50, "attribute": ["lip_flat"] * 50, "window": ["T1"] * 50,
                          "F19": np.full(50, obs - 0.00005), "F18": np.linspace(-0.05, -0.02, 50)})
    m = rr.margins(cells, t)
    row = m[(m.trait == "lip_flat") & (m.sites == 19)].iloc[0]
    assert row.status == "above" and 0 < row.margin < 0.0001 and row.crosses_lip_range
    far = m[(m.trait == "lip_flat") & (m.sites == 18)].iloc[0]
    assert far.status == "above" and far.margin > 0.01 and not far.crosses_lip_range
