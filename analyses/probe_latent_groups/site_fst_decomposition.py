"""Decomposition of the record's site-level F_ST into per-class variance and the weighting term (2026-10-19).

ht - hs, with hs the site-weighted mean within-site Gini-Simpson and ht the
Gini-Simpson of the sherd-weighted pool, equals the sum over classes of the
site-weighted variance of the share plus the sum over classes of
(site-mean share)^2 - (pooled share)^2. The second term is the
composition-size correlation: it is zero when a site's shares are
independent of its total, as in the copying simulator. Descriptive;
procedure for the figures logged under the per-class resemblance entry.
Usage: .venv/bin/python analyses/probe_latent_groups/site_fst_decomposition.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "stan" / "basin_composition.json"
OUT = ROOT / "output" / "findings" / "site_fst_decomposition.md"


def main() -> int:
    d = json.loads(DATA.read_text()); y = np.asarray(d["y"], float); cls = d["_classes"]
    n = y.sum(1); P = y / n[:, None]; pb = P.mean(0); pooled = y.sum(0) / y.sum()
    hs = float(np.mean(1 - (P ** 2).sum(1))); ht = float(1 - (pooled ** 2).sum()); ht_site = float(1 - (pb ** 2).sum())
    var = P.var(0); wt = pb ** 2 - pooled ** 2
    assert abs((ht - hs) - (var.sum() + wt.sum())) < 1e-12
    L = ["# Site-level F_ST of the record: per-class variance and the weighting term", "",
         (f"Produced by `analyses/probe_latent_groups/site_fst_decomposition.py` from `data/stan/basin_composition.json` "
          f"({y.shape[0]} sites, {int(y.sum())} sherds). hs {hs:.4f}, ht {ht:.4f}, 1 - hs/ht = {1 - hs / ht:.4f}; with a site-mean "
          f"pool ht {ht_site:.4f} and 1 - hs/ht = {1 - hs / ht_site:.4f}. ht - hs = {ht - hs:.4f} = sum of per-class site-weighted "
          f"variances {var.sum():.4f} + sum of (site-mean share)^2 - (pooled share)^2, {wt.sum():.4f}. Descriptive."), "",
         "| class | site-mean share | pooled share | site-weighted variance | (site mean)^2 - pooled^2 |", "|---|---|---|---|---|"]
    L += [f"| {c} | {pb[k]:.4f} | {pooled[k]:.4f} | {var[k]:.4f} | {wt[k]:+.4f} |" for k, c in enumerate(cls)]
    L.append("")
    OUT.write_text("\n".join(L)); print("\n".join(L[2:4])); print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
