"""Does the composition-size pattern, or the common classes' departure, depend on which collection the sherds came from? (D36 closing step)

The analysis matrix takes each of the 28 basin assemblages from one source
(`scripts/build_analysis_matrix.py`, author's ruling of 2026-09-21): 18 are
the 1940 Phillips-Ford-Griffin surface collections alone (PFGData.xlsx); 9
are those collections plus Lipo's 1996 to 1997 collections at the same site,
as the cpl-pfg tab of mainfort-pfg-cpl.xlsx sums them: six large-count
additions (Belle Meade +1,455 decorated sherds, Beck +1,154, Nickel +1,120,
Rose Mound +738, Cramor Place +286, Castile Landing +220) and three small
ones (Hollywood +25, of which 24 Painted; Barton Ranch +9; Neeleys Ferry
+4); Holden Lake is Lipo's collection alone. The six large additions are at
six of the eight largest collections, so collection size and collection
source are confounded. Both parts of each mixed row exist in the raw files,
so the nine sites give paired comparisons of two collections of the same
site by different collectors fifty years apart.

Readings, descriptive throughout (rule 18):

1. The per-class Spearman correlation of share with collection size, over
   all 28 rows (the record's +0.42 for Parkin_Punctated and -0.57 for
   Painted), within the 18 rows with no addition, within the 21 without a
   large-count addition, and, matched on the same 27 sites, on the current
   counts and on the PFG parts alone (Codex CC-4).
2. At each of the nine mixed sites and pooled, per class: counts and shares
   in the PFG part and in the Lipo addition, and the log odds ratio Lipo
   minus PFG (that order). The interval is the 5th to 95th percentile of the
   contrast under independent multinomial sampling of each part with a
   symmetric Dirichlet(alpha) prior over the ten classes, alpha = 1/2 (the
   ten-class Jeffreys prior; marginally Beta(y_c + 1/2, N - y_c + 9/2)),
   with the same medians under alpha 0.1 and 1 beside it (rule 20b). The
   pooled contrast is between two collection-size-weighted mixtures and is
   labeled so; the per-site rows are the paired reading (CC-4). The intervals
   describe uncertainty in the latent class proportions of each collection
   under that model only; they do not include source attribution, vessel
   clustering, unequal sampling of site areas, recovery or classification
   differences, or occupation-time differences between the two visits.
3. The per-class moment F and occupancy of the per-type table recomputed on
   the 18 rows with no addition, and matched on the 27 sites with a PFG part
   of at least 75 decorated sherds, on the current counts and on the PFG
   parts alone (Holden Lake has no PFG part and is outside the 27).

Rule 20 for the prior: it leans toward equal class shares (1/K) in each
part and moves a small collection further than a large one. Its effect on
the contrast between the two parts is not a shrinkage toward zero: pulling
one part toward 1/K can attenuate, amplify or reverse a contrast, and the
per-site alpha columns show all three (Neeleys Ferry's Parkin_Punctated
contrast is +1.43, -0.41 and -1.02 under alpha 0.1, 0.5 and 1, from 4 of 4
added sherds against 760 of 1,262); the prior predictive gives
each class a median share of about 0.05 with a 90 percent interval from
under 0.001 to about 0.37; sensitivity is the alpha table; recovery is tested on constructed
counts with a known contrast (tests/per_type/test_composition_by_collection.py).

The decomposition is side-effect free (CC-1): it reads the cpl-pfg tab and
the PFG site list of mainfort-pfg-cpl.xlsx directly rather than running the
reconciliation script, whose matching it reproduces (the name-to-grid alias
table is copied from it and the agreement is tested in a sandbox). PFGData's
composite site numbers are summed by base number with the one ruled
exception: Kent Place's "13-N-4/C" repeats collection C of "13-N-4/B,C,D,E"
(build script docstring), so only the latter is its PFG part. Guards on
every run: each Lipo-tab assemblage maps to exactly one tab id; the raw PFG
rows chosen for one assemblage carry disjoint collection letters and no raw
row serves two assemblages; only Holden Lake may lack a PFG part; every
addition is non-negative (the addition is defined as matrix row minus PFG
part, so that identity holds by construction and is not a check). The source
rows used are written into the findings as a provenance table.
Usage: .venv/bin/python analyses/per_type/composition_by_collection.py
Output: output/findings/composition_by_collection.md
"""
from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses" / "per_type", ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from mls_emergence.dataio.basin import load_curated

OUT = ROOT / "output" / "findings" / "composition_by_collection.md"
WORKBOOK = ROOT / "data" / "raw" / "mainfort-pfg-cpl.xlsx"
MIN_SHERDS = 75
COMMON = 3
SEED = 2026101960
DRAWS = 20000
ALPHA = 0.5
ALPHAS = (0.1, 0.5, 1.0)
LARGE = 100          # decorated sherds: the six ruled additions are 220 or more, the three small ones 25 or fewer; a count-size label only
NO_PFG_ALLOWED = ("Holden_Lake",)
# Name-to-grid aliases of scripts/reconcile_ceramic_matrix.py (C_GRID), copied so that no script with side effects is run; the
# agreement is tested in tests/per_type/test_composition_by_collection.py.
C_GRID = {"kentplace": "13-N-4", "lakecormorant": "13-P-8", "castilelanding": "13-N-21", "nickel": "13-N-15",
          "rosemound": "12-N-3", "holdenlake": "Holden_Lake"}
# Composite PFGData rows that repeat a collection already inside another row of the same base number (build script docstring).
PFG_ROW_OVERRIDE = {"13-N-4": "13-N-4/B,C,D,E"}


def norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def tab_ids(name: str, tab_ids_all: list[str], g2n: dict) -> list[str]:
    """The cpl-pfg tab ids an assemblage name matches, as the reconciliation script matches them."""
    k = norm(name)
    gid = C_GRID.get(k)
    return [gid] if gid else [i for i in tab_ids_all if norm(g2n.get(i, i)) == k]


ALL_LETTERS = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZ")


def collection_letters(site_id: str) -> frozenset:
    """The collection letters of a composite PFGData site number ('13-N-4/B,C,D,E' -> {B,C,D,E}).

    A row without a suffix is the whole site, so it is taken to contain every
    component and overlaps any component row of the same base (Codex WB-7).
    """
    if "/" not in site_id:
        return ALL_LETTERS
    return frozenset(re.split(r"[,&]", site_id.split("/", 1)[1]))


def one_tab_id(name: str, tab_ids_all: list[str], g2n: dict) -> str:
    """The single cpl-pfg id an assemblage maps to; raises when the match is not exactly one id occurring exactly once in the tab."""
    ids = tab_ids(name, tab_ids_all, g2n)
    if len(ids) != 1:
        raise ValueError(f"{name}: expected one cpl-pfg id, found {ids}")
    if tab_ids_all.count(ids[0]) != 1:
        raise ValueError(f"{name}: cpl-pfg id {ids[0]} occurs {tab_ids_all.count(ids[0])} times in the tab")
    return ids[0]


def load_sources() -> tuple[pd.DataFrame, list[str], dict]:
    """The cpl-pfg tab, its ids, and the PFG site number to name map, read directly (no script is run)."""
    tab = pd.read_excel(WORKBOOK, sheet_name="cpl-pfg"); tab["id"] = tab["Assemblages"].astype(str).str.strip()
    pfg = pd.read_excel(WORKBOOK, sheet_name="PFG"); pfg["Site Number"] = pfg["Site Number"].astype(str).str.strip()
    g2n = {r["Site Number"]: str(r["Site Name"]).strip() for _, r in pfg.dropna(subset=["Site Name"]).iterrows()}
    return tab, tab["id"].tolist(), g2n


def decompose() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.Series, pd.DataFrame]:
    """Matrix rows, PFG part, Lipo addition, collection class, and the provenance table (source rows used) of each basin assemblage.

    Classes: PFG (no addition), PFG+Lipo large (addition of at least LARGE sherds), PFG+Lipo small, Lipo (no PFG part).
    """
    bam = importlib.import_module("build_analysis_matrix")
    counts, _coords = load_curated(); TY = list(counts.columns)
    P = bam.pfg_table(TY); P["base"] = P["id"].str.split("/").str[0]; P["row"] = np.arange(len(P))
    tab, ids_all, g2n = load_sources()
    tabset = tab.set_index("id")[TY].apply(pd.to_numeric, errors="coerce").fillna(0)
    pfg_part = pd.DataFrame(0.0, index=counts.index, columns=TY); cls = {}; prov = []; used_rows: dict[int, str] = {}
    for name in counts.index:
        ids = tab_ids(name, ids_all, g2n)
        if ids:                                                 # a Lipo-tab row: PFG part from PFGData by grid id
            gid = one_tab_id(name, ids_all, g2n)
            if not np.allclose(tabset.loc[gid].to_numpy(float), counts.loc[name].to_numpy(float)):
                raise ValueError(f"{name}: the matrix row is not the cpl-pfg row {gid}")
            rows = P[P["id"] == PFG_ROW_OVERRIDE[gid]] if gid in PFG_ROW_OVERRIDE else P[P["base"] == gid]
            if gid in PFG_ROW_OVERRIDE and len(rows) != 1:
                raise ValueError(f"{name}: the ruled PFG row {PFG_ROW_OVERRIDE[gid]} is not present exactly once")
            letters = [collection_letters(i) for i in rows["id"]]
            for i in range(len(letters)):
                for j in range(i + 1, len(letters)):
                    if letters[i] & letters[j]:
                        raise ValueError(f"{name}: PFG rows {rows['id'].tolist()} share a collection letter")
            for r_ix in rows["row"]:
                if r_ix in used_rows:
                    raise ValueError(f"{name}: PFG row {P.loc[r_ix, 'id']} already serves {used_rows[r_ix]}")
                used_rows[r_ix] = name
            part = rows[TY].sum().to_numpy(float)
            if part.sum() == 0 and name not in NO_PFG_ALLOWED:
                raise ValueError(f"{name}: no PFG rows found and it is not a ruled Lipo-only assemblage")
            pfg_part.loc[name] = part
            extra = (counts.loc[name].to_numpy(float) - part).sum()
            cls[name] = ("Lipo" if part.sum() == 0 else "PFG" if extra == 0 else "PFG+Lipo large" if extra >= LARGE else "PFG+Lipo small")
            prov.append({"assemblage": name, "cpl-pfg id": gid, "PFG rows": "; ".join(rows["id"]) or "none", "PFG part": int(part.sum()),
                         "addition": int(extra)})
        else:                                                   # a PFGData row: the matrix row is the PFG part
            pfg_part.loc[name] = counts.loc[name].to_numpy(float); cls[name] = "PFG"
            prov.append({"assemblage": name, "cpl-pfg id": "", "PFG rows": "matrix row (PFGData.xlsx)", "PFG part": int(counts.loc[name].sum()),
                         "addition": 0})
    addition = counts.astype(float) - pfg_part
    if (addition < -1e-9).any().any():
        bad = addition[(addition < -1e-9).any(axis=1)]
        raise ValueError(f"a PFG part exceeds its matrix row, so the decomposition is wrong:\n{bad}")
    if not np.allclose(pfg_part + addition, counts.astype(float)):
        raise ValueError("PFG part plus addition does not reproduce the matrix")
    return counts.astype(int), pfg_part.astype(int), addition.round().astype(int), pd.Series(cls, name="collection"), pd.DataFrame(prov)


def _check_counts(y: np.ndarray, what: str) -> np.ndarray:
    y = np.asarray(y, float)
    if y.size == 0 or not np.isfinite(y).all() or (y < 0).any() or not np.equal(y, np.rint(y)).all():
        raise ValueError(f"{what} must be finite non-negative integer counts")
    return y


def spearman_by_class(counts: pd.DataFrame) -> pd.DataFrame:
    """Per class: Spearman of share with collection size over the given rows, with a reason where it is undefined."""
    y = _check_counts(counts.to_numpy(), "counts")
    if y.ndim != 2 or (y.sum(1) <= 0).any():
        raise ValueError("every row needs a positive total")
    n = y.sum(1); P = y / n[:, None]; occ = (y > 0).sum(0); rows = []
    for j, cl in enumerate(counts.columns):
        reason = ("fewer than 3 rows" if len(y) < 3 else "class at fewer than 2 rows" if occ[j] < 2
                  else "constant totals" if np.ptp(n) == 0 else "constant shares" if np.ptp(P[:, j]) == 0 else "")
        rho = float(spearmanr(P[:, j], n)[0]) if not reason else float("nan")
        rows.append({"class": cl, "rho": rho, "n_rows": int(len(y)), "rows_present": int(occ[j]), "undefined": reason})
    return pd.DataFrame(rows)


def log_odds_ratio(a: np.ndarray, b: np.ndarray, rng: np.random.Generator, draws: int = DRAWS, alpha: float = ALPHA) -> pd.DataFrame:
    """Per class, log odds in a minus log odds in b (a = Lipo addition, b = PFG part): independent Dirichlet(y + alpha) posteriors.

    Positive means the class is a larger share of a than of b. Both parts must
    be non-empty count vectors of the same length with at least two classes;
    a class absent from both is reported with NaN.
    """
    a = _check_counts(a, "a"); b = _check_counts(b, "b")
    if a.shape != b.shape or a.ndim != 1 or a.shape[0] < 2 or a.sum() <= 0 or b.sum() <= 0:
        raise ValueError("a and b must be non-empty count vectors of the same length with at least two classes")
    if not (isinstance(draws, (int, np.integer)) and draws > 0) or not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("draws must be a positive integer and alpha a positive number")
    pa = rng.dirichlet(a + alpha, size=int(draws)); pb = rng.dirichlet(b + alpha, size=int(draws))
    lo = np.log(pa / (1 - pa)) - np.log(pb / (1 - pb))
    q = np.quantile(lo, [0.05, 0.5, 0.95], axis=0)
    out = pd.DataFrame({"a_count": a.astype(int), "b_count": b.astype(int), "a_share": a / a.sum(), "b_share": b / b.sum(),
                        "q05": q[0], "median": q[1], "q95": q[2]})
    out.loc[(a == 0) & (b == 0), ["q05", "median", "q95"]] = np.nan
    return out


def per_class_table(counts: pd.DataFrame) -> pd.DataFrame:
    pct = importlib.import_module("per_class_table")
    f = pct.moment_f(counts.to_numpy(float)); occ = (counts.to_numpy(float) > 0).sum(0)
    return pd.DataFrame({"class": counts.columns, "F": f, "sites": occ.astype(int), "n_rows": len(counts)})


def prior_predictive(rng: np.random.Generator, alpha: float = ALPHA, k: int = 10, draws: int = DRAWS) -> tuple[float, float, float]:
    p = rng.dirichlet(np.full(k, alpha), size=draws)[:, 0]
    return tuple(float(x) for x in np.quantile(p, [0.05, 0.5, 0.95]))


def _lor_cell(r) -> str:
    return f"{r['median']:+.2f} [{r['q05']:+.2f}, {r['q95']:+.2f}]" if np.isfinite(r["median"]) else "absent from both parts"


def main() -> int:
    counts, pfg, add, cls, prov = decompose(); TY = list(counts.columns)
    large = cls.index[cls == "PFG+Lipo large"].tolist(); small = cls.index[cls == "PFG+Lipo small"].tolist()
    mixed = sorted(large + small, key=lambda a: -add.loc[a].sum()); only = cls.index[cls == "PFG"].tolist()
    lipo_only = cls.index[cls == "Lipo"].tolist(); unlarge = only + small
    matched = [a for a in counts.index if a not in lipo_only and pfg.loc[a].sum() >= MIN_SHERDS]
    below = [a for a in counts.index if a not in lipo_only and pfg.loc[a].sum() < MIN_SHERDS]
    rng = np.random.default_rng(SEED); pp = prior_predictive(rng)
    L = ["# Composition by collection: the 1940 PFG surface collections against Lipo's later additions", "",
         (f"Produced by `analyses/per_type/composition_by_collection.py`. Of the 28 basin assemblages, {len(only)} are the 1940 "
          f"Phillips-Ford-Griffin collections alone; {len(large)} are those plus a large-count addition from Lipo's 1996 to 1997 collections ("
          + ", ".join(f"{a} +{int(add.loc[a].sum())}" for a in mixed if a in large) + f"); {len(small)} carry a small one ("
          + ", ".join(f"{a} +{int(add.loc[a].sum())}" for a in small) + f"); {', '.join(lipo_only)} is Lipo's alone. Large and small "
          f"are count-size labels at {LARGE} sherds, not a judgment of compositional weight: Hollywood's 25 added sherds are 24 Painted, the "
          "largest change in any site's Painted share among the nine. Decomposition checked on this run: one cpl-pfg id per assemblage, "
          "the matrix row equal to that id's row, disjoint collection letters within and across assemblages, every addition non-negative "
          "(the addition is defined as the matrix row minus the PFG part, so that identity is not a check). Descriptive readings under rule 18. Intervals in section 2 are 5th to 95th "
          f"percentiles of the contrast under independent multinomial sampling of each part with a symmetric Dirichlet({ALPHA}) prior over "
          f"the ten classes ({DRAWS:,} draws, seed {SEED}); they describe the latent class proportions of each collection under that model "
          "and nothing about why the two collections differ (source attribution, vessel clustering, which parts of a site were walked, "
          "recovery or classification differences, or the fifty years between the visits are all outside them). Prior (rule 20): it leans "
          "toward equal class shares (1/K) in each part and moves a small collection further than a large one; its effect on the contrast "
          "between the parts is not a shrinkage toward zero, since pulling one part toward 1/K can attenuate, amplify or reverse a contrast, "
          "and the per-site alpha columns show all three (Neeleys Ferry's Parkin_Punctated contrast changes sign across the alphas); "
          f"prior predictive share per class median {pp[1]:.3f}, 90 percent interval "
          f"[{pp[0]:.4f}, {pp[2]:.3f}] (a distribution of latent shares, not of counts at these collection sizes); sensitivity under alpha "
          f"{ALPHAS[0]} and {ALPHAS[2]} is shown beside every per-site contrast; recovery on constructed counts, including ten classes at 25 "
          "against 200 sherds, is in the tests."), ""]
    # 1. share against size
    S = {"all 28 rows": spearman_by_class(counts), f"{len(only)} rows, no addition": spearman_by_class(counts.loc[only]),
         f"{len(unlarge)} rows, no large addition": spearman_by_class(counts.loc[unlarge]),
         f"{len(matched)} matched sites, current counts": spearman_by_class(counts.loc[matched]),
         f"{len(matched)} matched sites, PFG parts": spearman_by_class(pfg.loc[matched])}
    L += ["## 1. Share against collection size", "",
          ("Spearman correlation of each class's site share with the site's decorated-sherd total. The last two columns hold the site "
           f"set fixed (the {len(matched)} sites whose PFG part has at least {MIN_SHERDS} sherds) and change only the collection: the "
           "current counts against the PFG parts alone. An undefined entry says why."), "",
          "| class | " + " | ".join(S) + " |", "|---|" + "---|" * len(S)]
    for j, cl in enumerate(TY):
        cells = [(f"{t.rho.iloc[j]:+.2f}" if np.isfinite(t.rho.iloc[j]) else f"undefined ({t.undefined.iloc[j]})") + f" ({t.rows_present.iloc[j]})" for t in S.values()]
        L.append(f"| {cl} | " + " | ".join(cells) + " |")
    L += ["", "Rows present in parentheses. " + f"Totals of the {len(only)} rows with no addition: {int(counts.loc[only].sum(1).min())} to "
          f"{int(counts.loc[only].sum(1).max())} sherds; of the six large-addition rows: {int(counts.loc[large].sum(1).min())} to "
          f"{int(counts.loc[large].sum(1).max())}; PFG parts of the matched sites: {int(pfg.loc[matched].sum(1).min())} to {int(pfg.loc[matched].sum(1).max())}.", ""]
    # 2. paired comparison
    L += ["## 2. The nine mixed sites: PFG part against Lipo addition", "",
          ("Log odds ratio is the Lipo addition minus the PFG part, so a positive value means the class is a larger share of Lipo's "
           "collection than of the 1940 collection at the same site. The per-site rows are the paired reading. The pooled rows contrast "
           "two collection-size-weighted mixtures (the pooled PFG parts against the pooled additions), which can differ even when no site "
           "differs within itself; they are not a site-adjusted effect."), ""]
    for label, group in (("Pooled over the six large-addition sites", large), ("Pooled over all nine mixed sites", mixed)):
        pooled = log_odds_ratio(add.loc[group].sum(0).to_numpy(), pfg.loc[group].sum(0).to_numpy(), rng)
        L += [f"### {label} (size-weighted)", "", "| class | PFG parts: count, share | Lipo additions: count, share | log odds ratio, Lipo minus PFG, median [5%, 95%] |", "|---|---|---|---|"]
        for j, cl in enumerate(TY):
            r = pooled.iloc[j]
            L.append(f"| {cl} | {int(r['b_count'])}, {r['b_share']:.3f} | {int(r['a_count'])}, {r['a_share']:.3f} | {_lor_cell(r)} |")
        L.append("")
    L += ["### Per site, the three common classes, with the prior sensitivity", "",
          "| site | addition | PFG part total | Lipo addition total | " + " | ".join(f"{TY[j]}: PFG share, Lipo share; log odds ratio at alpha {ALPHA} [5%, 95%]; medians at alpha {ALPHAS[0]}, {ALPHAS[2]}" for j in range(COMMON)) + " |",
          "|---|---|---|---|" + "---|" * COMMON]
    for a in mixed:
        rs = {al: log_odds_ratio(add.loc[a].to_numpy(), pfg.loc[a].to_numpy(), rng, alpha=al) for al in ALPHAS}
        cells = []
        for j in range(COMMON):
            r = rs[ALPHA].iloc[j]
            cells.append(f"{r['b_share']:.3f}, {r['a_share']:.3f}; {_lor_cell(r)}; {rs[ALPHAS[0]].iloc[j]['median']:+.2f}, {rs[ALPHAS[2]].iloc[j]['median']:+.2f}")
        L.append(f"| {a} | {cls[a].split()[-1]} | {int(pfg.loc[a].sum())} | {int(add.loc[a].sum())} | " + " | ".join(cells) + " |")
    L += ["", "### Per site, all ten classes: counts in the PFG part and in the Lipo addition", "",
          "| site | " + " | ".join(TY) + " |", "|---|" + "---|" * len(TY)]
    for a in mixed:
        L.append(f"| {a} | " + " | ".join(f"{int(pfg.loc[a, c])} + {int(add.loc[a, c])}" for c in TY) + " |")
    L += ["", "Cells are PFG part + Lipo addition.", ""]
    # 3. per-class table
    T = {"28 rows, current": per_class_table(counts), f"{len(only)} rows, no addition": per_class_table(counts.loc[only]),
         f"{len(matched)} matched sites, current": per_class_table(counts.loc[matched]), f"{len(matched)} matched sites, PFG parts": per_class_table(pfg.loc[matched])}
    L += ["## 3. The per-class table on PFG material alone", "",
          (f"Moment F (common classes) and sites present (rare classes) as in the per-type table. The matched columns hold the {len(matched)} "
           f"sites with a PFG part of at least {MIN_SHERDS} sherds and change only the collection ({', '.join(lipo_only)} has no PFG part"
           + (f"; below the cutoff: {', '.join(below)}" if below else "; every PFG part clears the cutoff") + ")."), "",
          "| class | " + " | ".join(T) + " |", "|---|" + "---|" * len(T)]
    for j, cl in enumerate(TY):
        cells = [(f"F {t.F.iloc[j]:.3f}" if j < COMMON else f"{int(t.sites.iloc[j])} of {int(t.n_rows.iloc[j])} sites") for t in T.values()]
        L.append(f"| {cl} | " + " | ".join(cells) + " |")
    L += ["", "PFG parts of the nine mixed sites: " + ", ".join(f"{a} {int(pfg.loc[a].sum())}" for a in mixed) + ".", ""]
    # provenance
    L += ["## Provenance: source rows used", "", "| assemblage | cpl-pfg id | PFG rows | PFG part | addition |", "|---|---|---|---|---|"]
    for _, r in prov.iterrows():
        L.append(f"| {r['assemblage']} | {r['cpl-pfg id']} | {r['PFG rows']} | {r['PFG part']} | {r['addition']} |")
    L.append("")
    OUT.write_text("\n".join(L)); print("\n".join(L)); print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
