"""100_rim_attributes.py - the attribute scale: do rim attributes show what the types do not?

EXPLORATORY. Reported in the Discussion (Limits) and in the supplement section
"Rim attributes: a first look at a finer scale". The paper's Discussion argues that every
result rests on the scale of measurement, the decorated type, and that
attributes could separate two explanations the types cannot: copying restricted
at the phase lines, which should put a step at the lines in EVERY attribute and
make the attributes break in the same places, and variation arising locally,
which need not. Rim attributes are recorded independently of decoration on many
of the same assemblages.

THE DATA. Table 1 of Mainfort (2003b), "Late Period Ceramic Rim Attribute
Variation in the Central Mississippi Valley", Southeastern Archaeology
22(1):33-46: counts of six rim attributes on 7,195 rims from 55 late-period
assemblages. By the author's ruling (2026-10-04) the PUBLISHED counts are the
record. Each attribute is a separate yes-or-no observation on a rim. The file
Mainfort sent the author (acknowledge him) also carries three things the paper
does not print: the paste of the rims, lip form, and the degree of interior
bevel. Those are reported in a separate part and labeled unpublished.
`scripts/build_mainfort_rim_counts.py` builds the count table and documents
what each column is.

THREE QUESTIONS.
  A. HOW DOES EACH ATTRIBUTE VARY ACROSS THE REGION? For every assemblage with
     at least 45 rims (the paper's own minimum): does the difference between
     two assemblages grow with the distance between them, and does the
     attribute run with latitude or longitude?
  B. DO THE ATTRIBUTES VARY TOGETHER? Rank correlations between attributes
     across assemblages. A boundary on learning would move attributes together;
     independent local variation would not.
  C. DO THE PHASE LINES DIVIDE THEM? For the assemblages of the primary
     analysis that have rim data, under the same Phillips phases: the paper's
     phase-line comparison (analysis 86) on each attribute, on the six
     attributes taken together, and on the decorated types of the same
     assemblages.

UNCERTAINTY (rule 18). Each assemblage's frequency of an attribute is drawn
from Beta(count + 1/2, N - count + 1/2), or Dirichlet for more than two
classes, and every statistic is reported as a posterior median with a 95
percent credible interval. No test.

COORDINATES. The analysis matrix's where the assemblage is in it, else the
Mainfort (2003a) site file's, else the rough coordinate in Mainfort's rim file
(18 sites, all outside the three phases of the primary analysis; they enter
parts A and B only).

Output: output/findings/rim_attributes.md
Usage: PYTHONPATH=src .venv/bin/python analyses/100_rim_attributes.py
"""
from __future__ import annotations

import importlib
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analyses"))
sys.path.insert(0, str(ROOT / "src"))

COUNTS = ROOT / "data" / "processed" / "mainfort_rim_counts.csv"
SITES_2003B = ROOT / "data" / "processed" / "mainfort2003_sites.csv"
OUT_MD = ROOT / "output" / "findings" / "rim_attributes.md"
MIN_RIMS = 45            # Mainfort (2003b): analyses include assemblages with 45 or more rims
MIN_RIMS_PHASES = 25     # part C keeps smaller assemblages of the primary analysis; counts carry the uncertainty
DRAWS = 2000
SEED = 100000

# Rim-file name -> analysis-matrix name. Matched by hand; the file's names are abbreviations.
TO_MATRIX = {
    "Woodlyn": "Woodlyn", "Belle Meade": "Belle_Meade", "Starkley": "Starkley", "Beck": "Beck",
    "Clay Hill": "Clay_Hill", "Irby": "Irby", "Walls/Harris": "Walls", "Vernon Paul": "Vernon_Paul",
    "Barton Ranch": "Barton_Ranch", "Big Eddy": "Big_Eddy", "Williamson": "Williamson",
    "NeeleysFerry": "Neeleys_Ferry", "Fortune": "Fortune", "CastileLndg": "Castile_Landing",
    "Mound Place": "Mound_Place", "Davis": "Davis", "Parkin": "Parkin", "Kent": "Kent_Place",
    "RoseMound": "Rose_Mound", "Cheatam": "Cheatham", "Notgrass": "Notgrass",
    "UpNodena": "Upper_Nodena", "Carson Lake": "Carson_Lake",
}
TO_2003B = {"Castile Landing": "Castile Lg.", "Neeleys Ferry": "Neeley's Ferry",
            "Richardson's Landing": "Richardson's Lg.", "40LA17": "40LA7"}

# attribute -> (count columns; None is the complement, N minus the named counts), class labels,
# and the class shown where one percentage stands for the attribute (None shows every class)
PUBLISHED = {
    "interior bevel": (["beveled", None], ["beveled", "not beveled"], "beveled"),
    "beaded applique strip": (["beaded", None], ["present", "absent"], "present"),
    "exterior notching": (["notched_exterior", None], ["present", "absent"], "present"),
    "notched top": (["notched_top", None], ["present", "absent"], "present"),
    "band of nodes": (["noded", None], ["present", "absent"], "present"),
    "widely spaced nodes": (["wide_nodes", None], ["present", "absent"], "present"),
}
FILE_ONLY = {
    "paste: Bell share of rims in the two pastes": (["paste_nf", "paste_bell"], ["Neeley's Ferry", "Bell"], "Bell"),
    "paste, three classes": (["paste_nf", "paste_bell", "paste_other"], ["Neeley's Ferry", "Bell", "other"], None),
    "lip form": (["lip_flat", "lip_round"], ["flat", "round"], "flat"),
    "interior bevel, by degree": (["bevel_sharp", "bevel_moderate", "bevel_slight", "bevel_none"],
                                  ["sharp", "moderate", "slight", "none"], None),
}
ATTR = {**PUBLISHED, **FILE_ONLY}


def load_rims() -> pd.DataFrame:
    if not COUNTS.exists():
        raise FileNotFoundError(f"{COUNTS} is missing; run scripts/build_mainfort_rim_counts.py")
    return pd.read_csv(COUNTS)


def counts_for(d: pd.DataFrame, attr: str) -> tuple[np.ndarray, list[str]]:
    """Counts (assemblages x classes) for one attribute."""
    cols, labs, _ = ATTR[attr]
    n = d["N"].to_numpy(float)
    named = sum(d[c].to_numpy(float) for c in cols if c is not None)
    m = np.column_stack([d[c].to_numpy(float) if c is not None else n - named for c in cols])
    if (m < 0).any():
        raise ValueError(f"negative count for {attr}: the named counts exceed N")
    return m, labs


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Rank correlation; nan when either side does not vary."""
    ra, rb = rankdata(a), rankdata(b)
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def partial_spearman(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Rank correlation of a with b after the rank of c is removed from both."""
    ra, rb, rc = (rankdata(v) - (len(v) + 1) / 2.0 for v in (a, b, c))
    if rc.std() == 0:
        return spearman(a, b)
    ea = ra - rc * (ra @ rc) / (rc @ rc)
    eb = rb - rc * (rb @ rc) / (rc @ rc)
    if ea.std() == 0 or eb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ea, eb)[0, 1])


def q(a) -> tuple[float, float, float]:
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return (float(np.median(a)), float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5)))


def best_coordinates(d: pd.DataFrame, matrix_names, matrix_xy) -> tuple[np.ndarray, list[str]]:
    """(latitude, longitude) per rim assemblage and where each came from."""
    s2 = pd.read_csv(SITES_2003B, comment="#").set_index("site")
    xy, src = [], []
    for r in d.itertuples(index=False):
        mn = TO_MATRIX.get(r.file_name)
        k2 = TO_2003B.get(r.site, r.site)
        if mn in matrix_names:
            xy.append(matrix_xy[matrix_names.index(mn)]); src.append("analysis matrix")
        elif k2 in s2.index:
            xy.append([float(s2.loc[k2, "latitude"]), float(s2.loc[k2, "longitude"])]); src.append("Mainfort 2003a site file")
        else:
            xy.append([float(r.latitude_rough), float(r.longitude_rough)]); src.append("rim file (rough)")
    xy = np.array(xy, float)
    if not ((xy[:, 0] > 34.0) & (xy[:, 0] < 37.0) & (xy[:, 1] > -91.5) & (xy[:, 1] < -89.0)).all():
        bad = [d["site"].iloc[i] for i in np.flatnonzero(~((xy[:, 0] > 34.0) & (xy[:, 0] < 37.0)))]
        raise ValueError(f"a coordinate falls outside the central Mississippi Valley: {bad}")
    return xy, src


def compare_together(tables, pidx, pts, rev, t74, n_draws, n_alt, seed, per_draw=20):
    """Analysis 86's F_ST comparison on several attributes at once: the mean
    F_ST over attributes for the phases against the same mean for same-size
    divisions around random centers and compact ones."""
    k = len(np.unique(pidx)); sizes = np.bincount(pidx); slot = np.repeat(np.arange(k), sizes)
    rng = np.random.default_rng(seed)
    alt, opt = [], []
    for _ in range(n_alt):
        a = t74.assign_exact(pts, pts[rng.choice(len(pts), size=k, replace=False)], slot)
        alt.append(a); opt.append(t74.local_search(pts, a, sizes))
    rng = np.random.default_rng(seed + 1)
    f_ph, w_alt, w_opt, m_alt, m_opt = [], [], [], [], []
    for _ in range(n_draws):
        mds = [np.array([rng.dirichlet(r + 0.5) for r in m]) * m.sum(1)[:, None] for m in tables]
        mean_f = lambda lab, mds=mds: float(np.mean([rev.fst_by(md, lab) for md in mds]))
        fp = mean_f(pidx)
        fa = np.array([mean_f(alt[i]) for i in rng.choice(n_alt, per_draw)])
        fo = np.array([mean_f(opt[i]) for i in rng.choice(n_alt, per_draw)])
        f_ph.append(fp); w_alt.append(np.mean(fp > fa)); w_opt.append(np.mean(fp > fo))
        m_alt.append(np.median(fa)); m_opt.append(np.median(fo))
    return {"phase": q(f_ph), "alt": q(m_alt), "opt": q(m_opt), "p_alt": float(np.mean(w_alt)), "p_opt": float(np.mean(w_opt))}


def main() -> int:
    mf = importlib.import_module("make_figures")
    m35 = importlib.import_module("35_basin_pullout")
    ph = importlib.import_module("36_canonical_phase_map")
    rev = importlib.import_module("47_revision_analysis")
    t74 = importlib.import_module("74_phase_partition_test")
    t86 = importlib.import_module("86_partition_posterior")

    d = load_rims()
    full_counts, full_xy = m35.load_full()
    matrix_names = [str(i) for i in full_counts.index]
    xy_all, src = best_coordinates(d, matrix_names, full_xy[["Latitude", "Longitude"]].to_numpy(float))
    d["coordinate_source"] = src
    rng = np.random.default_rng(SEED)

    # ---------------- A and B: the region
    keep = (d["N"] >= MIN_RIMS).to_numpy()
    R = d[keep].reset_index(drop=True)
    xyR = xy_all[keep]
    ptsR = t74.km_xy(xyR)
    iu = np.triu_indices(len(R), 1)
    distR = np.sqrt(((ptsR[:, None] - ptsR[None]) ** 2).sum(-1))[iu]
    n_src = pd.Series(np.array(src)[keep]).value_counts().to_dict()

    L = ["# The attribute scale: Mainfort's rim attributes (exploratory)", "",
         ("Produced by `analyses/100_rim_attributes.py`. Counts are Table 1 of Mainfort (2003b), *Southeastern "
         f"Archaeology* 22(1):33-46, as printed: {len(d)} assemblages, {int(d['N'].sum()):,} rims, six attributes, "
         "each a separate yes-or-no observation on a rim. Paste, lip form and bevel degree come from the file "
         "Mainfort sent the author and are not in the paper; they are in the last part. Every statistic is a "
         f"posterior median with a 95 percent credible interval over {DRAWS:,} draws of each assemblage's "
         "frequencies (Beta or Dirichlet, counts + 1/2)."), "",
         f"## A. How each attribute varies across the region ({len(R)} assemblages with at least {MIN_RIMS} rims)", "",
         "Coordinates: " + ", ".join(f"{v} from the {k}" for k, v in n_src.items())
         + f". Distances are straight-line; the pairs span {distR.min():.0f} to {distR.max():.0f} km.", "",
         ("| attribute | percent of rims: median assemblage (lowest to highest) | rank correlation of difference "
         "with distance | with latitude (north positive) | with longitude (east positive) |"),
         "|---|---|---|---|---|"]
    freq_draws = {}
    for attr in PUBLISHED:
        m, _labs = counts_for(R, attr)
        c, n = m[:, 0], m.sum(1)
        P = rng.beta(c + 0.5, n - c + 0.5, size=(DRAWS, len(R)))
        freq_draws[attr] = P
        rd = [spearman(np.abs(p[:, None] - p[None])[iu], distR) for p in P]
        rlat = [spearman(p, xyR[:, 0]) for p in P]
        rlon = [spearman(p, xyR[:, 1]) for p in P]
        obs = 100.0 * c / n
        f = lambda t: f"{t[0]:+.2f} ({t[1]:+.2f} to {t[2]:+.2f})"
        L.append(f"| {attr} | {np.median(obs):.1f} ({obs.min():.1f} to {obs.max():.1f}) | {f(q(rd))} | "
                 f"{f(q(rlat))} | {f(q(rlon))} |")

    L += ["", ("A positive correlation with distance means assemblages farther apart differ more in that attribute. "
          "A positive correlation with latitude means the attribute is commoner to the north."), "",
          "## B. Do the attributes vary together?", "",
          "Rank correlation between two attributes' frequencies across the same assemblages.", "",
          "| attribute pair | rank correlation (95% CrI) |", "|---|---|"]
    conc = []
    for a, b in combinations(PUBLISHED, 2):
        r = q([spearman(pa, pb) for pa, pb in zip(freq_draws[a], freq_draws[b])])
        conc.append((a, b, r))
        L.append(f"| {a} and {b} | {r[0]:+.2f} ({r[1]:+.2f} to {r[2]:+.2f}) |")
    clear = [(a, b, r) for a, b, r in conc if r[1] > 0 or r[2] < 0]

    # ---------------- C: the phase lines
    counts, coords = mf._load_curated()
    names28 = [str(i) for i in counts.index]
    data = rev.load_sets()["basin"]
    if names28 != [str(n) for n in data["names"]]:
        raise RuntimeError("assemblage order differs between loaders")
    labels28 = np.asarray(ph.assign_primary_phases(names28, coords.to_numpy(float))[0])
    d["matrix_name"] = d["file_name"].map(TO_MATRIX)
    sub = d[d["matrix_name"].isin(names28) & (d["N"] >= MIN_RIMS_PHASES)].copy()
    idx28 = np.array([names28.index(nm) for nm in sub["matrix_name"]])
    labels = labels28[idx28]
    phases = sorted(set(labels))
    pidx = np.array([phases.index(l) for l in labels])
    pts = t74.km_xy(coords.to_numpy(float)[idx28])
    dist = data["d"][np.ix_(idx28, idx28)]
    missing = [nm for nm in names28 if nm not in set(sub["matrix_name"])]

    def phase_table(attrs):
        out = ["| attribute | class | " + " | ".join(phases) + " |", "|---|---|" + "---:|" * len(phases)]
        for attr in attrs:
            m, labs = counts_for(sub, attr)
            pooled = np.array([m[pidx == g].sum(0) for g in range(len(phases))])
            pct = 100.0 * pooled / pooled.sum(1, keepdims=True)
            for lab in (labs if ATTR[attr][2] is None else [ATTR[attr][2]]):
                jj = labs.index(lab)
                out.append(f"| {attr} | {lab} | " + " | ".join(f"{pct[g, jj]:.1f}" for g in range(len(phases))) + " |")
        return out

    def line_rows(named_tables, start):
        out, res = [], {}
        for k, (attr, m) in enumerate(named_tables):
            r = t86.compare(m.astype(float), pidx, pts, rev, t74, DRAWS, 300, SEED + start + k, dist=dist)
            res[attr] = r
            out.append(f"| {attr} | {r['phase'][0]:.4f} ({r['phase'][1]:.4f} to {r['phase'][2]:.4f}) | "
                       f"{r['alt'][0]:.4f} / {r['opt'][0]:.4f} | {r['p_alt']:.2f} / {r['p_opt']:.2f} | "
                       f"{r['be_phase'][0]:+.1f} | {r['be_p_alt']:.2f} / {r['be_p_opt']:.2f} |")
            print(out[-1], flush=True)
        return out, res

    head = [("| attribute | F_ST, posterior median (95% CrI) | median F_ST, random-center / compact | "
            "P(lines separate better than random-center / compact) | boundary excess | P(larger than "
            "random-center / compact) |"), "|---|---|---|---|---:|---|"]
    L += ["", f"## C. Do the phase lines divide the attributes? ({len(sub)} assemblages of the primary analysis)", "",
          f"{len(sub)} of the {len(names28)} analyzed assemblages have at least {MIN_RIMS_PHASES} rims: "
          + ", ".join(f"{p} {int((labels == p).sum())}" for p in phases)
          + f" ({int(sub['N'].sum()):,} rims; {int(sub['N'].min())} to {int(sub['N'].max())} per assemblage, "
          f"{int((sub['N'] < MIN_RIMS).sum())} below Mainfort's minimum of {MIN_RIMS}). Without rim data: "
          + ", ".join(missing) + ".", "", "Percent of rims, pooled by phase:", ""]
    L += phase_table(PUBLISHED)
    L += ["", ("The paper's phase-line comparison (analysis 86), attribute by attribute. A probability near one "
          "half means the phase lines divide the attribute no better than arbitrary lines."), ""] + head
    dec = counts.to_numpy(float)[idx28]
    pub_tables = [(a, counts_for(sub, a)[0]) for a in PUBLISHED]
    rows_pub, res_pub = line_rows([("decorated types, same assemblages", dec)] + pub_tables, 0)
    L += rows_pub
    tog = compare_together([m for _, m in pub_tables], pidx, pts, rev, t74, DRAWS, 300, SEED + 50)
    L += [(f"| the six attributes together (mean F_ST) | {tog['phase'][0]:.4f} ({tog['phase'][1]:.4f} to "
          f"{tog['phase'][2]:.4f}) | {tog['alt'][0]:.4f} / {tog['opt'][0]:.4f} | {tog['p_alt']:.2f} / "
          f"{tog['p_opt']:.2f} | | |")]

    # ---------------- D: does any attribute go with what distance leaves unexplained?
    # The paper's open question is the difference between assemblages that
    # distance does not account for. If it were age, assemblages that differ in
    # a time-sensitive attribute (interior bevel is read as late: House 1991,
    # 1993; Mainfort 2003b) should also differ more in decoration at a given
    # distance. If it were the ware being made, it should go with paste.
    iu19 = np.triu_indices(len(sub), 1)
    d_river = dist[iu19]
    d_line = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))[iu19]
    explain = {"interior bevel (published)": counts_for(sub, "interior bevel")[0],
               "exterior notching (published)": counts_for(sub, "exterior notching")[0],
               "paste: Bell share of rims (file only)": counts_for(sub, "paste: Bell share of rims in the two pastes")[0][:, ::-1]}
    rngd = np.random.default_rng(SEED + 77)
    res_d = {k: {"raw": [], "river": [], "line": []} for k in explain}
    dec_dist = []
    for _ in range(DRAWS):
        pd_ = np.array([rngd.dirichlet(r + 0.5) for r in dec])
        diss = 0.5 * np.abs(pd_[:, None, :] - pd_[None, :, :]).sum(-1)[iu19]
        dec_dist.append((spearman(diss, d_river), spearman(diss, d_line)))
        for k, m in explain.items():
            f = rngd.beta(m[:, 0] + 0.5, m[:, 1] + 0.5)
            gap = np.abs(f[:, None] - f[None])[iu19]
            res_d[k]["raw"].append(spearman(diss, gap))
            res_d[k]["river"].append(partial_spearman(diss, gap, d_river))
            res_d[k]["line"].append(partial_spearman(diss, gap, d_line))
    dd = np.array(dec_dist)
    f3 = lambda t: f"{t[0]:+.2f} ({t[1]:+.2f} to {t[2]:+.2f})"
    L += ["", "## D. Does any attribute go with the difference in decoration that distance leaves?", "",
          (f"For the {len(sub)} assemblages of part C ({len(d_river)} pairs). The difference in decoration between "
          "two assemblages is half the summed absolute difference of their decorated-class proportions. It is "
          f"correlated with river distance at {f3(q(dd[:, 0]))} and with straight-line distance at "
          f"{f3(q(dd[:, 1]))}. Each row asks whether pairs that differ more in an attribute also differ more in "
          "decoration, first as it stands and then with distance removed."), "",
          ("| attribute | rank correlation of decoration difference with attribute difference | with river distance "
          "removed | with straight-line distance removed |"), "|---|---|---|---|"]
    for k, r in res_d.items():
        L.append(f"| {k} | {f3(q(r['raw']))} | {f3(q(r['river']))} | {f3(q(r['line']))} |")
    L += ["", ("Interior bevel is read as increasing through time. A clear positive value in its row would say "
          "that assemblages of different age differ more in decoration than distance alone predicts, which is "
          "the reading that the unexplained difference is chronological. A value near zero says bevel, and so "
          "age as bevel measures it, does not go with that difference."), ""]

    L += ["", "## The attributes in Mainfort's file that the paper does not print", "",
          ("Paste is the paste of the rim sherds (author); lip form and bevel degree are as recorded. Counts are "
          "the file's percentages times N. Use with that in mind."), "", "Percent of rims, pooled by phase:", ""]
    L += phase_table(FILE_ONLY)
    L += [""] + head
    rows_file, res_file = line_rows([(a, counts_for(sub, a)[0]) for a in FILE_ONLY], 20)
    L += rows_file

    p_pub = {a: r["p_alt"] for a, r in res_pub.items() if not a.startswith("decorated")}
    L += ["", "## Reading", "",
          "**A.** " + "; ".join(f"{a}: difference grows with distance" for a in PUBLISHED
                                 if q([spearman(np.abs(p[:, None] - p[None])[iu], distR)
                                       for p in freq_draws[a][:300]])[1] > 0) + ".",
          "", f"**B.** {len(clear)} of the {len(conc)} attribute pairs have a rank correlation whose interval "
          "excludes zero" + (": " + "; ".join(f"{a} and {b} ({r[0]:+.2f})" for a, b, r in clear) if clear else "") + ".",
          "", f"**C.** Against random-center divisions, the phase lines separate the six published attributes "
          f"with probabilities from {min(p_pub.values()):.2f} to {max(p_pub.values()):.2f}, and the six together "
          f"with {tog['p_alt']:.2f}; the decorated types on the same assemblages give "
          f"{res_pub['decorated types, same assemblages']['p_alt']:.2f}. Of the unpublished attributes, "
          + "; ".join(f"{a} {r['p_alt']:.2f}" for a, r in res_file.items()) + ".", "",
          ("A restriction on learning at the phase lines predicts high probabilities in every attribute and "
          "attributes that vary together. Variation that arises locally predicts probabilities scattered "
          "around one half and attributes that vary independently. Cautions: part C rests on "
          f"{len(sub)} assemblages; rim form and paste may answer to vessel function, cost and time as well "
          "as to who learned from whom; and the table gives each attribute's total, not the attributes of "
          "each rim, so it cannot say which attributes occur on the same vessels."), ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
