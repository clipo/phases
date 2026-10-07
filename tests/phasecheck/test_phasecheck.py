"""phasecheck: the general tool must give the paper's numbers on the paper's record,
tell a real boundary from a gradient, and refuse input it cannot use."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import phasecheck as pc
from phasecheck import divisions as dv
from phasecheck import geo
from phasecheck import measures as ms
from phasecheck.cli import main as cli_main

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = ROOT / "examples" / "phasecheck" / "st_francis_basin.csv"
RIVER = ROOT / "examples" / "phasecheck" / "st_francis_river_km.csv"
ORDER = "seriation_order"


# ---------------------------------------------------------------- the paper's record
def test_reproduces_the_papers_phase_boundary_comparison():
    """Expected values are the paper's, from analyses/86 (output/findings/partition_posterior.md),
    an implementation this package does not call. Same seed, so agreement is to the printed digit."""
    data = pc.read_table(EXAMPLE, distance=RIVER, order_col=ORDER)
    r = pc.boundary_comparison(data, draws=2000, n_alt=300, seed=86000)
    assert [round(v, 4) for v in r["fst"]] == [0.0121, 0.0101, 0.0144]
    assert round(r["fst_around"][0], 4) == 0.0112 and round(r["fst_compact"][0], 4) == 0.0103
    assert (round(r["p_fst_around"], 2), round(r["p_fst_compact"], 2)) == (0.57, 0.72)
    assert [round(v, 1) for v in r["excess"]] == [7.3, 4.6, 9.8]          # positive: more alike within a phase
    assert (round(r["p_excess_around"], 2), round(r["p_excess_compact"], 2)) == (0.58, 0.65)
    assert round(r["fst_direct"], 4) == 0.0123 and round(r["excess_direct"], 1) == 7.7


def test_reproduces_the_papers_recovery_fit_and_profiles():
    """Expected values from output/findings/phase_recovery.csv (all 28 scored), assemblage_phase_fit.md
    and phase_profiles.md."""
    data = pc.read_table(EXAMPLE, order_col=ORDER)
    rec = {r["weight"]: r for r in pc.recovery_by_location(data, kmeans_seeds=3)["rows"]}
    assert round(rec[0.0]["ward"], 4) == 0.7556 and round(rec[0.0]["kmeans"], 4) == 0.5087
    assert round(rec[0.25]["kmeans"], 4) == 0.7556
    assert round(rec[np.inf]["ward"], 4) == 0.3260
    assert rec[0.0]["ward_misplaced"] == ["Big_Eddy", "Castile_Landing"]
    fit = pc.own_phase_fit(data, draws=4000, seed=96)
    assert sum(r["p_own"] > 0.75 for r in fit) == 22
    assert sorted(r["name"] for r in fit if r["p_own"] < 0.25) == ["Belle_Meade", "Big_Eddy", "Holden_Lake", "Turnbow"]
    prof = pc.phase_profiles(data)
    assert round(prof["mean"], 1) == 19.4
    assert {k: round(v, 1) for k, v in prof["pairs"].items()} == {
        ("Kent", "Parkin"): 12.5, ("Kent", "Walls"): 24.2, ("Parkin", "Walls"): 21.6}


# ---------------------------------------------------------------- constructed records
def _constructed(kind: str, seed: int = 5):
    """36 assemblages scattered on a 60 km square, cut into three groups of 12 by arbitrary lines.

    kind="boundary": each group has its own class mix, the same everywhere inside it.
    kind="gradient": the class mix changes smoothly from west to east and ignores the lines.
    Same map, same lines, same sherd counts: only what the pottery follows differs.
    """
    rng = np.random.default_rng(seed)
    n, k = 36, 6
    east, north = rng.uniform(0, 60, n), rng.uniform(0, 60, n)
    lat, lon = 35.0 + north / 111.32, -90.0 + east / (111.32 * np.cos(np.deg2rad(35.27)))
    pts = geo.planar_km(lat, lon)
    idx = dv.assign_exact(pts, pts[rng.choice(n, 3, replace=False)], np.repeat(np.arange(3), 12))
    base = np.array([[.40, .25, .15, .10, .06, .04], [.10, .15, .40, .25, .04, .06], [.06, .04, .10, .15, .25, .40]])
    counts = np.empty((n, k), int)
    for i in range(n):
        if kind == "boundary":
            p = base[idx[i]]
        else:
            t = east[i] / 60.0
            p = (1 - t) * base[0] + t * base[2]
        counts[i] = rng.multinomial(300, p)
    df = pd.DataFrame(counts, columns=[f"class_{j + 1}" for j in range(k)])
    df.insert(0, "phase", [f"P{g + 1}" for g in idx]); df.insert(0, "longitude", lon)
    df.insert(0, "latitude", lat); df.insert(0, "name", [f"site_{i + 1:02d}" for i in range(n)])
    return df


def _uneven(df):
    """The same table with sherd totals that differ between assemblages, as real collections do."""
    out = df.copy()
    cols = [c for c in out.columns if c.startswith("class_")]
    out[cols] = out[cols].mul(1 + np.arange(len(out)) % 4, axis=0)
    return out


def _load(df, tmp_path, name="t.csv", **kw):
    p = tmp_path / name
    df.to_csv(p, index=False)
    return pc.read_table(p, **kw)


def test_a_real_boundary_is_found_and_a_gradient_is_not(tmp_path):
    """The discriminating pair, over six maps each. Lines that the pottery follows stand above nearly
    every alternative on every map. Arbitrary lines through a gradient stand anywhere, by how they
    happen to run, so across maps their standing averages near one half. A tool that called every
    spatial pattern a boundary would pass the first half and fail the second; one that never found a
    boundary would do the reverse."""
    seeds = range(5, 11)
    b = [pc.boundary_comparison(_load(_constructed("boundary", s), tmp_path), draws=150, n_alt=60, seed=1)
         for s in seeds]
    g = [pc.boundary_comparison(_load(_constructed("gradient", s), tmp_path), draws=150, n_alt=60, seed=1)
         for s in seeds]
    assert min(r["p_fst_around"] for r in b) > 0.9 and min(r["p_excess_around"] for r in b) > 0.9
    assert min(r["excess"][0] for r in b) > 50          # a large step at the lines
    assert 0.2 < np.mean([r["p_fst_around"] for r in g]) < 0.8
    assert max(abs(r["excess"][0]) for r in g) < 30     # no step of the boundary's size
    fit = pc.own_phase_fit(_load(_constructed("boundary"), tmp_path), draws=200, seed=2)
    assert all(r["p_own"] > 0.9 for r in fit)


def test_boundary_excess_direction():
    """Within minus between. Positive when same-group pairs are more alike; the opposite sign would
    mean assemblages are MORE alike across a line than on one side of it."""
    pts = np.column_stack([np.repeat(np.arange(6.0), 2), np.tile([0.0, 1.0], 6)])
    dist = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
    a, b = [90, 10], [10, 90]
    counts = np.array([a, b] * 6)                      # the class mix alternates north and south
    by_mix = np.tile([0, 1], 6)                        # groups that follow the mix
    across_mix = np.repeat([0, 1], 6)                  # groups that cut across it
    assert ms.boundary_excess(counts, dist, by_mix, min_pairs=1) > 100
    assert ms.boundary_excess(counts, dist, across_mix, min_pairs=1) < 0
    assert ms.profile_difference([100, 0], [0, 100]) == 100.0 and ms.profile_difference([60, 40], [60, 40]) == 0.0


def test_fst_and_agreement_on_known_cases():
    assert ms.cultural_fst(np.array([[50, 50], [50, 50]])) == pytest.approx(0.0)
    assert ms.cultural_fst(np.array([[100, 0], [0, 100]])) == pytest.approx(1.0)
    assert ms.adjusted_rand([0, 0, 1, 1], [1, 1, 0, 0]) == pytest.approx(1.0)     # labels do not matter
    assert ms.adjusted_rand([0, 0, 1, 1], [0, 1, 0, 1]) < 0.01
    d = geo.great_circle_km(np.array([35.0, 36.0]), np.array([-90.0, -90.0]))
    assert d[0, 1] == pytest.approx(111.2, abs=0.2) and d[0, 0] == 0.0


def test_alternative_divisions_keep_the_group_sizes():
    rng = np.random.default_rng(0)
    pts = rng.uniform(0, 50, (20, 2))
    idx = np.repeat([0, 1, 2], [9, 7, 4])
    around, compact = dv.alternatives(pts, idx, 15, rng)
    for a, c in zip(around, compact):
        assert sorted(np.bincount(a)) == [4, 7, 9] and sorted(np.bincount(c)) == [4, 7, 9]
        assert dv.spread(pts, c) <= dv.spread(pts, a) + 1e-9      # compact is never looser


# ---------------------------------------------------------------- input that must be refused
@pytest.mark.parametrize("change, message", [
    (lambda d: d.assign(class_1=-1), "negative values"),
    (lambda d: d.assign(class_1=d.class_1 / d.class_1.sum()), "whole numbers"),
    (lambda d: d.assign(latitude=d.longitude, longitude=d.latitude), "swapped"),
    (lambda d: d.assign(phase="one"), "single phase"),
    (lambda d: d.assign(name="same"), "unique"),
    (lambda d: d.drop(columns="phase"), "missing column"),
    (lambda d: d.assign(phase=["solo"] + list(d.phase[1:])), "at least two assemblages"),
    (lambda d: d.assign(latitude=np.linspace(20, 50, len(d))), "span about"),
    (lambda d: d.assign(class_2=[np.nan] + list(d.class_2[1:])), "blank or non-numeric"),
    # columns that are not sherd counts: each would give a plausible, wrong report if read as a class
    (lambda d: d.assign(site_number=np.arange(1, len(d) + 1)), "looks like an identifier"),
    (lambda d: d.assign(surveyed=np.linspace(1940, 1990, len(d)).astype(int)), "nearly constant"),
    (lambda d: d.assign(year=np.arange(len(d)) * 7 % 50 + 1940), "named like something other"),
    (lambda d: d.assign(mound=np.arange(len(d)) % 2), "presence flag"),
    (lambda d: _uneven(d).assign(all_sherds=_uneven(d).filter(like="class_").sum(axis=1)), "looks like a total"),
    (lambda d: d.assign(remarks="surface collection"), "holds text"),
    (lambda d: d.assign(**{c: (100 * d[c] / 300).round().astype(int) for c in d.filter(like="class_")}),
     "looks like percentages"),
])
def test_bad_input_is_refused_with_a_reason(tmp_path, change, message):
    with pytest.raises(ValueError, match=message):
        _load(change(_constructed("gradient")), tmp_path)


def test_extra_columns_can_be_left_out_or_classes_named(tmp_path):
    df = _constructed("gradient").assign(site_number=np.arange(1, 37), remarks="surface collection")
    a = _load(df, tmp_path, ignore_cols=["site_number", "remarks"])
    b = _load(df, tmp_path, class_cols=[f"class_{j}" for j in range(1, 7)])
    assert a.classes == b.classes == [f"class_{j}" for j in range(1, 7)]
    assert np.array_equal(a.counts, b.counts)
    with pytest.raises(ValueError, match="more than once"):
        _load(df, tmp_path, class_cols=["class_1", "class_1", "class_2"])


def test_rows_of_about_100_sherds_can_be_declared_counts(tmp_path):
    df = _constructed("gradient")
    cols = [c for c in df.columns if c.startswith("class_")]
    df[cols] = (100 * df[cols] / 300).round().astype(int)
    with pytest.raises(ValueError, match="looks like percentages"):
        _load(df, tmp_path)
    assert len(_load(df, tmp_path, totals_of_100_are_counts=True).names) == 36


def test_a_region_across_the_date_line_is_not_read_as_the_globe():
    lat = np.array([-17.0, -17.2, -16.8, -17.1])
    lon = np.array([179.8, -179.9, 179.6, -179.7])
    p = geo.planar_km(lat, lon)
    assert np.ptp(p[:, 0]) < 80 and np.ptp(p[:, 1]) < 60
    assert p[1, 0] > p[0, 0]            # -179.9 lies east of 179.8


def test_phase_named_na_and_windows_text_are_read(tmp_path):
    df = _constructed("gradient")
    df["phase"] = df["phase"].replace({"P1": "NA", "P2": "Peña"})
    p = tmp_path / "w.csv"
    df.to_csv(p, index=False, encoding="cp1252")
    data = pc.read_table(p)
    assert sorted(data.phase_names) == ["NA", "P3", "Peña"]


def test_boundary_excess_says_when_it_has_no_distance_control():
    """Three groups of two give no distance bin five pairs of each kind. The value is then the raw gap,
    and bins_used must say 0 so the report can say so."""
    rng = np.random.default_rng(3)
    pts = rng.uniform(0, 30, (6, 2))
    dist = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
    counts = rng.multinomial(200, [.5, .3, .2], size=6)
    value, used = ms.boundary_excess_detail(counts, dist, np.repeat([0, 1, 2], 2))
    assert used == 0 and value == ms.boundary_excess(counts, dist, np.repeat([0, 1, 2], 2))
    big = rng.multinomial(200, [.5, .3, .2], size=30)
    pts = rng.uniform(0, 30, (30, 2))
    dist = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
    assert ms.boundary_excess_detail(big, dist, np.repeat([0, 1, 2], 10))[1] >= 2


def test_fit_names_the_likeliest_phase_not_only_the_own_probability(tmp_path):
    """With several phases, a low probability for the own phase does not make another phase the best fit.
    The row must carry which phase wins most draws, and that can be the assemblage's own."""
    fit = pc.own_phase_fit(_load(_constructed("gradient"), tmp_path), draws=300, seed=4)
    for r in fit:
        assert 0 <= r["p_own"] <= r["p_likeliest"] <= 1
        assert (r["likeliest"] == r["phase"]) == (r["p_own"] == r["p_likeliest"])


def test_distance_matrix_must_match_the_assemblages(tmp_path):
    df = _constructed("gradient")
    names = list(df.name)
    good = pd.DataFrame(np.abs(np.subtract.outer(np.arange(36.0), np.arange(36.0))), index=names, columns=names)
    good.to_csv(tmp_path / "d.csv", index_label="name")
    assert _load(df, tmp_path, distance=tmp_path / "d.csv").dist_kind.startswith("supplied")
    good.iloc[:-1, :-1].to_csv(tmp_path / "short.csv", index_label="name")
    with pytest.raises(ValueError, match="no row or column"):
        _load(df, tmp_path, distance=tmp_path / "short.csv")
    lop = good.copy(); lop.iloc[0, 1] = 99.0
    lop.to_csv(tmp_path / "lop.csv", index_label="name")
    with pytest.raises(ValueError, match="not symmetric"):
        _load(df, tmp_path, distance=tmp_path / "lop.csv")


def test_minimum_count_drops_and_lists_small_assemblages(tmp_path):
    df = _constructed("gradient")
    df.loc[0, [c for c in df.columns if c.startswith("class_")]] = [3, 2, 1, 0, 0, 0]
    data = _load(df, tmp_path, min_count=50)
    assert data.dropped == [("site_01", 6)] and len(data.names) == 35


def test_command_line_writes_a_report(tmp_path):
    p = tmp_path / "in.xlsx"
    _constructed("boundary").to_excel(p, index=False)
    out = tmp_path / "rep"
    assert cli_main([str(p), "--out", str(out), "--draws", "120", "--alternatives", "30"]) == 0
    text = (out / "report.md").read_text()
    assert "## 2. Do the phase boundaries separate the assemblages better than other lines?" in text
    assert (out / "assemblage_fit.csv").exists() and (out / "results.json").exists()


def test_command_line_reports_bad_input_without_a_traceback(tmp_path, capsys):
    p = tmp_path / "bad.csv"
    _constructed("gradient").assign(phase="one").to_csv(p, index=False)
    assert cli_main([str(p), "--out", str(tmp_path / "rep")]) == 2
    assert "single phase" in capsys.readouterr().err
    assert cli_main([str(tmp_path / "no_such_file.csv")]) == 2
    assert "no such file" in capsys.readouterr().err
