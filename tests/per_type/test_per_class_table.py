"""The per-class table (spec 6.5): a class absent everywhere is reported absent, never nan, never dropped."""
import importlib
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "analyses" / "per_type", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
pct = importlib.import_module("per_class_table")


def test_absent_class_is_reported_absent_not_nan():
    y = np.array([[5, 0, 0], [3, 0, 2]])
    row = pct.per_class_rows(y, labels=np.array([0, 1]))
    assert row["occ_1"] == 0 and not np.isfinite(row["F_1"])
    assert pct.fmt_f(row["F_1"]) == "absent"


from mls_emergence.probe import rivals

POOLED = np.array([7702, 4504, 1224, 150, 215, 138, 49, 67, 34, 18], float) / 14101
LABELS_11 = np.array([1] * 11 + [0] * 17)


def test_planted_f_rises_on_a_homogeneous_profile():
    """Direction test (rule 21a): on a profile identical at every site, planting class 1 at m = 3 in the 11-site cluster
    raises its moment F from about zero to the value the planted share gives by hand. On a heterogeneous profile F can
    FALL when the marked cluster starts with less of the class (Codex RD-4), so a fall on generator records is not a
    software defect; this test uses the homogeneous case where the direction is not in doubt."""
    prof = np.tile(POOLED, (28, 1)); c, m = 1, 3.0
    planted = rivals.plant_marker(prof, LABELS_11, 1, c, m)
    y0 = np.rint(prof * 1e6).astype(int); y1 = np.rint(planted * 1e6).astype(int)
    f0 = pct.moment_f(y0)[c]; f1 = pct.moment_f(y1)[c]
    p = POOLED[c]; share = m * p / (1 + (m - 1) * p)                      # by hand, not through plant_marker
    shares = np.array([share] * 11 + [p] * 17); mean = shares.mean()
    expected = shares.var(ddof=1) / (mean * (1 - mean))
    assert abs(f0) < 1e-3
    assert f1 > 0.05 and abs(f1 - expected) < 1e-3, (f1, expected)


def test_planted_occupancy_expectation_is_nondecreasing():
    """Expected occupancy sum_i 1 - (1 - p_ic)^n_i cannot fall under m > 1; strictly rises when a marked site with a
    positive total has 0 < p < 1; equal when every marked site's share is 0 or 1."""
    def expected_occ(prof, totals, c):
        return float((1 - (1 - prof[:, c]) ** totals).sum())
    c = 2; totals = np.array([50, 50, 50, 50])
    prof = np.array([[.5, .3, .2], [.1, .6, .3], [.25, .25, .5], [.7, .2, .1]]); labels = np.array([0, 1, 1, 0])
    e1 = expected_occ(rivals.plant_marker(prof, labels, 1, c, 1.0), totals, c)
    e3 = expected_occ(rivals.plant_marker(prof, labels, 1, c, 3.0), totals, c)
    assert e3 > e1
    boundary = np.array([[.5, .3, .2], [.4, .6, .0], [.0, .0, 1.0], [.7, .2, .1]])   # marked rows have p_c = 0 and 1
    b1 = expected_occ(rivals.plant_marker(boundary, labels, 1, c, 1.0), totals, c)
    b3 = expected_occ(rivals.plant_marker(boundary, labels, 1, c, 3.0), totals, c)
    assert b3 == b1


# ---- regional-drift spec tests 7.9, 7.10, 7.14 (C3 and M in the table) ----

def _tiny_config():
    """A four-site configuration so the C3 and M paths run in a test without the record."""
    rng = np.random.default_rng(3)
    pts = rng.uniform(0, 30, (4, 2)); D = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
    pct.G.update(N=4, K=10, D=D, pooled=POOLED, totals=np.array([200, 300, 150, 400]), labels=np.array([1, 0, 0, 1]), records=2)


def test_m1_equals_the_neutral_record_exactly_and_pairs_across_markers():
    _tiny_config()
    op = {"kernel": "gaussian", "length": 6.0, "n_ind": 50, "innovation": 0.05, "mixing": 0.3}
    sc = {"L_d": 5, "L_R": 8, "W": 64, "S": 70, "burnin": 30}
    kids = pct.record_streams("primary", 3, 0)
    g3 = pct.generator_spec("C3", "0.1", sc["L_d"], sc["L_R"])
    y3, ex3 = pct.one_record_c3m(op, g3, sc, kids, marker=None)
    gm = pct.generator_spec("M", "0.1/1", sc["L_d"], sc["L_R"])
    y1, _ex1 = pct.one_record_c3m(op, gm, sc, kids, marker={"c": 1, "cluster": 1})
    assert np.array_equal(y3, y1)                                              # m = 1 on the shared stream is the neutral record
    gm3 = pct.generator_spec("M", "0.1/3", sc["L_d"], sc["L_R"])
    ya, exa = pct.one_record_c3m(op, gm3, sc, kids, marker={"c": 1, "cluster": 1})
    _yb, exb = pct.one_record_c3m(op, gm3, sc, kids, marker={"c": 2, "cluster": 0})
    for key in ("ends", "lengths", "targets_at_ends"):
        assert np.array_equal(exa[key], exb[key]) and np.array_equal(exa[key], ex3[key])   # history, spans and targets shared
    assert not np.array_equal(ya, y3)                                          # the marker alone changed the counts
    assert np.isfinite(ex3["profile_f_before"]).any() and "lost" in ex3 and ex3["lost"].shape == (10,)


def test_c3_and_m_settings_parse():
    from mls_emergence.transmission.spatial import kappa_for_window_f
    g = pct.generator_spec("C3", "0.1", 182, 292)
    assert g["W"] == 2336 and g["F_W"] == 0.1 and np.isclose(g["kappa"], kappa_for_window_f(0.1, 2336)) and g["m"] is None
    gm = pct.generator_spec("M", "0.1/2", 182, 292)
    assert gm["kappa"] == g["kappa"] and gm["m"] == 2.0 and gm["F_W"] == 0.1
    g0 = pct.generator_spec("C3", "0", 182, 292)
    assert g0["kappa"] is None and g0["F_W"] == 0.0
    import pytest
    for bad in ("0.1/2/3", "abc", "0.1/", "/2", "0.1/0", "1.5"):
        with pytest.raises(ValueError):
            pct.generator_spec("M", bad, 182, 292)
    with pytest.raises(ValueError):
        pct.generator_spec("C3", "0.2", 182, 292)        # not on the registered amplitude grid


def test_frozen_scales_and_lengths():
    p = pct.frozen_scales("primary"); s = pct.frozen_scales("secondary")
    assert (p["L_d"], p["L_R"], p["W"], p["S"], p["burnin"]) == (182, 292, 2336, 2600, 5000)
    assert (s["L_d"], s["L_R"], s["W"], s["S"], s["burnin"]) == (12, 40, 320, 340, 5000)
    assert p["S"] - p["W"] == 264 and p["S"] - p["W"] - p["L_d"] + 1 == 83
    assert s["S"] - s["W"] == 20 and s["S"] - s["W"] - s["L_d"] + 1 == 9


def test_record_streams_are_keyed_by_point_amplitude_and_record():
    a = pct.record_streams("primary", 3, 7); b = pct.record_streams("primary", 3, 7); c = pct.record_streams("secondary", 3, 7)
    assert a["root"] == b["root"] == 2026101950 + 3000 + 7 and c["root"] == 2026101950 + 100000 + 3000 + 7
    assert [k.spawn_key for k in a["children"]] == [(0,), (1,), (2,), (3,)]
    assert a["sim_seed"] == b["sim_seed"] != c["sim_seed"]


def test_moment_f_removes_binomial_sampling_variance():
    """Identical site shares at small totals: the corrected F is near 0; without the correction it is near 1/n;
    a correction divided by n instead of n - 1 leaves a bias of about 1/n^2 (0.0625 at n = 4)."""
    rng = np.random.default_rng(5); n = 4; sites = 4000
    y = rng.multinomial(n, [0.5, 0.5], size=sites)
    f = pct.moment_f(y)[0]
    P = y / n; uncorrected = P[:, 0].var(ddof=1) / (P[:, 0].mean() * (1 - P[:, 0].mean()))
    assert abs(f) < 0.02 and uncorrected > 0.2                   # the /n mutant gives about +0.06 here; no correction about +0.25


def test_generator_spec_rejects_a_setting_for_the_neutral_generator():
    import pytest
    with pytest.raises(ValueError):
        pct.generator_spec("neutral", "2", 182, 292)
    assert pct.generator_spec("neutral", "neutral", 182, 292)["kind"] == "neutral"


def test_paired_marker_record_leaves_unmarked_sites_bitwise_identical():
    """Per-site sherd streams (closing step, 2026-10-19): a planted site must not redraw the sites after it, so every unmarked
    site's counts in an M record equal the paired C3 record's exactly, and only marked sites differ."""
    rng = np.random.default_rng(4); n_sites = 16
    pts = rng.uniform(0, 40, (n_sites, 2)); D = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
    labels = np.array([1] * 8 + [0] * 8)                       # the marked sites come first in the draw order
    pct.G.update(N=n_sites, K=10, D=D, pooled=POOLED, totals=rng.integers(200, 1500, size=n_sites), labels=labels, records=2)
    op = {"kernel": "gaussian", "length": 6.0, "n_ind": 50, "innovation": 0.05, "mixing": 0.3}
    sc = {"L_d": 5, "L_R": 8, "W": 64, "S": 70, "burnin": 30}
    diffs = 0
    for r in range(6):
        kids = pct.record_streams("primary", 3, r)
        y3, _ = pct.one_record_c3m(op, pct.generator_spec("C3", "0.1", sc["L_d"], sc["L_R"]), sc, kids, marker=None)
        ym, _ = pct.one_record_c3m(op, pct.generator_spec("M", "0.1/3", sc["L_d"], sc["L_R"]), sc, kids, marker={"c": 0, "cluster": 1})
        assert np.array_equal(y3[labels == 0], ym[labels == 0])   # with one shared stream, most of these rows differ
        diffs += int(not np.array_equal(y3[labels == 1], ym[labels == 1]))
    assert diffs == 6                                           # the marked sites do change


def test_per_site_streams_are_the_sherds_childs_spawns():
    kids = pct.record_streams("primary", 0, 0)
    direct = np.random.default_rng(np.random.SeedSequence(kids["root"], spawn_key=(3, 5))).integers(0, 10**9, 4)
    spawned = np.random.default_rng(kids["children"][3].spawn(6)[5]).integers(0, 10**9, 4)
    assert np.array_equal(direct, spawned)


def test_moment_f_has_one_definition_in_the_signatures_package():
    from mls_emergence.signatures import variance as V
    assert pct.moment_f is V.moment_f
    import pytest
    with pytest.raises(ValueError):
        V.moment_f(np.array([[1, 0], [3, 2]]))                                   # a site with one observation
