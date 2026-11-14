"""The drift-extension run script (spec 2026-11-14, draft 3, sections 3, 5, 7 and 8).

Expected values are literals worked by hand, taken from the committed findings
files, or taken from the auditor's independent recomputation (Codex spec audit
2026-11-14, report line 186, and repair check, report line 95). None is
computed through the function under test. Rule 3: each test names the rival
it rejects.
"""
import importlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "analyses" / "per_type", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
de = importlib.import_module("drift_extension")

from mls_emergence.probe import spans


@pytest.fixture(scope="module")
def ctx():
    return de.context()


# ---- the grid, the scales and the keys ---------------------------------------------------------

def test_grid_sizes():
    assert de.AMPLITUDES == (0.0, 0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9)
    assert len(de.DESIGNS) == 63 and len(set(de.DESIGNS)) == 63
    assert len(de.DESIGNS) * len(de.AMPLITUDES) == 504
    assert de.D0 == ("Ld", "8LR", "c0") and de.D0 in de.DESIGNS


@pytest.mark.parametrize("point, window, W, S", [
    # W = multiple of L_R (292 primary, 40 secondary); S = W + L_d (182, 12) rounded up to 100, or to 10 below 1,000.
    ("primary", "8LR", 2336, 2600), ("primary", "2LR", 584, 770), ("primary", "32LR", 9344, 9600),
    ("secondary", "8LR", 320, 340), ("secondary", "2LR", 80, 100), ("secondary", "32LR", 1280, 1300)])
def test_scales_by_hand(point, window, W, S):
    sc = de.scales(point, window)
    assert (sc["W"], sc["S"], sc["burnin"]) == (W, S, 5000)
    # Rival rejected: a window fixed at 8 L_R whatever the window name.


def test_scales_raise_when_a_time_scale_file_disagrees(tmp_path, monkeypatch):
    f = tmp_path / "time_scales_primary.json"
    f.write_text(json.dumps({"L_d": 180, "L_R": 292, "burnin_used": 5000}))
    monkeypatch.setitem(de.TS_FILES, "primary", f)
    with pytest.raises(ValueError):
        de.scales("primary", "8LR")
    f.write_text(json.dumps({"L_d": 182, "L_R": 292, "burnin_used": 5000}))
    assert de.scales("primary", "8LR")["L_d"] == 182


def test_span_limits_by_hand():
    # (L_d, ceil(L_d / 4), 1): ceil(182 / 4) = ceil(45.5) = 46; ceil(12 / 4) = 3.
    assert [de.span_limit(s, 182) for s in de.SPANS] == [182, 46, 1]
    assert [de.span_limit(s, 12) for s in de.SPANS] == [12, 3, 1]


def test_kappa_is_recomputed_for_the_window():
    # F_W = 1 - (kappa / (kappa + 1))^W. At W = 584 and F_W = 0.3, kappa = 1 / expm1(-log(0.7) / 584).
    k2 = de.kappa("primary", "2LR", 4); k8 = de.kappa("primary", "8LR", 4)
    assert k2 == pytest.approx(1 / math.expm1(-math.log(0.7) / 584), rel=1e-12)
    assert 1 - (k2 / (k2 + 1)) ** 584 == pytest.approx(0.3, abs=1e-9)
    assert 1 - (k8 / (k8 + 1)) ** 2336 == pytest.approx(0.3, abs=1e-9)
    assert k8 > 3.9 * k2          # rival rejected: one kappa reused across windows (ratio 1)
    assert de.kappa("primary", "8LR", 0) is None


def test_seed_key_is_built_from_the_registered_tuple():
    k = de.seed_key("secondary", "32LR", 7, 49, "coupling")
    assert k.entropy == 2026111400 and tuple(k.spawn_key) == (1, 2, 7, 49, 4)
    assert tuple(de.seed_key("primary", "2LR", 3, 5, "sherds", site=27).spawn_key) == (0, 1, 3, 5, 3, 27)
    for bad in ({"a": 8}, {"r": 50}, {"r": -1}, {"child": "nope"}, {"point": "tertiary"}, {"window": "4LR"}):
        args = {"point": "primary", "window": "8LR", "a": 0, "r": 0, "child": "history"}; args.update(bad)
        with pytest.raises(ValueError):
            de.seed_key(**args)
    with pytest.raises(ValueError):
        de.seed_key("primary", "8LR", 0, 0, "sherds")            # sherds need a site
    with pytest.raises(ValueError):
        de.seed_key("primary", "8LR", 0, 0, "lengths", site=3)   # a site index on a record-level stream


def test_base_uniforms_are_pinned_literals():
    # Computed once under NumPy 2.5.4 with Generator.random on a fresh generator of the key, and recorded here.
    b = de.base_draws("primary", "8LR", 0, 0, 28)
    np.testing.assert_allclose(b["v"][:3], [0.09562106, 0.97733475, 0.93100093], atol=5e-9)
    b = de.base_draws("secondary", "32LR", 7, 49, 28)
    np.testing.assert_allclose(b["u"][:3], [0.83291756, 0.76402573, 0.55400354], atol=5e-9)
    assert len(b["v"]) == len(b["u"]) == len(b["w"]) == len(b["ends0"]) == 28


# ---- observation designs -----------------------------------------------------------------------

def test_c0_returns_the_window_draw_unchanged_and_needs_no_coupling_uniforms(ctx):
    sc = de.scales("primary", "8LR")
    expected = spans.ends_window(np.random.default_rng(np.random.SeedSequence(2026111400, spawn_key=(0, 0, 2, 7, 2))), 28, sc["S"], sc["W"])
    b = de.base_draws("primary", "8LR", 2, 7, 28)
    reduced = {"v": b["v"], "ends0": b["ends0"]}              # no u, no w: c0 must not read them
    ends, lengths = de.design_spans(ctx, "primary", "8LR", "Ld", "c0", reduced)
    assert ends.tolist() == expected.tolist()
    assert lengths.tolist() == (1 + np.floor(b["v"] * 182)).astype(int).tolist()


def test_latitude_designs_reassign_the_same_ends(ctx):
    b = de.base_draws("primary", "8LR", 2, 7, 28); lat = ctx["lat"]; order = np.argsort(lat)
    got = {c: de.design_spans(ctx, "primary", "8LR", "Ld", c, b)[0] for c in ("c0", "lat50N", "lat100N", "lat50S", "lat100S")}
    for c, e in got.items():
        assert sorted(e.tolist()) == sorted(got["c0"].tolist()), c
    assert (np.diff(got["lat100N"][order]) >= 0).all()          # later to the north
    assert (np.diff(got["lat100S"][order]) <= 0).all()          # later to the south
    assert got["lat50N"].tolist() != got["lat50S"].tolist() and got["lat50N"].tolist() != got["lat100N"].tolist()


def test_step_designs_put_the_named_cluster_late_and_share_the_uniforms(ctx):
    b = de.base_draws("primary", "8LR", 2, 7, 28); labels = ctx["labels"]
    large = labels == ctx["large_label"]
    assert int(large.sum()) == 17 and int((~large).sum()) == 11
    eL = de.design_spans(ctx, "primary", "8LR", "Ld", "stepL", b)[0]
    eS = de.design_spans(ctx, "primary", "8LR", "Ld", "stepS", b)[0]
    assert eL[large].min() > eL[~large].max()                    # stepL: the larger cluster is late
    assert eS[~large].min() > eS[large].max()                    # stepS: the larger cluster is early
    # By hand: S = 2600, W = 2336, window 264..2599, early half 264..1431 (1168 steps), late half 1432..2599 (1168).
    np.testing.assert_array_equal(eL[large], 1432 + np.floor(b["w"][large] * 1168).astype(int))
    np.testing.assert_array_equal(eS[large], 264 + np.floor(b["w"][large] * 1168).astype(int))
    # The same w places a site at the same offset in whichever half its cluster holds.
    np.testing.assert_array_equal(eL - eS, np.where(large, 1168, -1168))


def test_span_lengths_share_v_across_limits(ctx):
    b = de.base_draws("primary", "8LR", 2, 7, 28)
    L182, L46, L1 = (de.design_spans(ctx, "primary", "8LR", s, "c0", b)[1] for s in de.SPANS)
    np.testing.assert_array_equal(L46, 1 + np.floor(b["v"] * 46).astype(int))
    assert L1.tolist() == [1] * 28 and L182.max() > 46
    assert (np.diff(L182[np.argsort(b["v"])]) >= 0).all()         # the same v orders the sites under every limit


def test_unknown_names_and_misfits_raise(ctx):
    b = de.base_draws("primary", "8LR", 2, 7, 28)
    with pytest.raises(ValueError):
        de.design_spans(ctx, "primary", "8LR", "Ld", "lat75N", b)
    with pytest.raises(ValueError):
        de.design_spans(ctx, "primary", "8LR", "Ld/2", "c0", b)
    with pytest.raises(ValueError):
        de.check_fit(2336, 300, 2600)          # window plus span longer than the stored record
    de.check_fit(2336, 182, 2600)
    for bad in (-0.1, 1.0, 0.2, float("nan")):  # outside [0, 1), or not on the registered grid
        with pytest.raises(ValueError):
            de.amplitude_index(bad)
    assert de.amplitude_index(0.7) == 6
    three = dict(ctx); three["labels"] = np.arange(28) % 3
    with pytest.raises(ValueError):
        de.design_spans(three, "primary", "8LR", "Ld", "stepL", b)


def test_numpy_version_guard():
    with pytest.raises(RuntimeError):
        de.require_numpy("2.4.6")
    de.require_numpy("2.5.4")


# ---- records: pairing, order invariance, per-site sherd streams ----------------------------------

def test_no_drift_keeps_the_target_at_the_pooled_shares(ctx):
    rec, targets = de.simulate(ctx, "secondary", "2LR", 0, 0)
    assert rec.shape == (100, 28, 10) and targets.shape == (100, 10)
    np.testing.assert_allclose(targets, np.tile(ctx["pooled"], (100, 1)))
    _rec, drifted = de.simulate(ctx, "secondary", "2LR", 7, 0)
    assert np.abs(drifted[-1] - ctx["pooled"]).max() > 0.01      # F_W = 0.9 moves it; a path that ignored kappa would not


def test_lat100N_and_c0_share_the_history_and_differ_in_site_ends(ctx):
    rec_a, tg_a = de.simulate(ctx, "secondary", "8LR", 5, 3)
    rec_b, tg_b = de.simulate(ctx, "secondary", "8LR", 5, 3)
    assert np.array_equal(rec_a, rec_b) and np.array_equal(tg_a, tg_b)
    b = de.base_draws("secondary", "8LR", 5, 3, 28)
    e0 = de.design_spans(ctx, "secondary", "8LR", "Ld", "c0", b)[0]
    e1 = de.design_spans(ctx, "secondary", "8LR", "Ld", "lat100N", b)[0]
    assert e0.tolist() != e1.tolist()
    other = de.simulate(ctx, "secondary", "8LR", 5, 4)[0]
    assert not np.array_equal(rec_a, other)                       # a different record is a different history


@pytest.mark.parametrize("window", ["8LR", "2LR", "32LR"])
def test_observations_do_not_depend_on_the_order_of_designs(ctx, window):
    designs = [d for d in de.DESIGNS if d[1] == window]
    assert len(designs) == 21
    fwd = de.record_counts(ctx, "secondary", window, 5, 0, designs)
    rev = de.record_counts(ctx, "secondary", window, 5, 0, designs[::-1])
    for d in designs:
        alone = de.record_counts(ctx, "secondary", window, 5, 0, [d])[d]
        assert np.array_equal(fwd[d], alone) and np.array_equal(rev[d], alone), d
    # Rival rejected: a generator object carried from one design to the next, which makes a design's
    # counts depend on which designs came before it.
    assert not np.array_equal(fwd[designs[0]], fwd[designs[-1]])


def test_each_site_draws_its_sherds_from_its_own_stream(ctx):
    rec, targets = de.simulate(ctx, "secondary", "8LR", 5, 0)
    b = de.base_draws("secondary", "8LR", 5, 0, 28)
    ends, lengths = de.design_spans(ctx, "secondary", "8LR", "Ld", "c0", b)
    y0, _ = de.observe(ctx, rec, targets, ends, lengths, "secondary", "8LR", 5, 0)
    moved = ends.copy(); moved[0] = ends[0] - 7 if ends[0] - 7 - lengths[0] + 1 >= 0 else ends[0] + 7
    y1, _ = de.observe(ctx, rec, targets, moved, lengths, "secondary", "8LR", 5, 0)
    assert np.array_equal(y0[1:], y1[1:])                         # a change at site 0 changes no other site
    assert (y0.sum(1) == ctx["totals"]).all()
    # Each site's counts are the multinomial draw of a fresh generator of the key (point 1, window 0, amplitude 5,
    # record 0, child 3, site i), built here from the literal key and not through seed_key.
    prof = spans.span_profile(rec, ends, lengths)
    for i in (0, 1, 13, 27):
        g = np.random.default_rng(np.random.SeedSequence(2026111400, spawn_key=(1, 0, 5, 0, 3, i)))
        assert y0[i].tolist() == g.multinomial(int(ctx["totals"][i]), prof[i]).tolist(), i
    # Rival rejected: one stream consumed across the sites in order, under which only site 0 would match.


# ---- readouts and cells --------------------------------------------------------------------------

def test_record_readouts_match_the_auditors_independent_values(ctx):
    r = de.readouts(ctx, ctx["y"])
    # Codex repair check 2026-11-14, report line 95, computed without this code.
    assert r["F_0"] == pytest.approx(0.1646173519, abs=1e-9)
    assert r["F_1"] == pytest.approx(0.1172237789, abs=1e-9)
    assert r["F_2"] == pytest.approx(0.1986379181, abs=1e-9)
    assert r["lat_0"] == pytest.approx(+0.5402298851, abs=1e-9)   # Parkin Punctated: higher to the north
    assert r["lat_1"] == pytest.approx(+0.0574712644, abs=1e-9)
    assert r["lat_2"] == pytest.approx(-0.7750410509, abs=1e-9)   # painted types: higher to the south
    assert r["D_dist"] == pytest.approx(+0.3599509600, abs=1e-9)  # farther pairs differ more
    assert r["lat_reason_0"] == "" and r["D_reason"] == ""


def _cell(values, col="F_0"):
    rows = [{"point": "primary", "span": "Ld", "window": "8LR", "coupling": "c0", "F_W": 0.3, "record": i, col: v, "lost_0": i < 4}
            for i, v in enumerate(values)]
    return pd.DataFrame(rows)


def test_band_statuses_by_hand():
    # Values 0, 1, ..., 50 (51 values, so that both quantile positions are whole numbers and the endpoints are exact):
    # linear quantile at 0.1 is position 5 -> 5; at 0.9 is position 45 -> 45. Endpoints are inclusive.
    v = np.arange(51.0)
    for x, want in ((45.5, "above"), (20.0, "inside"), (4.5, "below"), (5.0, "inside"), (45.0, "inside"), (float("nan"), "unavailable")):
        row = de.stat_row(_cell(v), "F_0", x)
        assert row["status"] == want, (x, row)
    assert (row["lo"], row["median"], row["hi"], row["n_finite"], row["n_lost"]) == (5.0, 25.0, 45.0, 51, 4)
    # Rivals rejected: exclusive endpoints (5.0 would be below); a band of the minimum and maximum (4.5 would be inside).


def test_eligibility_boundary_is_45_finite_values():
    v = np.arange(50) / 100.0
    v44 = v.copy(); v44[:6] = np.nan
    v45 = v.copy(); v45[:5] = np.nan
    r44 = de.stat_row(_cell(v44), "F_0", 0.9); r45 = de.stat_row(_cell(v45), "F_0", 0.9)
    assert (r44["n_finite"], r44["status"]) == (44, "unavailable")
    assert (r45["n_finite"], r45["status"]) == (45, "above")


def test_a_lat_is_banded_on_absolute_values():
    # Signed values symmetric about zero: the band of the absolute value is not the band of the signed value.
    v = np.linspace(-0.49, 0.49, 50)
    signed = de.stat_row(_cell(v, "lat_0"), "lat_rank_0", 0.45)
    unsigned = de.stat_row(_cell(v, "lat_0"), "A_lat_0", 0.45)
    assert signed["lo"] < 0 < signed["hi"] and unsigned["lo"] > 0
    assert unsigned["median"] == pytest.approx(np.median(np.abs(v)))
    assert de.stat_row(_cell(v, "lat_0"), "A_lat_0", 0.0)["status"] == "below"   # a record with no trend sits below


def test_committed_table_parser_reads_the_committed_primary_cell():
    # output/findings/per_type_generators_C3_primary.md, row "C3 0.3": the committed cells of the three common classes.
    c = de.committed_cells("primary")
    assert c[0.3]["F_0"] == {"median": 0.035, "lo": 0.011, "hi": 0.077, "n_finite": 50, "status": "outside"}
    assert c[0.3]["F_1"]["hi"] == 0.111 and c[0.3]["F_2"]["median"] == 0.029
    assert c[0.3]["lost"][:3] == [2, 6, 32] and c[0.0]["lost"] == [0] * 10
    assert sorted(c) == [0.0, 0.01, 0.03, 0.1, 0.3]


def _synthetic_records(real):
    """A full point of per-record rows with planted values: at coupling c0 and F_W = 0.9 every compared statistic is
    spread around the record's own value (inside); everywhere else the simulated values are near zero (record above)."""
    import itertools
    near0 = np.linspace(0.001, 0.002, 50); spread = np.linspace(-0.05, 0.05, 50); rows = []
    for (span, window, coupling), (a, fw) in itertools.product(de.DESIGNS, enumerate(de.AMPLITUDES)):
        hit = coupling == "c0" and fw == 0.9
        for r in range(50):
            row = {"point": "primary", "span": span, "L": 1, "window": window, "W": 100, "S": 200, "coupling": coupling, "F_W": fw,
                   "kappa": 1.0, "record": r, "D_dist": real["D_dist"] + spread[r] if hit else near0[r]}
            for k in range(10):
                row[f"F_{k}"] = real[f"F_{k}"] + spread[r] if hit and k < 3 else near0[r]
                row[f"lat_{k}"] = real[f"lat_{k}"] + spread[r] if hit else (near0[r] if r % 2 else -near0[r])
                row[f"occ_{k}"] = 20; row[f"lost_{k}"] = r < 7; row[f"pfT_{k}"] = 0.1; row[f"pfB_{k}"] = 0.1
            rows.append(row)
    return pd.DataFrame(rows)


def test_cell_table_and_findings_section_on_planted_records(ctx):
    real = de.readouts(ctx, ctx["y"])
    t = de.cell_table(_synthetic_records(real), real)
    assert len(t) == 504 * len(de.stat_names()) and set(t["point"]) == {"primary"}
    hit = t[(t.coupling == "c0") & (t.F_W == 0.9) & t.stat.isin(["F_0", "F_1", "F_2", "D_dist", "A_lat_0", "A_lat_1", "A_lat_2"])]
    assert len(hit) == 9 * 7 and (hit.status == "inside").all()
    miss = t[(t.coupling != "c0") & t.stat.isin(["F_0", "F_1", "F_2"])]
    assert (miss.status == "above").all() and (miss.n_lost == 7).all()
    assert (t[t.stat.str.startswith("pf")].status == "not compared").all()
    repro = {"F_W": 0.3, "cells": {f"F_{k}": {"median": 0.0, "lo": 0.0, "hi": 0.0} for k in range(3)}}
    text = "\n".join(de.point_section("primary", t, ctx, real, repro))
    # Planted: L = S = the 9 cells at c0 and F_W 0.9 (3 span limits x 3 windows), resting on 3 x 50 histories.
    assert "**Statement (form c).**" in text and "at 9 cells" in text
    assert "(L): 9, resting on 150 distinct simulated histories" in text
    assert "| c0 | 72 | 72 | 9 | 9 | 9 |" in text and "| stepS | 72 | 72 | 0 | 0 | 0 |" in text
    assert "inside at F_W 0.9; above at 0, 0.01, 0.03, 0.1, 0.3, 0.5, 0.7; below at none; unavailable at none." in text


@pytest.mark.parametrize("point", ["secondary", "primary"])
def test_the_committed_c3_cell_is_reproduced_in_this_environment(ctx, point):
    """Under NumPy 2.5.4 the committed C3 cell at F_W 0.3 regenerates to the committed three decimals; under another
    version the multinomial draws differ (measured 2026-11-14 under 2.4.6 at the primary point), so the test is
    skipped there and the run itself refuses to start."""
    if np.__version__ != de.NUMPY_REQUIRED:
        pytest.skip(f"the committed tables were produced under NumPy {de.NUMPY_REQUIRED}")
    got = de.reproduce_committed(ctx, point)
    assert got["cells"]["F_0"]["n_finite"] == 50
    if point == "primary":      # literals from output/findings/per_type_generators_C3_primary.md, row C3 0.3
        assert (got["cells"]["F_0"]["median"], got["cells"]["F_0"]["lo"], got["cells"]["F_0"]["hi"]) == (0.035, 0.011, 0.077)
        assert got["lost"][:3] == [2, 6, 32]


def test_reproduction_check_raises_when_the_frozen_scales_are_wrong(ctx, monkeypatch):
    # With L_d moved from 12 to 14 at the secondary point the regenerated cell must not pass as the committed one.
    if np.__version__ != de.NUMPY_REQUIRED:
        pytest.skip(f"the committed tables were produced under NumPy {de.NUMPY_REQUIRED}")
    monkeypatch.setitem(de.FROZEN, "secondary", {"L_d": 14, "L_R": 40, "burnin": 5000})
    with pytest.raises(RuntimeError):
        de.reproduce_committed(ctx, "secondary")


# ---- a stochastic probe of one registered setting (spec section 8) -------------------------------

def test_probe_latitude_ordered_timing_raises_distance_decay_at_one_setting(ctx):
    """Primary point, F_W = 0.5, W = 8 L_R, L = 1, records 0 to 19: D_dist under lat100N minus D_dist under c0.

    One-sided. The other side would mean that, at this setting, ordering
    occupation dates by latitude did not raise D_dist relative to
    location-independent timing. It would say nothing about whether distance
    decay is present under either design, since both values can be positive.
    This checks the code path at one setting, not the design's premise.
    """
    a = de.amplitude_index(0.5); diffs = []
    for r in range(20):
        ys = de.record_counts(ctx, "primary", "8LR", a, r, [("1", "8LR", "c0"), ("1", "8LR", "lat100N")])
        d0 = de.readouts(ctx, ys[("1", "8LR", "c0")])["D_dist"]; d1 = de.readouts(ctx, ys[("1", "8LR", "lat100N")])["D_dist"]
        if np.isfinite(d0) and np.isfinite(d1):
            diffs.append(d1 - d0)
    assert len(diffs) >= 15, f"unavailable: only {len(diffs)} of 20 pairs have both values finite"
    assert np.median(diffs) > 0, f"median paired difference {np.median(diffs):+.3f} over {len(diffs)} pairs"
