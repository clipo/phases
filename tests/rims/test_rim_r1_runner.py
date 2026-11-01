"""Rim R1 runner (spec section 3): record keys, per-site streams, 28 populations transmitting with 19 observed, frozen scales,
binary and reference cells, output naming. Tests written before the runner."""
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "analyses" / "rims", ROOT / "analyses" / "per_type", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses", ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
r1 = importlib.import_module("rim_r1")


def _config(seed=0, n=12, observed=(0, 2, 3, 5, 7, 8, 10, 11)):
    rng = np.random.default_rng(seed)
    pts = rng.uniform(0, 40, (n, 2)); D = np.sqrt(((pts[:, None] - pts[None]) ** 2).sum(-1))
    pooled = np.array([7702, 4504, 1224, 150, 215, 138, 49, 67, 34, 18], float) / 14101
    r1.G.update(N=n, K=10, D=D, pooled=pooled, observed=np.array(observed), rim_totals=rng.integers(27, 278, len(observed)),
                sherd_totals=rng.integers(79, 1710, len(observed)), pfg_totals=rng.integers(79, 700, len(observed)),
                pooled_pfg=pooled, records=2)
    return {"L_d": 5, "L_R": 8, "W": 64, "S": 70, "burnin": 30}


def test_record_key_is_a_tuple_entropy_distinct_from_the_integer_roots():
    k = r1.record_streams("primary", "rim", 0, 0, 0)
    assert k["entropy"] == (2026101970, 0, 0, 0, 0, 0)
    assert [c.spawn_key for c in k["children"]] == [(0,), (1,), (2,), (3,)]
    assert k["sim_seed"] != np.random.default_rng(np.random.SeedSequence(2026101950).spawn(4)[0]).integers(0, 2**31 - 1)
    k2 = r1.record_streams("primary", "rim", 0, 1, 0)
    assert k2["sim_seed"] != k["sim_seed"]                                      # the two windows share nothing
    with pytest.raises(ValueError):
        r1.record_streams("tertiary", "rim", 0, 0, 0)


def test_binary_cell_observes_only_the_rim_sites_while_all_populations_transmit():
    sc = _config(); op = {"kernel": "gaussian", "length": 6.0, "n_ind": 50, "innovation": 0.05, "mixing": 0.3}
    kids = r1.record_streams("primary", "rim", 1, 0, 0)
    y, ex = r1.one_binary_record(op, 0.12, sc, 64, kids)
    assert y.shape == (len(r1.G["observed"]), 2) and np.array_equal(y.sum(1), r1.G["rim_totals"])
    assert ex["weights_shape"] == (12, 12) and ex["profile_shape"] == (12, 2)   # 12 populations transmit; 8 are observed
    f = r1.record_f_binary(y)
    assert np.isfinite(f)


def test_per_site_streams_in_the_binary_cell():
    sc = _config(); op = {"kernel": "gaussian", "length": 6.0, "n_ind": 50, "innovation": 0.05, "mixing": 0.3}
    kids = r1.record_streams("primary", "rim", 1, 0, 3)
    y1, _ = r1.one_binary_record(op, 0.12, sc, 64, kids)
    saved = r1.G["rim_totals"].copy(); r1.G["rim_totals"] = saved.copy(); r1.G["rim_totals"][0] = saved[0] + 50   # change one site's total
    y2, _ = r1.one_binary_record(op, 0.12, sc, 64, kids)
    r1.G["rim_totals"] = saved
    assert np.array_equal(y1[1:], y2[1:]) and not np.array_equal(y1[0], y2[0])   # other sites' draws are untouched


def test_reference_cell_gives_ten_class_counts_at_sherd_totals_and_common_class_f():
    sc = _config(); op = {"kernel": "gaussian", "length": 6.0, "n_ind": 50, "innovation": 0.05, "mixing": 0.3}
    kids = r1.record_streams("primary", "decoration_current", 0, 1, 0)
    y, _ex = r1.one_reference_record(op, r1.G["pooled"], r1.G["sherd_totals"], sc, 5, kids)
    assert y.shape == (8, 10) and np.array_equal(y.sum(1), r1.G["sherd_totals"])
    f = r1.record_f_common(y)
    assert f.shape == (3,)


def test_output_names_carry_kind_and_point_and_scales_come_frozen():
    assert r1.output_stem("rim", "primary") == "rim_r1_rim_primary" and r1.output_stem("decoration_pfg", "secondary") == "rim_r1_decoration_pfg_secondary"
    sc = r1.scales_for("primary")
    assert (sc["L_d"], sc["L_R"], sc["S"], sc["burnin"]) == (182, 292, 2600, 5000) and r1.window_steps("T1", sc) == 182 and r1.window_steps("8LR", sc) == 2336
    sc2 = r1.scales_for("secondary"); assert (sc2["S"], r1.window_steps("8LR", sc2)) == (340, 320)
    with pytest.raises(ValueError):
        r1.window_steps("4LR", sc)
