"""Ranked readouts (type-origination spec 2026-11-14, draft 1, section 5)."""
import numpy as np
import pytest

from mls_emergence.dataio.basin import geo_km, load_curated
from mls_emergence.signatures import ranked_readouts as rr


def test_rank_order_by_hand_with_a_tie_and_an_absent_class():
    # Pooled counts by column: 5, 9, 0, 9, 1 -> columns 1 and 3 tie at 9 (column order), then 0, then 4; column 2 is dropped.
    y = np.array([[2, 4, 0, 9, 1], [3, 5, 0, 0, 0]])
    assert rr.rank_order(y).tolist() == [1, 3, 0, 4]


def test_fewer_than_ten_classes_gives_zero_occupancy_at_missing_ranks():
    y = np.array([[50, 30, 15, 5, 0], [40, 40, 20, 0, 0], [60, 20, 10, 9, 1]])
    lat = np.array([1.0, 2.0, 3.0]); D = np.abs(lat[:, None] - lat[None, :])
    out = rr.ranked_readouts(y, lat, D)
    # Ranks by pooled count: columns 0 (150), 1 (90), 2 (45), 3 (14), 4 (1). Rank 4 is present at 2 sites, rank 5 at 1.
    assert out["K_obs"] == 5 and out["occ_r4"] == 2 and out["occ_r5"] == 1
    assert [out[f"occ_r{r}"] for r in range(6, 11)] == [0] * 5 and out["beyond10"] == 0.0
    two = rr.ranked_readouts(y[:, :2], lat, D)
    assert np.isnan(two["F_r3"]) and np.isnan(two["lat_r3"]) and two["K_obs"] == 2


def test_share_beyond_rank_ten():
    y = np.ones((2, 12), int); y[:, 0] = 89                       # 100 sherds per site, classes 11 and 12 hold 1 each
    assert rr.ranked_readouts(y, np.array([1.0, 2.0]), np.array([[0, 1], [1, 0.0]]))["beyond10"] == pytest.approx(0.02)


def test_the_record_by_rank():
    counts, coords = load_curated(); xy = coords.to_numpy(float)
    out = rr.ranked_readouts(counts.to_numpy(int), xy[:, 0], geo_km(xy))
    # Ranks 1 to 3 are Parkin Punctated, Barton/Kent/MPI and the painted types: the values of the Codex repair check
    # 2026-11-14 (report line 95), computed without this code.
    assert out["F_r1"] == pytest.approx(0.1646173519, abs=1e-9) and out["F_r3"] == pytest.approx(0.1986379181, abs=1e-9)
    assert out["lat_r1"] == pytest.approx(+0.5402298851, abs=1e-9) and out["lat_r3"] == pytest.approx(-0.7750410509, abs=1e-9)
    assert out["D_dist"] == pytest.approx(+0.3599509600, abs=1e-9)
    # Ranks 4 to 10 by pooled count (215, 150, 138, 67, 49, 34, 18: Ranch, Fortune, Walls, Rhodes, Wallace, Vernon Paul,
    # Hull), with the site counts of output/findings/drift_extension.md: 22, 15, 14, 15, 7, 10, 8.
    assert [out[f"occ_r{r}"] for r in range(4, 11)] == [22, 15, 14, 15, 7, 10, 8]
    assert out["K_obs"] == 10 and out["beyond10"] == 0.0
