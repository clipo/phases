"""Composition by collection (D36 closing step): the decomposition of the matrix into PFG parts and Lipo additions, and its summaries."""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT, ROOT / "analyses", ROOT / "analyses" / "probe_latent_groups", ROOT / "analyses" / "per_type", ROOT / "scripts", ROOT / "src"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
cbc = importlib.import_module("composition_by_collection")

# The author's ruling of 2026-09-21 (scripts/build_analysis_matrix.py docstring): Lipo's additions over the PFG counts.
RULED_ADDITIONS = {"Belle_Meade": 1455, "Beck": 1154, "Nickel": 1120, "Rose_Mound": 738, "Cramor_Place": 286, "Castile_Landing": 220}
# "A few sherds at three other sites" (same docstring), measured here.
MINOR_ADDITIONS = {"Barton_Ranch": 9, "Hollywood": 25, "Neeleys_Ferry": 4}


@pytest.fixture(scope="module")
def parts():
    return cbc.decompose()


def test_decomposition_reproduces_the_matrix_and_the_ruled_additions(parts):
    counts, pfg, add, cls, prov = parts
    assert len(counts) == 28 and (add >= 0).all().all()
    assert (pfg + add == counts).all().all()
    assert cls.value_counts().to_dict() == {"PFG": 18, "PFG+Lipo large": 6, "PFG+Lipo small": 3, "Lipo": 1}
    assert {a: int(add.loc[a].sum()) for a in cls.index[cls == "PFG+Lipo large"]} == RULED_ADDITIONS
    assert {a: int(add.loc[a].sum()) for a in cls.index[cls == "PFG+Lipo small"]} == MINOR_ADDITIONS
    assert len(prov) == 28 and prov.set_index("assemblage").loc["Kent_Place", "PFG rows"] == "13-N-4/B,C,D,E"
    assert prov.set_index("assemblage").loc["Nickel", "PFG rows"] == "13-N-15/B; 13-N-15/C,D,E"
    assert int(add.loc["Hollywood", "Painted"]) == 24                     # the one minor addition that could move a share
    assert (add.loc[cls.index[cls == "PFG"]].sum(1) == 0).all()
    assert (pfg.loc[[a for a in counts.index if a != "Holden_Lake"]].sum(1) >= cbc.MIN_SHERDS).all()   # all 27 PFG parts clear the cutoff
    assert cls.index[cls == "Lipo"].tolist() == ["Holden_Lake"] and int(pfg.loc["Holden_Lake"].sum()) == 0
    assert int(add.loc["Kent_Place"].sum()) == 0                       # collection C counted once (the ruled override)
    assert (pfg.loc[list(RULED_ADDITIONS)].sum(1) >= cbc.MIN_SHERDS).all()   # every PFG part stays above the cutoff


def test_spearman_by_class_direction_and_domain():
    # share of class 0 rises with total: rho positive; class 1 falls: negative; a class at one row only: NaN
    counts = pd.DataFrame({"a": [10, 30, 60, 100], "b": [40, 40, 40, 40], "c": [0, 0, 0, 5]})
    s = cbc.spearman_by_class(counts).set_index("class")
    assert s.loc["a", "rho"] > 0 and s.loc["b", "rho"] < 0 and np.isnan(s.loc["c", "rho"])   # a falling 'a' would mean shares and totals were paired wrongly
    two = cbc.spearman_by_class(counts.iloc[:2])
    assert two.rho.isna().all()


def test_log_odds_ratio_direction_and_zero_handling():
    rng = np.random.default_rng(0)
    r = cbc.log_odds_ratio(np.array([90, 10, 0]), np.array([50, 50, 0]), rng, draws=4000)
    assert r.loc[0, "median"] > 0 and r.loc[1, "median"] < 0                   # class 0 is a larger share of a: Lipo minus PFG positive
    assert r.loc[0, "q05"] > 0                                                 # 90 of 100 against 50 of 100 does not straddle zero
    assert np.isnan(r.loc[2, "median"])                                        # no sherds in either part
    same = cbc.log_odds_ratio(np.array([50, 50]), np.array([50, 50]), rng, draws=4000)
    assert same.loc[0, "q05"] < 0 < same.loc[0, "q95"]                        # probe point: equal shares straddle zero
    with pytest.raises(ValueError):
        cbc.log_odds_ratio(np.array([1, -1]), np.array([1, 1]), rng)


def test_decompose_writes_nothing(tmp_path):
    """CC-1: obtaining the decomposition must not write the legacy assemblage_sources.csv (or anything else under data/)."""
    import hashlib
    legacy = ROOT / "data" / "processed" / "assemblage_sources.csv"
    before = {p: (p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest()) for p in (ROOT / "data").rglob("*") if p.is_file()}
    cbc.decompose()
    after = {p: (p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest()) for p in (ROOT / "data").rglob("*") if p.is_file()}
    assert before == after and legacy in before


def test_matching_agrees_with_the_reconciliation_script(tmp_path):
    """The copied alias table and name matching give the same cpl-pfg id per assemblage as scripts/reconcile_ceramic_matrix.py,
    run in a sandbox so its own CSV write lands under tmp_path."""
    import contextlib
    import io
    import os
    import runpy
    (tmp_path / "data").mkdir(); (tmp_path / "data" / "processed").mkdir()
    os.symlink(ROOT / "data" / "raw", tmp_path / "data" / "raw")
    import shutil
    shutil.copy(ROOT / "data" / "processed" / "basin_members_curated.txt", tmp_path / "data" / "processed")   # an input the script reads
    here = os.getcwd(); os.chdir(tmp_path)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            g = runpy.run_path(str(ROOT / "scripts" / "reconcile_ceramic_matrix.py"))
    finally:
        os.chdir(here)
    assert (tmp_path / "data" / "processed" / "assemblage_sources.csv").exists()      # the write happened in the sandbox
    assert g["C_GRID"] == cbc.C_GRID
    tab, ids_all, g2n = cbc.load_sources()
    counts, *_ = cbc.decompose()
    for name in counts.index:
        k = g["norm"](name); gid = g["C_GRID"].get(k)
        theirs = [gid] if gid else [i for i in g["c"]["id"] if g["norm"](g["g2n"].get(i, i)) == k]
        assert cbc.tab_ids(name, ids_all, g2n) == theirs, name


def test_overlap_guard_rejects_shared_collection_letters():
    assert cbc.collection_letters("13-N-4/B,C,D,E") == frozenset("BCDE") and cbc.collection_letters("12-N-3/A&B") == frozenset("AB")
    assert cbc.collection_letters("13-O-7") == cbc.ALL_LETTERS                               # a whole-base row contains every component
    assert cbc.collection_letters("13-N-4/B,C,D,E") & cbc.collection_letters("13-N-4/C")       # the Kent overlap the override removes
    assert not (cbc.collection_letters("13-N-15/B") & cbc.collection_letters("13-N-15/C,D,E"))  # Nickel's rows are disjoint


def test_helpers_reject_invalid_counts_and_label_undefined_cases():
    rng = np.random.default_rng(0)
    for a, b in ((np.array([0, 0]), np.array([5, 5])), (np.array([1.5, 2]), np.array([5, 5])), (np.array([np.nan, 2]), np.array([5, 5])),
                 (np.array([3]), np.array([5])), (np.array([3, 2]), np.array([5, 5, 5]))):
        with pytest.raises(ValueError):
            cbc.log_odds_ratio(a, b, rng, draws=100)
    with pytest.raises(ValueError):
        cbc.log_odds_ratio(np.array([3, 2]), np.array([5, 5]), rng, draws=0)
    with pytest.raises(ValueError):
        cbc.log_odds_ratio(np.array([3, 2]), np.array([5, 5]), rng, draws=100, alpha=0.0)
    with pytest.raises(ValueError):
        cbc.spearman_by_class(pd.DataFrame({"a": [0, 1], "b": [0, 1]}))                      # an empty row
    with pytest.raises(ValueError):
        cbc.spearman_by_class(pd.DataFrame({"a": [1.5, 1, 2], "b": [1, 1, 1]}))              # fractional counts
    const = cbc.spearman_by_class(pd.DataFrame({"a": [5, 5, 5], "b": [5, 5, 5]})).set_index("class")
    assert const.loc["a", "undefined"] == "constant totals" and np.isnan(const.loc["a", "rho"])
    shares = cbc.spearman_by_class(pd.DataFrame({"a": [5, 10, 20], "b": [5, 10, 20]})).set_index("class")
    assert shares.loc["a", "undefined"] == "constant shares"


def test_log_odds_ratio_recovers_a_known_contrast():
    """Rule 20c: with counts drawn at known compositions, the posterior interval covers the true log odds ratio and the median is close."""
    rng = np.random.default_rng(7)
    pa = np.array([.6, .2, .1, .05, .05]); pb = np.array([.3, .4, .1, .1, .1])
    a = rng.multinomial(4000, pa); b = rng.multinomial(4000, pb)
    r = cbc.log_odds_ratio(a, b, np.random.default_rng(8), draws=20000)
    truth = np.log(pa / (1 - pa)) - np.log(pb / (1 - pb))
    assert np.all(r.q05.to_numpy() <= truth) and np.all(truth <= r.q95.to_numpy())
    assert np.all(np.abs(r["median"].to_numpy() - truth) < 0.12)
    # probe point: the prior alone (no data) would put every contrast near zero; a prior-dragged estimate at class 0 would miss +1.25
    assert r.loc[0, "q05"] > 0.9


def test_tab_id_must_occur_once_and_whole_base_rows_overlap_components():
    assert cbc.collection_letters("13-N-15") & cbc.collection_letters("13-N-15/B")       # a whole-base row contains every component
    import pandas as pd
    tab = pd.DataFrame({"id": ["13-N-4", "13-N-4", "11-N-1"]})
    with pytest.raises(ValueError):
        cbc.one_tab_id("Kent_Place", tab["id"].tolist(), {})                              # the alias maps to an id that occurs twice
    assert cbc.one_tab_id("Kent_Place", ["13-N-4", "11-N-1"], {}) == "13-N-4"


def test_log_odds_ratio_interval_coverage_at_small_ten_class_counts():
    """Rule 20c at the scale that matters: ten classes, 25 sherds against 200, the record's pooled shares as truth for both parts
    (true contrast 0 for every class). Over 300 replicates the 90 percent intervals of the common classes cover zero in a share
    consistent with 0.90 (binomial band 0.84 to 0.95); a prior that dragged the contrast would cover less."""
    rng = np.random.default_rng(21); draws = np.random.default_rng(22)
    p = np.array([7702, 4504, 1224, 150, 215, 138, 49, 67, 34, 18], float) / 14101
    covered = np.zeros(3); reps = 300
    for _ in range(reps):
        a = rng.multinomial(25, p); b = rng.multinomial(200, p)
        r = cbc.log_odds_ratio(a, b, draws, draws=2000)
        covered += ((r.q05.to_numpy()[:3] <= 0) & (0 <= r.q95.to_numpy()[:3])).astype(float)
    frac = covered / reps
    assert np.all(frac >= 0.84) and np.all(frac <= 0.97), frac
