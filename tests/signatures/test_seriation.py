import numpy as np
from mls_emergence.signatures import seriation as S


def test_unimodality_zero_for_battleship_order():
    M = np.array([[1,0,0],[3,1,0],[1,3,1],[0,1,3],[0,0,1]], float)
    assert S.unimodality_violation(M) == 0

def test_unimodality_positive_for_scrambled_order():
    M = np.array([[1,0,0],[3,1,0],[1,3,1],[0,1,3],[0,0,1]], float)
    scrambled = M[[0,3,1,4,2]]
    assert S.unimodality_violation(scrambled) > 0

def test_one_group_when_coseriable():
    coseriable = np.array([[3,0],[2,1],[1,2],[0,3]], float)
    assert S.n_seriation_groups(coseriable) == 1


# --- is_unimodal ---------------------------------------------------------

def test_is_unimodal_single_peak():
    assert S.is_unimodal(np.array([0.0, 0.3, 0.6, 0.3, 0.0]))
    assert S.is_unimodal(np.array([0.1, 0.2, 0.7]))          # monotone rise
    assert S.is_unimodal(np.array([0.7, 0.2, 0.1]))          # monotone fall

def test_is_unimodal_rejects_two_peaks():
    assert not S.is_unimodal(np.array([0.0, 0.5, 0.1, 0.5, 0.0]))

def test_is_unimodal_rejects_valley():
    # increasing toward both ends is forbidden (Lipo et al. 2015)
    assert not S.is_unimodal(np.array([0.6, 0.1, 0.6]))


def test_violation_score_rejects_valleys_and_matches_predicate():
    assert S.unimodality_violation(np.array([[8, 2], [2, 8], [8, 2]])) == 1
    rng = np.random.default_rng(4)
    for _ in range(200):
        f = rng.integers(0, 5, size=8)
        assert (S._col_violations(f) == 0) == S.is_unimodal(f)
        assert S._col_violations(f) == S._col_violations(f[::-1])

def test_is_unimodal_tolerance_absorbs_small_dip():
    seq = np.array([0.0, 0.4, 0.38, 0.6, 0.0])  # tiny dip at index 2
    assert not S.is_unimodal(seq, tol=0.0)
    assert S.is_unimodal(seq, tol=0.05)


# --- is_coseriable -------------------------------------------------------

def test_is_coseriable_given_order():
    M = np.array([[1,0,0],[3,1,0],[1,3,1],[0,1,3],[0,0,1]], float)
    assert S.is_coseriable(M, order=[0,1,2,3,4])
    assert not S.is_coseriable(M, order=[0,3,1,4,2])

def test_is_coseriable_any_order():
    # rows given scrambled; some ordering is valid
    M = np.array([[1,0,0],[3,1,0],[1,3,1],[0,1,3],[0,0,1]], float)
    assert S.is_coseriable(M[[2,0,4,1,3]])


# --- seriation_groups: the decisive structural cases ---------------------

def test_single_battleship_is_one_group():
    M = np.array([[1,0,0],[3,1,0],[1,3,1],[0,1,3],[0,0,1]], float)
    g = S.seriation_groups(M)
    assert g["n_groups"] == 1
    assert all(len(v) == 1 for v in g["membership"].values())

def test_two_disjoint_sequences_are_two_groups():
    A = np.array([[4,0,0,0],[2,2,0,0],[0,4,0,0]], float)      # types 0,1
    B = np.array([[0,0,4,0],[0,0,2,2],[0,0,0,4]], float)      # types 2,3
    M = np.vstack([A, B])
    g = S.seriation_groups(M)
    assert g["n_groups"] == 2
    assert g["multi"] == []   # clean partition, no bridge

def test_bridge_assemblage_belongs_to_two_groups():
    # Branching-lineage construction (cf. Lipo, Madsen & Dunnell 2015, Figs 7-8;
    # the "Parkin in two lineages" hook). A stem assemblage (row 0) leads to a
    # branch-point assemblage (row 1) that carries incipient amounts of BOTH
    # divergent decorated types. From there two lineages diverge: tip 2 (type t1)
    # and tip 3 (type t2). The tips share no type, so they can never be adjacent
    # in a seriation; they can only connect through the branch point.
    #
    #         t0 t1 t2
    #  row0  [ 4, 0, 0]   stem base   (overlaps only the branch point via t0)
    #  row1  [ 2, 1, 1]   branch point (t0 + both incipient tip types)  <- BRIDGE
    #  row2  [ 0, 4, 0]   tip of lineage A
    #  row3  [ 0, 0, 4]   tip of lineage B
    #
    # Each lineage is individually co-seriable, but the union {0,1,2,3} has NO
    # single valid order: the branch point would have to be adjacent to all three
    # other assemblages at once (the only assemblage each of them overlaps), which
    # a linear order cannot provide. The branch point therefore belongs to more
    # than one maximal group -- multi-membership, not a partition.
    M = np.array([[4, 0, 0],
                  [2, 1, 1],
                  [0, 4, 0],
                  [0, 0, 4]], float)

    # Each lineage seriates on its own; the full set does not.
    assert S.is_coseriable(M[[0, 1, 2]])
    assert S.is_coseriable(M[[0, 1, 3]])
    assert not S.is_coseriable(M)

    g = S.seriation_groups(M)
    # More than one maximal group, and at least one assemblage bridges them.
    assert g["n_groups"] > 1
    assert len(g["multi"]) >= 1
    # The branch point (row 1) is a bridge: it appears in more than one group.
    assert len(g["membership"][1]) >= 2
    assert 1 in g["multi"]
