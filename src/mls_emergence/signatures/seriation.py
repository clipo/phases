"""Deterministic frequency seriation (Signature 2).

Implements the deterministic core of iterative deterministic frequency
seriation (IDSS) following Lipo, Madsen & Dunnell (2015, PLOS ONE
10(4):e0124942). A set of assemblages (rows) is *co-seriable* when some
ordering of the rows renders every type column (relative frequencies)
unimodal: a single "battleship" peak (monotone non-decreasing then
non-increasing). IDSS recovers the set of MAXIMAL co-seriable orderings by
incremental construction (grow valid orderings by appending assemblages at
their ends) rather than by enumerating permutations.

The structurally important property reproduced here is that the maximal
co-seriable groups need not partition the assemblages: a single assemblage
may belong to two non-co-seriable groups at once (it "bridges" two otherwise
incompatible orderings). The count of maximal groups and the multi-membership
structure constitute Signature 2 (cf. the "Parkin in two lineages" hook in
Lipo, Madsen & Dunnell 2015).

Both deterministic criteria of the paper are enforced: joint unimodality and
continuity. Continuity here means adjacent assemblages must share at least one
co-present type (no turnover between neighbours with nothing in common) and no
per-type relative-frequency jump between neighbours may exceed ``cont`` (the
user-set continuity threshold of Lipo et al. 2015; default 1.0 forbids only a
complete turnover).

Scope note: this module implements only the DETERMINISTIC identity core. The
bootstrap confidence-interval significance layer of Lipo et al. (2015), which
distinguishes real coherence breakdown from sampling error, is OUT OF SCOPE
here. The ``tol`` parameter is a deterministic stand-in for that noise
tolerance: direction reversals smaller than ``tol`` (in relative frequency) are
treated as flat rather than as unimodality violations.
"""

from __future__ import annotations
import numpy as np


def _col_violations(freqs: np.ndarray) -> int:
    """Direction changes beyond a single peak in one frequency column.

    A perfectly unimodal column has at most one direction change (ascending
    then descending, or monotone in one direction).  Each change beyond the
    first counts as one violation.
    """
    diffs = np.sign(np.diff(freqs))
    diffs = diffs[diffs != 0]
    if len(diffs) < 2:
        return 0
    # Number of sign changes in the non-zero-difference sequence.
    n_changes = int((np.diff(diffs) != 0).sum())
    # A unimodal column has at most 1 change (rise then fall); violations are
    # any changes beyond that first one.
    return max(0, n_changes - 1)


def unimodality_violation(matrix: np.ndarray) -> int:
    """Total unimodality violations across type columns for the given row order."""
    freqs = matrix / matrix.sum(axis=1, keepdims=True)
    return int(sum(_col_violations(freqs[:, j]) for j in range(freqs.shape[1])))


def is_unimodal(freqs: np.ndarray, tol: float = 0.0) -> bool:
    """True if a 1-D frequency sequence has a single (battleship) peak.

    The sequence is unimodal when it rises (non-strictly) to a peak and then
    falls. A tolerance ``tol`` treats reversals smaller than ``tol`` as flat:
    a step counts as a rise only if it increases by more than ``tol``, and as a
    fall only if it decreases by more than ``tol``; intermediate steps are
    ignored. This absorbs sampling noise without admitting a genuine second
    peak.
    """
    f = np.asarray(freqs, dtype=float)
    if f.size < 3:
        return True
    # Reduce the sequence to its significant direction changes.
    seen_fall = False
    prev = f[0]
    for x in f[1:]:
        delta = x - prev
        if delta > tol:
            step = 1
        elif delta < -tol:
            step = -1
        else:
            step = 0
        if step == 1:
            if seen_fall:
                # Rising again after a fall -> second peak.
                return False
        elif step == -1:
            seen_fall = True
        prev = x
    return True


def _prep(matrix: np.ndarray):
    """Precompute per-row relative frequencies and presence mask (order-free)."""
    totals = matrix.sum(axis=1, keepdims=True)
    safe = np.where(totals == 0, 1.0, totals)
    freqs = matrix / safe
    present = matrix > 0
    zero_rows = (totals[:, 0] == 0)
    return freqs, present, zero_rows


def _order_is_valid(matrix: np.ndarray, order: list[int], tol: float = 0.0,
                    cont: float = 1.0, prep=None) -> bool:
    """True if a row ordering satisfies joint unimodality and continuity.

    Two conditions from Lipo, Madsen & Dunnell (2015) are enforced:

    1. Joint unimodality: every type column (relative frequencies) is unimodal
       (single peak), within ``tol``.
    2. Continuity: adjacent assemblages must share at least one co-present type
       (no abrupt turnover between neighbours that have nothing in common), and
       no per-type relative-frequency jump between neighbours may exceed
       ``cont``. The default ``cont=1.0`` forbids only a complete (100%)
       turnover; the adjacency-overlap requirement is always applied. The
       continuity threshold is the user-set parameter of the paper (it reports
       e.g. 0.30); smaller values admit fewer orders.

    ``prep`` is the optional output of :func:`_prep` (precomputed frequencies
    and presence mask) to avoid renormalizing on every call.
    """
    if len(order) <= 1:
        return True
    if prep is None:
        prep = _prep(matrix)
    freqs, present, zero_rows = prep
    if zero_rows[order].any():
        return False
    f = freqs[order]
    if not all(is_unimodal(f[:, j], tol) for j in range(f.shape[1])):
        return False
    p = present[order]
    co_present = p[:-1] & p[1:]
    if not co_present.any(axis=1).all():
        return False
    jumps = np.abs(np.diff(f, axis=0))
    if (jumps > cont).any():
        return False
    return True


def _canonical(order: list[int]) -> tuple[int, ...]:
    """Canonical key for an ordering and its reverse (a seriation = its reverse)."""
    rev = order[::-1]
    return tuple(order) if order <= rev else tuple(rev)


def _extend(matrix, order, remaining, tol, cont, found, visited, limits, prep):
    """Recursively grow ``order`` by appending remaining rows at its two ends.

    Records maximal valid orderings (nothing can be appended at either end) in
    ``found`` (a dict keyed by canonical form). This enumerates the
    deterministic solution set by incremental end-extension rather than by
    permuting all rows.

    ``visited`` memoizes the canonical form of every partial ordering already
    expanded, so a partial ordering reachable by several insertion sequences is
    expanded once. This turns the search over insertion *sequences* (factorial)
    into a search over distinct partial *orderings*, which keeps the
    construction tractable. ``limits`` is a dict with a ``max_visited`` budget
    and a ``truncated`` flag set when the budget is hit.

    Candidates are appended only at the two ENDS of the current order (not at
    interior positions). This is complete: any contiguous window of a valid
    order is itself a valid order (a contiguous slice of a unimodal sequence is
    unimodal, and adjacency overlap/continuity are local to neighbours), so
    every maximal order is reconstructed by growing outward from one of its
    adjacent pairs. Restricting growth to the ends bounds the partial orders to
    contiguous windows of maximal orders rather than the exponentially many
    interior-insertion paths, which is what the deterministic IDSS construction
    exploits.
    """
    if len(visited) >= limits["max_visited"]:
        limits["truncated"] = True
        return
    key = _canonical(order)
    if key in visited:
        return
    visited.add(key)
    extended = False
    for cand in remaining:
        for new_order in (order + [cand], [cand] + order):
            if _order_is_valid(matrix, new_order, tol, cont, prep):
                extended = True
                new_remaining = [r for r in remaining if r != cand]
                _extend(matrix, new_order, new_remaining, tol, cont, found,
                        visited, limits, prep)
                if limits["truncated"]:
                    return
    if not extended:
        # No assemblage can be appended at either end: this ordering is maximal
        # under end-extension. Record it; orderings that turn out to be subsets
        # of a larger group (reached from a different seed) are pruned once at
        # the end of ``seriation_solutions``.
        if key not in found:
            found[key] = list(order)


def seriation_solutions(matrix: np.ndarray, tol: float = 0.0,
                        max_solutions: int = 1000,
                        cont: float = 1.0,
                        max_visited: int = 200_000) -> list[list[int]]:
    """Maximal valid co-seriable orderings via deterministic iterative growth.

    Returns each maximal ordering as a list of row indices, reverse-deduplicated
    (an ordering and its reverse are the same seriation). Construction starts
    from every valid pair seed and grows by appending one assemblage at a time
    to either end while joint unimodality and continuity hold, until nothing can
    be appended.

    ``tol`` is the unimodality noise tolerance; ``cont`` is the continuity
    threshold (max adjacent per-type relative-frequency jump; adjacency overlap
    is always required). ``max_solutions`` caps the number of distinct maximal
    groups collected. ``max_visited`` caps the number of distinct partial
    orderings expanded (the search budget). Either cap being reached means the
    result is TRUNCATED and not guaranteed complete; a ``RuntimeWarning`` is
    emitted so the caller never silently trusts a partial enumeration.

    Scaling: the construction is exact and verified against brute-force
    enumeration on small cases. Cost is governed by the number of distinct valid
    partial orderings, which stays small when assemblages have narrow type
    overlap or the continuity threshold is tight (the regime of the worked LMV
    example in Lipo et al. 2015, where the largest solution held only four
    assemblages). The worst case (broad overlap and a permissive continuity
    threshold, so a single large group admits exponentially many sub-windows)
    grows quickly: clean, broadly overlapping synthetic battleships become slow
    around 16-20 assemblages under the default caps. This solver is therefore
    intended for small/medium sets (tens of assemblages) or for regional subsets
    of larger collections; it is not expected to enumerate ~233 assemblages in
    one pass. Tighten ``cont`` or raise the caps deliberately for larger runs.
    """
    matrix = np.asarray(matrix, dtype=float)
    n = matrix.shape[0]
    if n == 0:
        return []
    if n == 1:
        return [[0]]

    found: dict[tuple[int, ...], list[int]] = {}
    visited: set[tuple[int, ...]] = set()
    limits = {"max_visited": max_visited, "truncated": False}
    prep = _prep(matrix)

    # Seed from every valid unordered pair; the recursion appends at both ends,
    # so unordered pairs suffice as seeds.
    for i in range(n):
        for j in range(i + 1, n):
            seed = [i, j]
            if _order_is_valid(matrix, seed, tol, cont, prep):
                remaining = [r for r in range(n) if r not in (i, j)]
                _extend(matrix, seed, remaining, tol, cont, found,
                        visited, limits, prep)
                if limits["truncated"]:
                    break
        if limits["truncated"]:
            break

    # Handle isolates: rows that pair with nothing still form a singleton group.
    covered = set()
    for order in found.values():
        covered.update(order)
    for r in range(n):
        if r not in covered:
            key = _canonical([r])
            if key not in found:
                found[key] = [r]

    # Drop any maximal ordering that is a subset of a larger one (a strictly
    # larger group subsumes it; a smaller solution is not maximal as a group).
    solutions = list(found.values())
    sets = [frozenset(s) for s in solutions]
    keep = []
    for idx, s in enumerate(sets):
        if any(s < t for k, t in enumerate(sets) if k != idx):
            continue
        keep.append(solutions[idx])

    # Deduplicate distinct groups (by set) for the cap; distinct orderings of the
    # same set are one solution for capping purposes.
    distinct_groups = {frozenset(s) for s in keep}
    if limits["truncated"] or len(distinct_groups) > max_solutions:
        import warnings
        reason = ("search budget exhausted (max_visited reached)"
                  if limits["truncated"]
                  else f"more than max_solutions={max_solutions} groups found")
        warnings.warn(
            f"seriation_solutions may be incomplete: {reason}. "
            f"(max_solutions={max_solutions}, max_visited={max_visited}). "
            "Result is not guaranteed exhaustive; raise the caps or run on a "
            "subset.",
            RuntimeWarning,
            stacklevel=2,
        )
        if len(keep) > max_solutions:
            keep = keep[:max_solutions]
    return keep


def is_coseriable(matrix: np.ndarray, order: list[int] | None = None,
                  tol: float = 0.0, cont: float = 1.0) -> bool:
    """Test co-seriability.

    If ``order`` is given, test whether that specific ordering satisfies joint
    unimodality and continuity. If ``order`` is None, test whether ANY ordering
    of all rows is valid, using the iterative construction (a single maximal
    solution that contains every row witnesses full co-seriability).
    """
    matrix = np.asarray(matrix, dtype=float)
    if order is not None:
        return _order_is_valid(matrix, list(order), tol, cont)
    n = matrix.shape[0]
    if n <= 1:
        return True
    if n == 2:
        return _order_is_valid(matrix, [0, 1], tol, cont)
    for sol in seriation_solutions(matrix, tol=tol, cont=cont):
        if len(sol) == n:
            return True
    return False


def seriation_groups(matrix: np.ndarray, tol: float = 0.0,
                     cont: float = 1.0) -> dict:
    """Maximal co-seriable groups and their (possibly overlapping) membership.

    Returns a dict with:
      - ``n_groups``: number of maximal co-seriable groups (a group = the SET of
        assemblages in a maximal solution; distinct orderings over the same set
        count as one group).
      - ``membership``: mapping row index -> sorted list of group ids it belongs
        to. A row in more than one group is a bridge assemblage.
      - ``multi``: sorted list of row indices that belong to more than one group.

    Groups need not partition the rows; overlap is the Signature-2 structure.
    """
    matrix = np.asarray(matrix, dtype=float)
    n = matrix.shape[0]
    solutions = seriation_solutions(matrix, tol=tol, cont=cont)

    # Collapse orderings to their row sets; distinct orders of the same set are
    # one group. Then drop any set that is a subset of a larger group.
    seen: list[frozenset[int]] = []
    for sol in solutions:
        s = frozenset(sol)
        if s not in seen:
            seen.append(s)
    groups: list[frozenset[int]] = []
    for s in seen:
        if any(s < t for t in seen if t is not s):
            continue
        groups.append(s)

    membership: dict[int, list[int]] = {r: [] for r in range(n)}
    for gid, g in enumerate(groups):
        for r in sorted(g):
            membership[r].append(gid)
    multi = sorted(r for r, gids in membership.items() if len(gids) > 1)

    return {
        "n_groups": len(groups),
        "membership": membership,
        "multi": multi,
    }


def n_seriation_groups(matrix: np.ndarray, tol: float = 0.0,
                       cont: float = 1.0) -> int:
    """Number of maximal co-seriable groups (Signature 2 scalar).

    Re-implemented on the deterministic IDSS solver; returns
    ``seriation_groups(matrix, tol, cont)["n_groups"]``. Default ``tol=0.0``
    keeps the signature compatible with existing callers.
    """
    return seriation_groups(matrix, tol=tol, cont=cont)["n_groups"]
