"""Record keys, named random streams, and paired accumulation spans (spec 3.3, 3.4; R2-4).

One base key per scenario and record, with no span-setting term. Named child
streams are spawned from it by a fixed name-to-integer table, never by the
order in which code happens to consume one stream, so adding a stream later
cannot shift the others. The history, partition, span lengths and S6 offsets
come from setting-independent streams and are shared by T0 and T1 of a
record; the T1 end jitter, subsampling, folds, Stan seeds and conditional
draws come from setting-dependent streams.
"""
from __future__ import annotations

import numpy as np

from mls_emergence.transmission.spatial import sample_record

BASE = 2026101600
SCENARIOS = ("S0", "S1", "S2", "S3", "S4", "S5", "S6", "S0-stress")

# Setting-independent streams (shared by T0 and T1).
STREAMS = {"history": 1, "partition": 2, "lengths": 3, "offsets": 4}
# Setting-dependent streams: the integer is offset by 100 * t (t = 0 for T0, 1 for T1).
SETTING_STREAMS = {"ends": 11, "subsample": 12, "folds": 13, "stan": 14, "conditional_a": 15, "conditional_b": 16}


def base_key(scenario: int | str, record: int) -> int:
    """2026101600 + 10000 s + r, with s the scenario index (0..7) and r the record index."""
    s = SCENARIOS.index(scenario) if isinstance(scenario, str) else int(scenario)
    if not 0 <= s < len(SCENARIOS) or record < 0 or record >= 10000:
        raise ValueError(f"scenario index {s} or record {record} out of range")
    return BASE + 10000 * s + int(record)


def stream(key: int, name: str, setting: int | None = None) -> np.random.Generator:
    """A named child generator of the base key. Setting-dependent names need `setting` (0 or 1)."""
    if name in STREAMS:
        if setting is not None:
            raise ValueError(f"stream {name!r} is shared across settings; do not pass a setting")
        sid = STREAMS[name]
    elif name in SETTING_STREAMS:
        if setting not in (0, 1):
            raise ValueError(f"stream {name!r} depends on the span setting; pass setting 0 (T0) or 1 (T1)")
        sid = SETTING_STREAMS[name] + 100 * setting
    else:
        raise ValueError(f"unknown stream name {name!r}")
    return np.random.default_rng(np.random.SeedSequence(int(key), spawn_key=(sid,)))


def draw_lengths(rng: np.random.Generator, n_sites: int, L_d: int) -> np.ndarray:
    """Per-site span lengths, uniform on the integers 1..L_d (shared by T0 and T1)."""
    if L_d < 1:
        raise ValueError("L_d must be at least 1")
    return rng.integers(1, int(L_d) + 1, size=n_sites)


def ends_t0(n_sites: int, S: int) -> np.ndarray:
    """T0: every span ends at the last stored step."""
    return np.full(n_sites, int(S) - 1, dtype=int)


def ends_t1(rng: np.random.Generator, n_sites: int, S: int, L_d: int) -> np.ndarray:
    """T1: each span's end falls uniformly over the last L_d stored steps."""
    if L_d > S:
        raise ValueError("L_d exceeds the stored length")
    return rng.integers(int(S) - int(L_d), int(S), size=n_sites)


def ends_window(rng: np.random.Generator, n_sites: int, S: int, W: int) -> np.ndarray:
    """C1 (per-type spec 6.1): each span's end falls uniformly over the last W stored steps.

    T1 is the case W = L_d. W is swept in units of the regional time scale and
    is not identified (D28, D30).
    """
    if isinstance(W, bool) or not isinstance(W, (int, np.integer)) or W < 1 or W > S:
        raise ValueError("window exceeds the stored length, is below 1, or is not an integer")
    return rng.integers(int(S) - int(W), int(S), size=n_sites)


def area_offsets(rng: np.random.Generator, areas: np.ndarray, L_R: int) -> np.ndarray:
    """S6: one offset per area, uniform on 0..L_R, expanded to a per-site vector."""
    areas = np.asarray(areas)
    uniq = sorted(set(areas.tolist()))
    off = {a: int(rng.integers(0, int(L_R) + 1)) for a in uniq}
    return np.array([off[a] for a in areas], dtype=int)


def ranks_for_ends(ends: np.ndarray, lengths: np.ndarray, S: int) -> np.ndarray:
    """The rank argument that makes sample_record place a span of length w ending at e.

    sample_record puts the end at (w - 1) + rint(rank (S - w)); the inverse is
    rank = (e - w + 1) / (S - w), with rank 0 when S == w. The realized ends
    are checked afterwards by sample_spans.
    """
    ends, lengths = np.asarray(ends, int), np.asarray(lengths, int)
    if (lengths < 1).any() or (ends < lengths - 1).any() or (ends > S - 1).any():
        raise ValueError("a span does not fit: need w >= 1 and w - 1 <= e <= S - 1")
    denom = S - lengths
    return np.where(denom > 0, (ends - lengths + 1) / np.where(denom > 0, denom, 1), 0.0)


def span_profile(record: np.ndarray, ends: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Each site's class profile averaged over its own span (end e, length w), the quantity sample_record draws from (R2 needs it)."""
    f = np.asarray(record, float); ends = np.asarray(ends, int); lengths = np.asarray(lengths, int)
    if (lengths < 1).any() or (ends - lengths + 1 < 0).any() or (ends >= len(f)).any():
        raise ValueError("spans must lie inside the stored record")
    p = np.array([f[e - w + 1:e + 1, i].mean(0) for i, (e, w) in enumerate(zip(ends, lengths))])
    return p / p.sum(1, keepdims=True)


def sample_spans(record: np.ndarray, ends: np.ndarray, lengths: np.ndarray, totals: np.ndarray,
                 rng: np.random.Generator) -> np.ndarray:
    """Counts for every site from its own span (end e, length w), through sample_record.

    Asserts that the realized slots equal the drawn ends (the audit's inversion
    check, 1,010,255 trials), so a rounding slip cannot pass silently.
    """
    S = len(record)
    ranks = ranks_for_ends(ends, lengths, S)
    lengths = np.asarray(lengths, int)
    realized = lengths - 1 + np.rint(ranks * (S - lengths)).astype(int)
    if not np.array_equal(realized, np.asarray(ends, int)):
        raise AssertionError("realized span ends differ from the drawn ends")
    return sample_record(record, ranks, totals, rng, window=lengths)


# --- Drift-extension timing (spec 2026-11-14, draft 3, sections 3.3 and 5) --------------------
# Three deterministic functions of their arguments. None holds or receives a random generator:
# the caller draws the base uniforms once per record, so every observation design of a record
# reads the same uniforms whatever the order in which designs are generated.

def _unit_uniforms(x, name: str) -> np.ndarray:
    """A 1-D array of finite values on [0, 1), or a ValueError naming the argument."""
    a = np.asarray(x, float)
    if a.ndim != 1 or a.size == 0 or not np.isfinite(a).all() or (a < 0).any() or (a >= 1).any():
        raise ValueError(f"{name} must be a non-empty 1-D array of finite values on [0, 1)")
    return a


def _is_int(x) -> bool:
    return isinstance(x, (int, np.integer)) and not isinstance(x, (bool, np.bool_))


def lengths_from_uniforms(v, L: int) -> np.ndarray:
    """Span lengths 1 + floor(v L): uniform integers on 1..L when v is uniform on [0, 1).

    The same v under two limits gives lengths in the same order, which is how
    the three span-length limits of a record are paired.
    """
    v = _unit_uniforms(v, "v")
    if not _is_int(L) or L < 1:
        raise ValueError("the span-length limit L must be an integer of at least 1")
    return 1 + np.floor(v * int(L)).astype(int)


def ends_coupled(ends, latitude, c: float, north_late: bool, u) -> np.ndarray:
    """Reassign a record's span ends to sites along latitude with coupling c.

    r is the latitude rank (1 = southernmost, n = northernmost), q = (r - 0.5) / n
    for a north-late design and 1 - (r - 0.5) / n for a south-late one, and the
    score is z = c q + (1 - c) u. The end values, sorted ascending, are assigned
    to the sites in ascending z (stable sort: equal scores keep site order).
    The multiset of ends is unchanged. Sign: under north_late the later ends go
    to the northern sites, so the rank correlation of latitude with the score
    is positive (+0.700 in expectation at c = 0.5 and 28 sites); under
    south-late it is negative. c = 0 orders by u alone, c = 1 by latitude alone.
    """
    ends = np.asarray(ends); lat = np.asarray(latitude, float)
    if ends.ndim != 1 or not np.issubdtype(ends.dtype, np.integer):
        raise ValueError("ends must be a 1-D integer array")
    u = _unit_uniforms(u, "u"); n = len(ends)
    if lat.shape != (n,) or u.shape != (n,):
        raise ValueError("ends, latitude and u must have one entry per site")
    if not np.isfinite(lat).all() or len(np.unique(lat)) != n:
        raise ValueError("latitudes must be finite and distinct: the latitude rank is undefined under ties")
    if not isinstance(north_late, (bool, np.bool_)):
        raise ValueError("north_late must be True (later in the north) or False (later in the south)")
    if isinstance(c, (bool, np.bool_)) or not np.isfinite(c) or c < 0 or c > 1:
        raise ValueError("the coupling c must lie in [0, 1]")
    r = np.empty(n); r[np.argsort(lat, kind="stable")] = np.arange(1, n + 1)
    q = (r - 0.5) / n
    if not north_late:
        q = 1.0 - q
    z = float(c) * q + (1.0 - float(c)) * u
    out = np.empty(n, dtype=ends.dtype)
    out[np.argsort(z, kind="stable")] = np.sort(ends, kind="stable")
    return out


def ends_step(w, labels, late_label, S: int, W: int) -> np.ndarray:
    """Span ends for two clusters that occupy disjoint halves of the window.

    The window is the last W stored steps, S - W .. S - 1. Its early half is
    the first floor(W / 2) of them and its late half the remaining
    W - floor(W / 2). A site's end is h0 + floor(w h), with h0 the first step
    of its cluster's half and h that half's length, so every end of the late
    cluster is later than every end of the early one.
    """
    w = _unit_uniforms(w, "w"); labels = np.asarray(labels)
    if labels.shape != w.shape:
        raise ValueError("w and labels must have one entry per site")
    uniq = np.unique(labels)
    if len(uniq) != 2:
        raise ValueError(f"a step design needs exactly two cluster labels, got {len(uniq)}")
    if late_label not in uniq.tolist():
        raise ValueError(f"late_label {late_label!r} is not one of the two cluster labels {uniq.tolist()}")
    if not _is_int(S) or not _is_int(W) or W < 2 or W > S:
        raise ValueError("W must be an integer from 2 (each half needs a step) to the stored length S")
    early = int(W) // 2; start = int(S) - int(W)
    late = labels == late_label
    h0 = np.where(late, start + early, start)
    h = np.where(late, int(W) - early, early)
    return (h0 + np.floor(w * h)).astype(int)
