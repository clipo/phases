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
