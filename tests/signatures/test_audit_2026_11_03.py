"""Pins for two statements the 2026-11-03 Codex audit corrected (F11, F17)."""
from __future__ import annotations

import importlib
import inspect
import sys
from pathlib import Path

import numpy as np

from mls_emergence.dataio.basin import load_curated
from mls_emergence.signatures.variance import cultural_fst

ROOT = Path(__file__).resolve().parents[2]


def test_two_sample_plugin_fst_leading_order_is_one_over_2n_not_one_over_n():
    # Expected value derived by hand, not through cultural_fst: E[H_T - H_S] = H/(2n) and
    # E[H_T] = H(1 - 1/(2n)), so the ratio of expectations is about 1/(2n). The rival, 1/n, is
    # the statement the supplement used to make; at n = 400 the two differ by a factor of two,
    # far outside the tolerance, so this fails if the supplement's old value were true.
    p, n = np.array([0.5, 0.3, 0.15, 0.05]), 400
    rng = np.random.default_rng(20261103)
    vals = [cultural_fst(np.vstack([rng.multinomial(n, p), rng.multinomial(n, p)])) for _ in range(20000)]
    assert abs(np.mean(vals) - 1 / (2 * n)) < 0.05 / (2 * n)
    assert np.mean(vals) < 0.75 / n


def test_a3_producer_selects_the_line_from_phase_labels_where_the_rival_rule_differs():
    # On the real sites both rules give one line, so only a fixture where they differ can show which
    # rule the producer uses. Labels place Parkin on the left half; the k-means rival would follow
    # the vertical split instead, so selecting from k-means returns the other line and fails here.
    for p in (ROOT / "analyses", ROOT / "analyses" / "probe_latent_groups"):
        sys.path.insert(0, str(p))
    a3 = importlib.import_module("a3_contracts")
    left, top = np.array([0, 0, 0, 0, 1, 1, 1, 1]), np.array([0, 1, 0, 1, 0, 1, 0, 1])
    by_x, by_y = {"regions": left + 1, "name": "x"}, {"regions": top + 1, "name": "y"}
    labels = np.where(left == 0, "Parkin", "Kent")
    assert a3.choose_line(labels, [by_y, by_x])["name"] == "x"
    assert a3.pl.best_line_for(top + 1, [by_y, by_x])["name"] == "y"  # the rival's answer differs
    assert "choose_line(labels, cands)" in inspect.getsource(a3.main)  # main uses this selection


def test_a3_line_equals_the_kmeans_best_line_on_the_real_sites():
    # The Methods say the same line also best matches the k-means split. A different line under the
    # two rules would make that clause false. This does not show which rule the producer uses; the
    # test above does.
    for p in (ROOT / "analyses", ROOT / "analyses" / "probe_latent_groups"):
        sys.path.insert(0, str(p))
    a3 = importlib.import_module("a3_contracts")
    m36 = importlib.import_module("36_canonical_phase_map")
    counts, coords = load_curated(); xy = coords.to_numpy(float)
    labels, _ = m36.assign_primary_phases([str(i) for i in counts.index], xy); labels = np.asarray(labels)
    cands = a3.pl.candidate_lines(xy); g1 = a3.pl.kmeans2_regions(xy)
    assert len(cands) == 30
    parkin = np.where(labels == "Parkin", 1, 2)
    assert np.array_equal(a3.pl.best_line_for(parkin, cands)["regions"], a3.pl.best_line_for(g1, cands)["regions"])
    # k-means leaves one Parkin-area site with the other group: eleven of twelve, not all.
    in_parkin_cluster = g1 == g1[list(counts.index).index("Parkin")]
    assert (labels == "Parkin").sum() == 12 and (in_parkin_cluster & (labels == "Parkin")).sum() == 11
    assert in_parkin_cluster.sum() == 11
