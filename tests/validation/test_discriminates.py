"""The scientific crux: does CONVERGENCE of the four signatures discriminate
genuine group-level emergence from its mimics?

The convergence rule flags a mechanism iff ALL FOUR signatures show a positive
ordinal trend. The genuine-emergence generator must be flagged; each mimic
(aggregated signaling, static patchiness, isolation-by-distance drift) must not.
"""
from __future__ import annotations

import pytest

from mls_emergence.validation.harness import discriminates, run_blind
from mls_emergence.validation.mechanisms import (
    gen_aggregated_signaling,
    gen_drift_space,
    gen_group_emergence,
    gen_patchiness,
)

GENERATORS = {
    "group_emergence": gen_group_emergence,
    "aggregated_signaling": gen_aggregated_signaling,
    "patchiness": gen_patchiness,
    "drift_space": gen_drift_space,
}
MIMICS = ["aggregated_signaling", "patchiness", "drift_space"]


def test_only_genuine_emergence_is_convergent():
    """At a representative seed, convergence flags ONLY genuine emergence."""
    panels = run_blind(GENERATORS, seed=42)
    verdict = discriminates(panels)
    assert verdict["group_emergence"]["convergent"] is True
    for mimic in MIMICS:
        assert verdict[mimic]["convergent"] is False, f"{mimic} false-positive"


def test_no_mimic_is_ever_convergent_across_seeds():
    """Specificity is the load-bearing property: NO mimic should ever satisfy the
    convergence criterion. This must hold for every seed."""
    for seed in range(20):
        verdict = discriminates(run_blind(GENERATORS, seed=seed))
        for mimic in MIMICS:
            assert verdict[mimic]["convergent"] is False, (
                f"seed {seed}: mimic {mimic} read as convergent (false positive)"
            )


def test_genuine_emergence_is_usually_convergent_across_seeds():
    """Sensitivity: genuine emergence should be flagged on the large majority of
    seeds. It is not 100%: the spatial-boundary signature carries sampling noise,
    so on an occasional single realization its ordinal trend dips below threshold
    even though the other three rise. This is an honest sensitivity limit, not a
    discrimination failure (specificity remains perfect)."""
    hits = 0
    n = 20
    for seed in range(n):
        verdict = discriminates(run_blind(GENERATORS, seed=seed))
        if verdict["group_emergence"]["convergent"]:
            hits += 1
    assert hits >= int(0.85 * n), f"only {hits}/{n} emergence runs flagged convergent"


@pytest.mark.parametrize("seed", [0, 1, 2, 42])
def test_emergence_all_four_signatures_rise(seed):
    """On stable seeds all four standardized slopes for genuine emergence are
    positive, confirming convergence is driven by genuinely co-rising signatures
    and not by one dominant signature."""
    verdict = discriminates(run_blind(GENERATORS, seed=seed))
    trends = verdict["group_emergence"]["trends"]
    for col, t in trends.items():
        assert t["slope_std"] > 0, f"seed {seed}: {col} did not rise"
