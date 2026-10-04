"""The Dirichlet-multinomial density degrades at large concentration (F8).

Pins two things at once: that PyMC's density is exact where this project
actually samples, and that it is NOT exact far out, so the test cannot quietly
stop discriminating if a future PyMC changes the implementation.
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
for p in (str(ROOT), str(ROOT / "analyses")):
    if p not in sys.path:
        sys.path.insert(0, p)


@pytest.fixture(scope="module")
def mod():
    return importlib.import_module("58_dm_numerical_floor")


COUNTS = np.array([21, 17, 13, 13, 27, 28, 21, 22, 19, 19])
PI = np.full(10, 0.1)


@pytest.mark.parametrize("a0", [1e0, 1e2, 1e4, 1e6])
def test_density_is_exact_where_this_project_samples(mod, a0):
    """a0 = (1 - F)/F. The basin posterior sits at F ~ 0.064, i.e. a0 ~ 15, and
    its 95% interval spans a0 ~ 6 to 45. This range must be exact."""
    ex = mod.dm_logpmf_exact(COUNTS, a0 * PI)
    pv = mod.dm_logpmf_pymc(COUNTS, a0 * PI)
    assert abs(pv - ex) / abs(ex) < 1e-9, f"a0={a0:g}: {pv} vs exact {ex}"


def test_density_does_degrade_far_out_so_this_test_discriminates(mod):
    """The probe. If PyMC were exact everywhere, the assertions above would pass
    on any implementation and prove nothing about the region that matters.

    Measured 2026-09-02: at a0 = 1e17 PyMC returns about +192 where the exact
    value is about -26.25. A POSITIVE log pmf is impossible for a discrete
    distribution, which is how visible the failure is once you look."""
    a0 = 1e17
    ex = mod.dm_logpmf_exact(COUNTS, a0 * PI)
    pv = mod.dm_logpmf_pymc(COUNTS, a0 * PI)
    assert ex < 0, "the exact log pmf must be negative"
    assert abs(pv - ex) / abs(ex) > 1.0, (
        "PyMC's density no longer degrades at a0=1e17. That is good news, but it "
        "means the exactness assertions above are no longer discriminating and "
        "this test needs a new probe point.")


def test_exact_reference_shares_no_gamma_arithmetic(mod):
    """Rule 3, and ../mataa's numerical_floor lesson: an oracle that shares a
    computation with its subject is not an oracle. The reference must reproduce
    a case computable by hand."""
    # Two classes, one sherd each, alpha = (1, 1): DM reduces to a known value.
    # P(y=(1,1) | n=2, alpha=(1,1)) = 2!/(1!1!) * B(2,2)/B(1,1)... computed here
    # from first principles rather than through the implementation under test.
    from math import log
    exact_by_hand = log(2) + log(1) + log(1) - (log(2) + log(3))
    got = mod.dm_logpmf_exact(np.array([1, 1]), np.array([1.0, 1.0]))
    assert abs(got - exact_by_hand) < 1e-12
