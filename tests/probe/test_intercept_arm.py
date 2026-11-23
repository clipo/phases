"""The intercept-prior arm of the validation record writers (a1_contracts.py, a3_contracts.py; 2026-11-23)."""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "analyses", ROOT / "analyses" / "probe_latent_groups"):
    sys.path.insert(0, str(_p))

BASE10 = {"N": 28, "y": [[1]], "_provenance": {}, "prior_alpha_sd": 10.0, "prior_scale_sd": 2.0}


@pytest.mark.parametrize("module", ["a1_contracts", "a3_contracts"])
def test_each_arm_has_its_own_records_directory_and_pinned_values(module):
    m = importlib.import_module(module)
    main, alt = m.intercept_arm(10), m.intercept_arm(2)
    assert main["out"] == ROOT / "data" / "stan" and main["pin"].name == "a0_pinned_nuisances.json"
    assert alt["out"] == ROOT / "data" / "stan" / "alphaN0-2" and alt["pin"].name == "a0_pinned_nuisances_alphaN0-2.json"
    # A shared directory or file would let one arm's records overwrite the other's.
    assert main["out"] != alt["out"] and main["pin"] != alt["pin"]
    with pytest.raises(ValueError):
        m.intercept_arm(5)       # only the two priors the paper reports are accepted


@pytest.mark.parametrize("module", ["a1_contracts", "a3_contracts"])
def test_records_carry_the_arms_prior_and_refuse_values_pinned_under_the_other(module):
    m = importlib.import_module(module)
    alt, main = m.intercept_arm(2), m.intercept_arm(10)
    f = m.arm_fields(BASE10, {"prior_alpha_sd": 2}, alt)
    # The basin contract says 10; a record of the Normal(0, 2) arm that kept 10 would be fitted under the wrong prior.
    assert f["prior_alpha_sd"] == 2.0 and "y" not in f and "_provenance" not in f and f["prior_scale_sd"] == 2.0
    assert m.arm_fields(BASE10, {"prior_alpha_sd": 10}, main)["prior_alpha_sd"] == 10.0
    with pytest.raises(ValueError):
        m.arm_fields(BASE10, {"prior_alpha_sd": 10}, alt)      # values pinned from the Normal(0, 10) basin fit
    with pytest.raises(ValueError):
        m.arm_fields(BASE10, {"prior_alpha_sd": 2}, main)      # values pinned from the Normal(0, 2) basin fit
    with pytest.raises(ValueError):
        m.arm_fields(BASE10, {}, alt)                          # pinned values that do not state their prior
    with pytest.raises(ValueError):
        m.arm_fields(dict(BASE10, prior_alpha_sd=2.0), {"prior_alpha_sd": 10}, main)   # a basin contract not at the paper's prior
