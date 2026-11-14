"""Statement forms of the type-origination run (spec 2026-11-14, draft 1, section 7; Codex build audit 2026-11-14, F4).

The form is read on R, all seven rare ranks inside their bands. The
five-of-seven threshold R5 was chosen after the drift-extension run and must
not decide the form. Each table below is built by hand, one row per cell.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analyses" / "per_type"))
import type_origination as to


def _sets(rows):
    s = pd.DataFrame(rows, columns=["L", "n_rare", "K", "S3"])
    s["R"] = s["n_rare"] == 7; s["R5"] = s["n_rare"] >= 5
    return s


def test_five_of_seven_does_not_complete_the_intersection():
    # One cell has level, exactly five rare ranks, richness and space. No cell has level with all seven.
    # Spec rule (R): L and R is empty, so form a. The rival rule (R5) would report form c.
    s = _sets([(True, 5, True, True), (False, 7, True, True), (True, 0, False, False)])
    form, text = to.statement(s)
    assert form == "a" and "At no tested origination intensity" in text
    assert "at 1 cells" in to.r5_summary(s) and "1 have both" in to.r5_summary(s)   # R5 is still reported, separately


def test_form_b_when_all_seven_are_covered_without_richness_or_space():
    # Two cells have level and all seven; one fails richness, the other fails both. A third cell has
    # level, five rare ranks, richness and space, which under the rival rule (R5) would give form c.
    s = _sets([(True, 7, False, True), (True, 7, False, False), (True, 5, True, True)])
    form, text = to.statement(s)
    assert form == "b" and "at 2 cells" in text and "richness is outside at 2 of them and space at 1" in text


def test_form_c_needs_all_seven_with_richness_and_space():
    s = _sets([(True, 7, True, True), (True, 6, True, True)])
    form, text = to.statement(s)
    assert form == "c" and "at 1 cells" in text                    # the six-of-seven cell is not counted
