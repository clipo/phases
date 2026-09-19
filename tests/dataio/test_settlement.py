import pytest
from pathlib import Path
from mls_emergence.dataio import settlement, pfg

LMV = Path("data/LMVData.xlsx"); PFG = Path("data/raw/PFGData.xlsx")

@pytest.mark.data
@pytest.mark.skipif(not LMV.exists(), reason="LMVData gitignored/absent")
def test_lmv_has_coordinates_and_area():
    df = settlement.load_lmv(LMV)
    for col in ["Number", "Area", "Northing", "Easting", "Zone"]:
        assert col in df.columns
    assert len(df) > 2000

@pytest.mark.data
@pytest.mark.skipif(not LMV.exists(), reason="LMVData gitignored/absent")
def test_join_reports_coverage():
    counts = pfg.load_pfg_counts(PFG)
    joined, unmatched = settlement.join_pfg_to_lmv(counts, settlement.load_lmv(LMV))
    assert "Easting" in joined.columns and "Northing" in joined.columns
    assert isinstance(unmatched, list)


def test_normalize_grid_zero_o_substitution():
    from mls_emergence.dataio import settlement as S
    assert S.normalize_grid("15-0-10") == "15-O-10"
    assert S.normalize_grid("10-P-1") == "10-P-1"           # unchanged
    assert S.normalize_grid("12-N-3/A&B") == "12-N-3"        # compound -> base
    assert S.normalize_grid(" 13-n-4/B,C,D,E ") == "13-N-4"  # case+space+compound


@pytest.mark.data
@pytest.mark.skipif(not LMV.exists(), reason="LMVData gitignored/absent")
def test_join_coverage_improves_after_reconciliation():
    counts = pfg.load_pfg_counts(PFG)
    joined, unmatched = settlement.join_pfg_to_lmv(counts, settlement.load_lmv(LMV))
    matched = len(counts) - len(unmatched)
    assert matched >= 245   # up from 233; aim to recover most of the 33
