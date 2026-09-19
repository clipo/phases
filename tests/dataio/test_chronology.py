import pytest
from pathlib import Path
from mls_emergence.dataio import chronology

LMV = Path("data/LMVData.xlsx")

@pytest.mark.data
@pytest.mark.skipif(not LMV.exists(), reason="LMVData not present (gitignored)")
def test_period_assignments_have_late_mississippian_column():
    df = chronology.load_period_assignments(LMV)
    assert {"Number", "Period", "Terminal Period"}.issubset(df.columns)
    assert len(df) > 2000  # zones 15 + 16 combined
