from pathlib import Path
from mls_emergence.dataio import pfg

PFG = Path("data/raw/PFGData.xlsx")

def test_counts_shape_and_index():
    df = pfg.load_pfg_counts(PFG)
    assert df.shape[0] == 266          # assemblages (SherdData has 266 non-empty rows; xlsx max_row=393 is a blank-formatting artifact)
    assert "Parkin Punctated" in df.columns
    assert "Sherd Total" not in df.columns   # dropped
    assert "Site Name" not in df.columns     # moved to attrs
    assert df.index.name == "Site Number"
    assert (df.sum(axis=1) > 0).all()        # no empty assemblages
