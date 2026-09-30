"""66_innovation_boundary_highrep.py - the innovation-boundary sweep at 1,200 realizations.

Analysis 65 sweeps three departures at 250 realizations per cell. The
innovation-boundary departure (its departure C) is the one the manuscript
leans on, and its individual realizations are very widely spread, so the
manuscript quotes it at 1,200 realizations per cell. Until 2026-09-22 that run
was invoked by hand (`65_other_departures.py --only C --reps 1200 --tag c_highrep`) and had no
entry in the pipeline order, so a full rerun left it computed at a superseded
calibration cell (rule 15). This wrapper runs exactly that invocation so the
runner and MANIFEST carry it. Outputs, written by 65's own code:
`output/other_departures_c_highrep.json` and
`output/findings/other_departures_c_highrep.md`.

Usage: .venv/bin/python analyses/66_innovation_boundary_highrep.py
"""
import datetime
import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT_MD = ROOT / "output" / "findings" / "other_departures_c_highrep.md"
OUT_JSON = ROOT / "output" / "other_departures_c_highrep.json"


def stamp():
    """Append the provenance line that lets MANIFEST credit this entry point.

    65 performs the writes; this wrapper only fixes the arguments. The manifest
    builder credits a script with an output only when the script itself writes
    it, so the wrapper adds one dated line to the report it caused.
    """
    if not (OUT_MD.exists() and OUT_JSON.exists()):
        raise RuntimeError("65 --only C --reps 1200 did not leave its two outputs")
    OUT_MD.write_text(OUT_MD.read_text(encoding="utf-8").rstrip("\n")
                      + f"\n\nRun through `analyses/66_innovation_boundary_highrep.py` "
                        f"(65 --only C --reps 1200 --tag c_highrep) on {datetime.date.today().isoformat()}; "
                        f"the JSON beside it is `output/other_departures_c_highrep.json`.\n",
                      encoding="utf-8")


if __name__ == "__main__":
    # --tag routes the outputs to the *_c_highrep names; without it, 65
    # overwrites the 250-realization full-grid files (it did so once, 2026-09-22).
    sys.argv = [str(HERE / "65_other_departures.py"), "--only", "C", "--reps", "1200", "--tag", "c_highrep"]
    runpy.run_path(str(HERE / "65_other_departures.py"), run_name="__main__")
    stamp()
