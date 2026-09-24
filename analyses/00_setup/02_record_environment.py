"""Record the Python environment. Companion to 02_record_environment.R.

The rebuilt stack is newer than the one the original analyses ran under
(pyproject pins lower bounds only), so any comparison against a previously
committed number is a comparison across a version change as well. Writing the
versions down is what makes that attributable (Verification Regime rules 1, 14).

Usage: .venv/bin/python analyses/00_setup/02_record_environment.py
"""
from __future__ import annotations

import datetime as _dt
import platform
import subprocess
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "output" / "findings" / "environment_python.md"

KEY = ["pymc", "pytensor", "arviz", "numpy", "scipy", "pandas", "matplotlib",
       "networkx", "geopandas", "shapely", "pyproj", "cartopy", "openpyxl",
       "xlrd", "h5py", "h5netcdf", "pytest"]


def main() -> None:
    rows = []
    for p in KEY:
        try:
            rows.append((p, version(p)))
        except PackageNotFoundError:
            rows.append((p, "NOT INSTALLED"))
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                            capture_output=True, text=True).stdout.strip()

    L = ["# Python environment", "",
         f"Recorded {_dt.date.today().isoformat()} by "
         "`analyses/00_setup/02_record_environment.py`.", "",
         "| | |", "|---|---|",
         f"| platform | `{platform.platform()}` |",
         f"| machine | `{platform.machine()}` |",
         f"| python | {sys.version.split()[0]} |"]
    L += [f"| {p} | {v} |" for p, v in rows]
    L += ["", "## Why this is recorded", "",
          "`pyproject.toml` pins lower bounds only, so a rebuild does not "
          "reproduce the original environment. Any discrepancy against a "
          "previously committed number has to be attributed to a version "
          "change or to a defect before it is interpreted. Measured instance: "
          "the hierarchical convergence model's divergence count at the "
          "library default moved from at-or-below 5 to 7 across this rebuild "
          "(`output/findings/convergence_model_geometry.md`), while the "
          "Balding-Nichols fit reproduced its committed values exactly.", "",
          "<details><summary>Full pip freeze</summary>", "",
          "```", freeze, "```", "", "</details>", ""]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
