"""A3 under the inherited and the corrected intercept prior: a descriptive comparison (D25 follow-up).

Reads the six A3 findings tables written by validate_a3.R and
a3_pinned_range.R (free range, range pinned at 6 km and at 12 km; the
Normal(0, 2) copies kept as *_alphaN0-2.md, the Normal(0, 10) runs under the
plain names) and reports, per run: gate passes and rungs, and the posterior
medians of omega pooled by record type, gate-passing fits only. The tables
are parsed rather than the JSON summaries because the Normal(0, 2) JSON was
overwritten by the rerun; only its markdown copies remain.

It also reports the between-site contrast SD of the field that generated the
own-generator and step-against-slope records, sigma * sqrt(c(rho)) with
c(rho) = (N - sum_ij R_ij / N) / (N - 1), against A0's posterior median of
the same quantity under each intercept prior (a0_intercept_prior.json). The
generator centers the field across sites (a3_contracts.simulate_a3), so this
contrast SD is the field's whole contribution to the records.

Descriptive summaries only; no inference rests on them. Exploratory design
validation, not evidence about the record.
Usage: .venv/bin/python analyses/probe_latent_groups/a3_prior_comparison.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FINDINGS = ROOT / "output" / "findings"
OUTD = ROOT / "output" / "probe_latent_groups"
OUT = FINDINGS / "a3_prior_comparison.md"
RUNS = (("a3_validation", "free range"), ("a3_pinned_rho6", "range pinned at 6 km"),
        ("a3_pinned_rho12", "range pinned at 12 km"))
PRIORS = (("_alphaN0-2", "Normal(0, 2), inherited"), ("", "Normal(0, 10), corrected"))
GROUPS = (("own-generator, omega 0", r"own_(G1|L)_om0_r"), ("own-generator, omega 0.25", r"own_(G1|L)_om0\.25_r"),
          ("own-generator, omega 0.5", r"own_(G1|L)_om0\.5_r"), ("own-generator, omega 1", r"own_(G1|L)_om1_r"),
          ("smooth slope, rho 30 km, no barrier", r"slope_rho30_r"), ("smooth slope, rho 60 km, no barrier", r"slope_rho60_r"),
          ("simulator S0, no barrier", r"sim_S0_r"), ("simulator S1 and S2, barrier", r"sim_S[12]_r"))
CELL = re.compile(r"^([\d.]+) \[([\d.]+), ([\d.]+)\]$")


def read_table(path: Path) -> dict[str, dict]:
    """Record rows of one findings table: omega median and 90 percent interval, gate, rung."""
    if not path.exists():
        raise FileNotFoundError(f"{path} missing; run validate_a3.R / a3_pinned_range.R first")
    header, rows = None, {}
    for line in path.read_text().splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cells and cells[0] == "record":
            header = cells
        elif header and len(cells) == len(header) and re.match(r"(own|slope|sim)_", cells[0]):
            r = dict(zip(header, cells))
            m = CELL.match(r["omega post"])
            if m is None:
                raise ValueError(f"{path.name}: unparsed omega cell {r['omega post']!r} for {cells[0]}")
            rows[cells[0]] = {"median": float(m[1]), "q5": float(m[2]), "q95": float(m[3]),
                              "gate": r["gate"] == "pass", "rung": int(r["rung"])}
    if header is None or not rows:
        raise ValueError(f"{path.name}: no record table found")
    return rows


def contrast_sd(D: np.ndarray, sigma: float, rho: float) -> float:
    """SD of an exponential-kernel field's deviations from its site mean (contrast scale of D25)."""
    N = len(D); R = np.exp(-D / rho)
    return float(sigma * np.sqrt((N - R.sum() / N) / (N - 1)))


def main() -> int:
    L = ["# A3 under the inherited and the corrected intercept prior", "",
         ("Produced by `analyses/probe_latent_groups/a3_prior_comparison.py` from the findings tables of "
         "`validate_a3.R` and `a3_pinned_range.R`. Posterior medians of omega (SD in logits of a region-level shift), "
         "pooled by record type over gate-passing fits only. Descriptive; exploratory design validation, not evidence "
         "about the record. Prior on omega half-Normal(0, 2): median 1.35."), ""]
    for stem, run in RUNS:
        L += [f"## {run}", "", "| records | " + " | ".join(p for _, p in PRIORS) + " |", "|---|---|---|"]
        tabs = {tag: read_table(FINDINGS / f"{stem}{tag}.md") for tag, _ in PRIORS}
        gate = []
        for tag, _ in PRIORS:
            t = tabs[tag]; fails = sorted(k for k, v in t.items() if not v["gate"])
            gate.append(f"{sum(v['gate'] for v in t.values())} of {len(t)} pass "
                        f"(rung 1 in {sum(v['rung'] == 1 for v in t.values())})" + (f"; fail: {', '.join(fails)}" if fails else ""))
        L.append("| gate | " + " | ".join(gate) + " |")
        for label, pat in GROUPS:
            cells = []
            for tag, _ in PRIORS:
                v = sorted(x["median"] for k, x in tabs[tag].items() if re.match(pat, k) and x["gate"])
                if not v:
                    raise ValueError(f"{stem}{tag}: no gate-passing records for {label}")
                cells.append(f"mean {np.mean(v):.2f}, range {v[0]:.2f} to {v[-1]:.2f} (n {len(v)})")
            L.append(f"| {label} | " + " | ".join(cells) + " |")
        L.append("")

    D = np.asarray(json.loads((ROOT / "data" / "stan" / "basin_composition.json").read_text())["D"], float)
    pin = json.loads((OUTD / "a0_pinned_nuisances.json").read_text())
    truth = json.loads((ROOT / "data" / "stan" / "a3_truth.json").read_text())
    a0 = json.loads((OUTD / "a0_intercept_prior.json").read_text())
    L += ["## Contrast scale of the generating field", "",
          "| field | sigma | rho (km) | between-site contrast SD (logits) |", "|---|---|---|---|",
          (f"| own-generator records (A0 point under Normal(0, 2), `a0_pinned_nuisances.json`) | {pin['sigma']:.2f} | "
           f"{pin['rho']:.1f} | {contrast_sd(D, pin['sigma'], pin['rho']):.2f} |")]
    for rr in sorted({v["rho"] for k, v in truth.items() if k.startswith("slope_")}):
        s = {v["sigma"] for k, v in truth.items() if k.startswith("slope_") and v["rho"] == rr}
        if len(s) != 1:
            raise ValueError(f"slope records at rho {rr} do not share one sigma: {s}")
        s = s.pop(); L.append(f"| step-against-slope records | {s:.2f} | {rr:.1f} | {contrast_sd(D, s, rr):.2f} |")
    for v in a0.values():
        L.append(f"| A0 posterior median, {v['label']} | {v['sigma'][0]:.2f} | {v['rho'][0]:.1f} | "
                 f"{np.sqrt(v['field_contrast_var'][0]):.2f} (sqrt of the median of sigma^2 c(rho)) |")
    L.append("")
    OUT.write_text("\n".join(L) + "\n"); print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
