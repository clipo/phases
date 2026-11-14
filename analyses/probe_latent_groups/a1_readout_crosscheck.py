"""Stan generated quantities against the numpy readouts on the same draws (rule 3).

Reads the 20 draws validate_a1.R exports (delta, pi, resp, z_draw and Stan's
Delta, Delta_w, coherence, coherence_hard) and recomputes each with
mls_emergence.probe.readouts. Any disagreement above 1e-6 fails.
Usage: .venv/bin/python analyses/probe_latent_groups/a1_readout_crosscheck.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from mls_emergence.probe import readouts as ro

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "output" / "probe_latent_groups" / "a1_crosscheck_draws.json"


def main() -> int:
    d = json.loads(SRC.read_text())
    edges = np.asarray(d["edges"], int) - 1
    worst = {"Delta": 0.0, "Delta_w": 0.0, "coherence": 0.0, "coherence_hard": 0.0}
    for s in d["draws"]:
        delta = np.asarray(s["delta"], float); pi = np.asarray(s["pi"], float); resp = np.asarray(s["resp"], float)
        z = np.asarray(s["z_draw"], int) - 1
        q = ro.comembership(resp)
        mine = {"Delta": ro.delta_contrast(delta, pi), "Delta_w": ro.delta_contrast(delta, resp.mean(axis=0)),
                "coherence": ro.coherence(q, edges), "coherence_hard": ro.coherence(ro.hard_comembership(z), edges)}
        for k, v in mine.items():
            worst[k] = max(worst[k], abs(v - float(s[k])))
    ok = all(v < 1e-6 for v in worst.values())
    print("largest |Stan - numpy| over", len(d["draws"]), "draws:", {k: f"{v:.2e}" for k, v in worst.items()}, "| PASS" if ok else "| FAIL")
    out = ROOT / "output" / "findings" / "a1_readout_crosscheck.md"
    out.write_text("# Stan generated quantities against numpy readouts\n\n"
                   f"Produced by `analyses/probe_latent_groups/a1_readout_crosscheck.py` on {len(d['draws'])} exported draws of the A1 null fit. "
                   "Largest absolute differences: " + ", ".join(f"{k} {v:.2e}" for k, v in worst.items()) + f". {'PASS' if ok else 'FAIL'} at tolerance 1e-6.\n")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
