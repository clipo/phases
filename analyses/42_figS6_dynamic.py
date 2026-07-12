"""42_figS6_dynamic.py - assemble the merged dynamic-sufficiency figure (Figure S6).

Three panels, built from artifacts saved by earlier scripts (so it runs after the
ABC-SMC and SBC scripts, which `run_all.sh` reaches only at 38-39):

  A  ABC-SMC posterior of the transmission bias b (from 38_abc_smc_transmission.py,
     output/abc_smc_posterior.npz): whole-sequence, early-half, late-half.
  B  Simulation-based calibration of the ABC-SMC posterior: the rank histogram of
     the true b within its posterior (from 39_abc_smc_validation.py,
     output/abc_smc_sbc_ranks.npz). A flat histogram is calibration.
  C  Tempo-and-mode Akaike weights (from 20_tempo_mode_ews.py,
     output/tempo_akaike.npz).

Writes figures/figS6_dynamic.{png,pdf,svg,tif} directly (no manual renaming).

Usage: .venv/bin/python analyses/42_figS6_dynamic.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "analyses"))

import os  # noqa: E402
os.environ.setdefault("MLS_FIG_COLOR", "1")  # supplement figure is online-only; render in color

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from figstyle import OI_BLUE, OI_ORANGE, OI_VERMIL, save  # noqa: E402

OUT = ROOT / "output"


def _need(name):
    p = OUT / name
    if not p.exists():
        raise SystemExit(
            f"missing {p}; run the upstream script first "
            "(38/39 for the ABC-SMC posterior and SBC ranks, 20 for the tempo weights)."
        )
    return np.load(p, allow_pickle=True)


def main():
    post = _need("abc_smc_posterior.npz")
    sbc = _need("abc_smc_sbc_ranks.npz")
    tempo = _need("tempo_akaike.npz")

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(7, 2.6))

    # Panel A: ABC-SMC posterior of b
    axA.hist(post["post_b"], bins=25, density=True, color=OI_BLUE, alpha=0.6, label="whole")
    axA.hist(post["early"], bins=25, density=True, histtype="step",
             color=OI_ORANGE, lw=1.3, label="early half")
    axA.hist(post["late"], bins=25, density=True, histtype="step",
             color=OI_VERMIL, lw=1.3, label="late half")
    axA.axvline(0, color="0.3", ls="--", lw=1.0)
    axA.set_xlabel("transmission bias b (0 = neutral)")
    axA.set_ylabel("posterior density")
    axA.legend(frameon=False, fontsize=6)

    # Panel B: SBC rank histogram for b (flat = calibrated)
    ranks = sbc["ranks"]
    b_ranks = ranks[:, 1]
    n = len(b_ranks)
    axB.hist(b_ranks, bins=10, range=(0, 1), color=OI_BLUE)
    axB.axhline(n / 10.0, color="0.3", ls="--", lw=1.0)  # uniform expectation
    axB.set_xlabel("SBC rank of true b")
    axB.set_ylabel("count")
    axB.set_yticks([])

    # Panel C: tempo-and-mode Akaike weights
    models = [str(m) for m in tempo["models"]]
    type_means = tempo["type_means"]
    div_w = tempo["div_w"]
    xpos = np.arange(len(models))
    axC.bar(xpos - 0.2, type_means, width=0.4, color=OI_BLUE, label="type trajectories (mean)")
    axC.bar(xpos + 0.2, div_w, width=0.4, color=OI_ORANGE, label="diversity trajectory")
    axC.set_xticks(xpos)
    axC.set_xticklabels(models, fontsize=7)
    axC.set_ylabel("Akaike weight")
    axC.set_ylim(0, 1)
    axC.legend(frameon=False, fontsize=6)

    save(fig, "figS6_dynamic")
    print("wrote figures/figS6_dynamic.{png,pdf,svg,tif}")


if __name__ == "__main__":
    main()
