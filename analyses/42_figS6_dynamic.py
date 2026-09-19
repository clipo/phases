"""42_figS6_dynamic.py - assemble the merged dynamic-sufficiency figure (Figure S6).

Three panels, built from artifacts saved by earlier scripts (so it runs after the
ABC-SMC and SBC scripts, which `run_all.sh` reaches only at 38-39):

  A  ABC-SMC posterior of the transmission bias b (from 38_abc_smc_transmission.py,
     output/abc_smc_posterior.npz): whole-sequence, early-half, late-half.
  B  Simulation-based calibration of the ABC-SMC posterior: the rank histogram of
     the true b within its posterior (from 39_abc_smc_validation.py,
     output/abc_smc_sbc_ranks.npz). A flat histogram is calibration.
  C  Tempo and mode as parameters rather than model selection: the posteriors
     of the mean-reversion rate and the directional drift (from
     54_tempo_mode_posterior.py, output/tempo_posterior.npz). Replaced the
     Akaike weights over four models on 2026-09-02; unbiased motion, stasis and
     mean reversion are limits of one model, so two parameters say what a
     four-way comparison was approximating.

Writes figures/figS5_dynamic.{png,pdf,svg,tif} directly (no manual renaming).

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
from figstyle import OI_BLUE, OI_ORANGE, OI_VERMIL, panel_label, save  # noqa: E402

OUT = ROOT / "output"


# Which script writes each input, so a stale product can be named and re-made.
PRODUCER = {
    "abc_smc_posterior.npz": "38_abc_smc_transmission.py",
    "abc_smc_sbc_ranks.npz": "39_abc_smc_validation.py",
    "tempo_posterior.npz": "54_tempo_mode_posterior.py",
}


def _need(name):
    """Load an upstream product, refusing one older than the script that writes it.

    This figure plots saved posteriors rather than refitting, so a product left
    behind by an edited generator is invisible: the figure regenerates cleanly
    and prints numbers nobody can reproduce. That happened here. The committed
    figS6 showed mu -0.092 [-0.175, -0.015] while re-running 54 gave
    -0.093 [-0.180, -0.007], and the manuscript had been quoting the older fit.
    The same class of defect as the stale closure_posterior.npz that
    53.load_posterior now guards against.
    """
    p = OUT / name
    if not p.exists():
        raise SystemExit(
            f"missing {p}; run the upstream script first "
            "(38/39 for the ABC-SMC posterior and SBC ranks, 54 for the tempo posteriors)."
        )
    producer = ROOT / "analyses" / PRODUCER[name]
    if producer.exists() and p.stat().st_mtime < producer.stat().st_mtime - 1:
        raise SystemExit(
            f"{p.name} is older than {producer.name}, which writes it. Re-run\n"
            f"    .venv/bin/python analyses/{producer.name}\n"
            "before regenerating this figure; plotting the stale product would "
            "print numbers the current code does not produce."
        )
    return np.load(p, allow_pickle=True)


def main():
    post = _need("abc_smc_posterior.npz")
    sbc = _need("abc_smc_sbc_ranks.npz")
    # tempo_posterior.npz replaces tempo_akaike.npz: parameters, not weights.
    tempo = _need("tempo_posterior.npz")

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(7, 2.6))
    panel_label(axA, "A")
    panel_label(axB, "B")
    panel_label(axC, "C")

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

    # Panel C: tempo and mode as PARAMETERS, not Akaike weights over four models.
    # Rule 18 replaced the model selection with the posteriors of the two
    # parameters that distinguish the four cases, because unbiased motion,
    # stasis and mean reversion are limits of one model
    # (docs/FREQUENTIST_INVENTORY.md item E). The caption changed on 2026-09-02
    # and this panel is drawn from analyses/54_tempo_mode_posterior.py's saved
    # posteriors so the two cannot drift apart.
    alpha, mu = tempo["alpha"], tempo["mu"]
    axC2 = axC.twiny()
    axC.hist(np.log10(alpha), bins=40, color=OI_BLUE, alpha=0.65, density=True)
    axC.axvline(np.log10(np.median(alpha)), color=OI_BLUE, ls=":", lw=1.0)
    axC.set_xlabel("log$_{10}$ mean-reversion rate $\\alpha$", fontsize=7, color=OI_BLUE)
    axC.set_ylabel("posterior density")
    axC.set_yticks([])
    axC2.hist(mu, bins=40, color=OI_ORANGE, alpha=0.55, density=True)
    axC2.axvline(0.0, color="0.4", lw=0.8)
    axC2.set_xlabel("directional drift $\\mu$ per bin", fontsize=7, color=OI_ORANGE)
    # The verdict is computed, not typed: mu is resolved when its 95 percent
    # interval excludes zero; alpha is unresolved when its interval spans more
    # than two orders of magnitude (the unbiased-walk and stasis limits).
    _mlo, _mhi = np.percentile(mu, [2.5, 97.5])
    _alo, _ahi = np.percentile(alpha, [2.5, 97.5])
    verdict = ("$\\mu$ " + ("resolved" if (_mlo > 0 or _mhi < 0) else "unresolved")
               + "; $\\alpha$ " + ("unresolved" if _ahi / max(_alo, 1e-12) > 100 else "resolved"))
    axC.text(0.02, 0.96,
             f"$\\alpha$ {np.median(alpha):.2f} "
             f"[{np.percentile(alpha, 2.5):.2f}, {np.percentile(alpha, 97.5):.1f}]\n"
             f"$\\mu$ {np.median(mu):+.3f} "
             f"[{np.percentile(mu, 2.5):+.3f}, {np.percentile(mu, 97.5):+.3f}]\n"
             + verdict,
             transform=axC.transAxes, va="top", fontsize=5.8)

    save(fig, "figS5_dynamic")
    print("wrote figures/figS5_dynamic.{png,pdf,svg,tif}")


if __name__ == "__main__":
    main()
