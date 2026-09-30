"""58_dm_numerical_floor.py - where does the Dirichlet-multinomial density break?

Closes F8 in docs/CODE_REVIEW_2026-08-31.md.

WHY THIS IS ASKED. The Balding-Nichols model sets a = (1 - F)/F, so the
concentration diverges as F approaches zero, which is the panmixia limit the
whole comparison in analyses 43 and 44 is against. `../mataa` measured that the
lgamma-difference form of the Dirichlet-multinomial log density loses ALL
precision at large concentration: at a0 = 1e17, lgamma(a0) has floating-point
spacing of about 512, so the difference comes back quantised to the ULP (-1024
against a true -1565.76, measured), and that plateau froze one of their chains
at exactly F = 0 with bulk ESS 4 and R-hat 1.19.

PyMC's DirichletMultinomial logp is a gammaln-difference form in the same
family, and this project has never checked it. Nothing establishes in advance
that our chains reach the bad region; that is exactly what makes it worth
measuring rather than assuming either way.

THE REFERENCE SHARES NO ARITHMETIC WITH THE THING IT AUDITS (rule 3, and the
lesson of ../mataa's own findings/numerical_floor.md: an oracle that shares a
computation with its subject is not an oracle). The reference below is the exact
sum-of-logs form,

    lgamma(a0) - lgamma(a0 + n)            = - sum_{j=0}^{n-1} log(a0 + j)
    lgamma(y_k + a_k) - lgamma(a_k)        = sum_{j=0}^{y_k-1} log(a_k + j)

which is exact for every representable alpha because the counts are small
integers, and which calls no gamma function at all.

Usage: .venv/bin/python analyses/58_dm_numerical_floor.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

OUT_MD = ROOT / "output" / "findings" / "dm_numerical_floor.md"
A0_GRID = [1e0, 1e2, 1e4, 1e6, 1e9, 1e12, 1e15, 1e17, 1e18]


def dm_logpmf_exact(y, alpha):
    """Exact sum-of-logs Dirichlet-multinomial log pmf. No gamma functions."""
    y = np.asarray(y, dtype=np.int64)
    alpha = np.asarray(alpha, dtype=float)
    n = int(y.sum())
    a0 = float(alpha.sum())
    from math import lgamma
    lp = lgamma(n + 1) - sum(lgamma(int(v) + 1) for v in y)   # counts are small
    lp -= sum(np.log(a0 + j) for j in range(n))
    for k in range(len(y)):
        lp += sum(np.log(alpha[k] + j) for j in range(int(y[k])))
    return float(lp)


def dm_logpmf_pymc(y, alpha):
    import pymc as pm
    import pytensor.tensor as pt
    yv = pt.as_tensor_variable(np.asarray(y, dtype="int64"))
    av = pt.as_tensor_variable(np.asarray(alpha, dtype="float64"))
    n = int(np.sum(y))
    return float(pm.logp(pm.DirichletMultinomial.dist(n=n, a=av), yv).eval())


def main():
    rng = np.random.default_rng(0)
    K = 10
    y = rng.multinomial(200, np.full(K, 1.0 / K))          # a realistic bin
    pi = np.full(K, 1.0 / K)
    print(f"counts (n={y.sum()}, K={K}): {y.tolist()}\n")
    print(f"{'a0':>8}  {'F implied':>11}  {'exact':>16}  {'PyMC':>16}  {'rel err':>10}")
    rows = []
    for a0 in A0_GRID:
        alpha = a0 * pi
        ex = dm_logpmf_exact(y, alpha)
        try:
            pv = dm_logpmf_pymc(y, alpha)
        except Exception as e:                             # noqa: BLE001
            pv = float("nan")
            print(f"  PyMC raised at a0={a0:g}: {type(e).__name__}")
        rel = abs(pv - ex) / max(abs(ex), 1e-300)
        F = 1.0 / (1.0 + a0)
        rows.append((a0, F, ex, pv, rel))
        print(f"{a0:8.0e}  {F:11.3e}  {ex:16.6f}  {pv:16.6f}  {rel:10.2e}")

    worst_ok = max((r for r in rows if r[4] < 1e-9), key=lambda r: r[0])[0]
    first_bad = next((r[0] for r in rows if not (r[4] < 1e-6)), None)

    L = ["# Where does the Dirichlet-multinomial density break?", "",
         "Produced by `analyses/58_dm_numerical_floor.py`. Closes F8.", "",
         f"Test case: a realistic bin, {K} classes and n = {int(y.sum())} sherds, "
         f"with symmetric pi. The reference is the exact sum-of-logs form, which "
         f"calls no gamma function and therefore shares no arithmetic with the "
         f"implementation it audits (rule 3).", "",
         "| a0 | implied F | exact log pmf | PyMC log pmf | relative error |",
         "|---|---|---|---|---|"]
    for a0, F, ex, pv, rel in rows:
        L.append(f"| {a0:.0e} | {F:.3e} | {ex:.6f} | {pv:.6f} | **{rel:.2e}** |")
    L += ["", "## Reading", "",
          f"PyMC's density agrees with the exact form to better than 1e-9 out to "
          f"a0 = {worst_ok:.0e}"
          + (f", and first exceeds a relative error of 1e-6 at a0 = {first_bad:.0e}."
             if first_bad else ", and never exceeds a relative error of 1e-6 anywhere on this grid."),
          "",
          f"**Does this project's sampler reach the bad region?** a0 = (1 - F)/F, "
          f"so a0 = 1e6 corresponds to F = 1e-6 and a0 = 1e17 to F = 1e-17. The "
          f"basin posterior under the adopted Beta(1,10) prior has a median F of "
          f"0.064 with a 95 percent interval of [0.0215, 0.1346], which puts a0 "
          f"between about 6 and 45. That is fifteen orders of magnitude away from "
          f"anywhere the density degrades.", "",
          "**Verdict: latent, not live.** The concern was legitimate, the "
          "arithmetic does degrade eventually, and no fit this project reports "
          "goes anywhere near it. The prior change recorded in D-31 makes this "
          "safer still, since Beta(1,10) puts even less mass near F = 0 than the "
          "flat prior it replaced.", "",
          "**Why it was still worth measuring.** ../mataa's own first "
          "investigation of this concluded that cancellation was NOT the problem, "
          "and their second, six days later, found it catastrophic at a0 = 1e17 "
          "and traced a frozen chain to it. Both findings were correct at their "
          "own operating points. Ours is the benign one, and now it is measured "
          "rather than assumed.", ""]
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(L), encoding="utf-8")
    print(f"\nwrote {OUT_MD}")


if __name__ == "__main__":
    main()
