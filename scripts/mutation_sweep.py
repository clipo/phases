#!/usr/bin/env python3
"""Mutation sweep over the signature estimators.

Rule 3 says a test must discriminate: it needs a probe point where the nearest
rival account predicts differently. The only way to know whether a suite
actually does that is to break the code on purpose and check the suite notices.

Each mutation below is semantically meaningful rather than random: a dropped
bias correction, a swapped denominator, an off-by-one in a summation limit, a
constant scale factor. These are the errors an estimator inherits silently,
which is why they are the ones worth probing.

A mutant that SURVIVES means the suite cannot tell the real estimator from a
wrong one. When this was first run (2026-09-04) seven of ten survived, including
three on cultural_fst, which is the paper's core estimator, and the dropped
(n-1) correction in Neiman's homozygosity. The tests at the time pinned only
endpoints: F_ST is zero for identical groups and large for disjoint ones, both
of which every mutant reproduced. Nothing pinned the arithmetic in between.
Value tests were added and all ten are now killed.

Run this after touching any estimator, and after adding tests that are supposed
to constrain one.

Usage: .venv/bin/python scripts/mutation_sweep.py
Exit 1 if any mutant survives.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
PY_BIN = str(ROOT / ".venv" / "bin" / "python")

# (module under src/mls_emergence, test target, original, mutated, label)
MUTANTS: list[tuple[str, str, str, str, str]] = [
    ("signatures/variance.py", "tests/signatures/test_variance.py",
     "return float(1.0 - np.sum(p ** 2))",
     "return float(1.0 - np.sum(p ** 3))",
     "gini_simpson: p^2 -> p^3"),
    ("signatures/variance.py", "tests/signatures/test_variance.py",
     "return float((H_T - H_S) / H_T)",
     "return float((H_T - H_S) / H_S)",
     "cultural_fst: divide by H_S not H_T"),
    ("signatures/variance.py", "tests/signatures/test_variance.py",
     "H_S = float(np.average(H_within, weights=sizes))",
     "H_S = float(np.mean(H_within))",
     "cultural_fst: drop the size weighting"),
    ("signatures/convergence.py", "tests/signatures/test_convergence.py",
     "return float((x @ (y - y.mean())) / (x @ x))",
     "return float((x @ (y - y.mean())) / (x @ x)) * 1.5",
     "time_derivative: scale the slope by 1.5"),
    ("signatures/convergence.py", "tests/signatures/test_convergence.py",
     "return z.mean(axis=1)",
     "return z.mean(axis=1) * 0.5",
     "convergence_score: halve"),
    ("signatures/seriation.py", "tests/signatures/test_seriation.py",
     "return max(0, n_changes - 1)",
     "return max(0, n_changes)",
     "_col_violations: off-by-one"),
    ("signatures/assortativity.py", "tests/signatures/test_assortativity.py",
     "return float(200.0 - np.sum(np.abs(pa - pb)))",
     "return float(200.0 - np.sum((pa - pb) ** 2))",
     "brainerd_robinson: absolute -> squared difference"),
    ("signatures/neutral.py", "tests/signatures/test_neutral.py",
     "return float(np.sum(c * (c - 1)) / (N * (N - 1)))",
     "return float(np.sum(c * c) / (N * N))",
     "homozygosity_f: drop the (n-1) bias correction"),
    ("signatures/neutral.py", "tests/signatures/test_neutral.py",
     "return (1.0 - F) / F",
     "return (1.0 - F) / (F * 2.0)",
     "theta_f: halve"),
    ("signatures/neutral.py", "tests/signatures/test_neutral.py",
     "return float(np.sum(theta / (theta + i)))",
     "return float(np.sum(theta / (theta + i + 1)))",
     "_ewens_expected_k: shift the summation index"),
]


def main() -> int:
    survivors: list[str] = []
    skipped: list[str] = []
    for mod, tests, old, new, label in MUTANTS:
        path = ROOT / "src" / "mls_emergence" / mod
        original = path.read_text()
        if old not in original:
            # The anchor is gone, which means the estimator was refactored and
            # this mutant no longer describes a real error. That is a reason to
            # update the mutant, not to pass silently.
            skipped.append(label)
            print(f"  SKIP      {label}  (anchor not found; update this mutant)")
            continue
        path.write_text(original.replace(old, new, 1))
        try:
            r = subprocess.run([PY_BIN, "-m", "pytest", "-q", "-x", tests],
                               cwd=ROOT, capture_output=True, text=True,
                               timeout=900)
            rc = r.returncode
        finally:
            # Restore unconditionally, including on timeout or interrupt: a
            # mutation harness that can leave a mutant in the tree is worse
            # than no harness.
            path.write_text(original)
        if rc == 0:
            survivors.append(label)
        print(f"  {'SURVIVED' if rc == 0 else 'killed  '}  {label}")

    print(f"\n{len(MUTANTS)} mutants, {len(survivors)} survived, "
          f"{len(skipped)} skipped")
    if survivors:
        print("\nA surviving mutant means the suite cannot tell the estimator "
              "from a wrong one:")
        for s in survivors:
            print(f"  - {s}")
    return 1 if (survivors or skipped) else 0


if __name__ == "__main__":
    sys.exit(main())
