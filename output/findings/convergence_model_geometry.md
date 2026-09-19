# The hierarchical convergence model's default target_accept is too loose

Date: 2026-08-31. Measured with `.venv/bin/python` on the rebuilt environment
(pymc 6.3.1, pytensor 3.3.0, arviz 1.3.0, numpy 2.4.6).
Surfaced by `tests/inference/test_convergence_model.py::test_non_centered_model_samples_without_divergences`
failing on a clean install.

**Status: real, bounded, and it does NOT touch the manuscript's number.**

## What was measured

`sample_convergence` on the test's own synthetic panel (three signatures with a
true slope of +1.0, seriation slope +1.0, noise 0.15, seed 5), 500 draws,
1000 tune, 2 chains, sweeping `target_accept`:

| target_accept | divergences | tau bulk ESS | tau tail ESS | tau R-hat | treedepth >= 10 |
|---|---|---|---|---|---|
| **0.9 (the library default)** | **7 / 1000 (0.70%)** | 208 | 180 | **1.0112** | 0 |
| 0.95 | 1 / 1000 (0.10%) | 237 | 292 | 0.9999 | 0 |
| 0.99 | 0 / 1000 (0.00%) | 283 | 396 | 1.0076 | 0 |
| 0.999 | 0 / 1000 (0.00%) | 319 | 418 | 1.0072 | 0 |

At the default the model produces divergences and an R-hat above the 1.01
threshold **on synthetic data with a strong, clean, correctly-specified
signal**. That is the easiest case the model will ever see.

## Scope: the reported result is not affected

`analyses/40_hierarchical_convergence.py:108` passes `target_accept=0.99`
explicitly, so the fit behind the manuscript's P(all four signatures rise) ran
in the clean region of this table. The defect is in the library default that
every other caller inherits, including the test.

## Why this is a defect and not a tuning preference

Three things, in increasing order of importance.

1. **Rule 4, no silent parameters.** The default 0.9 in
   `inference/convergence_model.py:sample_convergence` is contradicted by the
   only production call site, which overrides it. A default that every real
   caller must override is not a default.
2. **The test encodes the symptom as acceptable.** It is named
   `test_non_centered_model_samples_without_divergences` and asserts
   `ndiv <= 5`. A test that tolerates five divergences is not testing what its
   name says, and it is why this went unnoticed: it passed until a version bump
   moved 5 to 7.
3. **Rule 17.** The model is already non-centered in `z_b`, which is the
   standard funnel fix, and it still funnels in `tau` at 0.9. Raising
   `adapt_delta` is a legitimate response only when the geometry is sound and
   the step size is simply too coarse. The evidence here says it is: divergences
   go to exactly zero by 0.99 and ESS improves monotonically, which is the
   signature of a fixable step size rather than a pathological posterior. Had
   divergences persisted at 0.999, the correct response would have been to
   reparameterize, not to keep climbing.

## Disposition

Open. The fix is to make 0.99 the default with a comment pointing here, tighten
the test to assert zero divergences plus an R-hat threshold, and add the
rule-16 diagnostic set to `convergence_summary` so a caller cannot report this
model's posterior without seeing its divergence count. Scheduled as item 1.

## A note on the environment

These versions are newer than those the original analysis ran under (the repo
pins only lower bounds). The divergence count moved from at-or-below 5 to 7,
which is what made the test fail. This is a reminder that the reproduction of
the frozen PyMC baseline in item 0 task 1 step 4 is comparing across a library
version change as well as across implementations, and any discrepancy has to be
attributed before it is interpreted.
