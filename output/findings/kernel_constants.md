# The inherited Matern range constants are wrong for kernels 2 and 3

Date: 2026-08-31. Found by `tests/testthat/test-kernels.R`, written under
Verification Regime rule 19 (a collaborator's committed result is re-derived
before anything here depends on it).

**Status: corrected in this repository. Affects `../mataa`, which is where the
code came from, and where it may be load-bearing.**

## What is wrong

`stan/spatial/kernels.stanfunctions` was copied verbatim from `../mataa`
(branch `stan-rebuild`) on 2026-08-31. Its header states that the practical
range is "the distance at which correlation falls to 0.05,
`range_const(kernel) * rho`", with constants 3, 4.743865 and 5.918649, and that
the half-distance `half_const(kernel) * rho` is where correlation falls to 0.5.

Those constants for the Matern 3/2 and 5/2 kernels are the roots in the
kernel's **argument** `x`, where `x = sqrt(3) d / rho` or `sqrt(5) d / rho`.
Using them as multipliers on `rho` omits the sqrt factor that converts argument
to distance. Kernel 1 escapes because its argument *is* `d / rho`.

Measured, by evaluating the kernels themselves at the quoted distances:

| constant | claimed correlation | actual correlation | error |
|---|---|---|---|
| `range_const(2) = 4.743865` | 0.05 | **0.00249** | 20x too far out |
| `range_const(3) = 5.918649` | 0.05 | **0.00013** | 385x too far out |
| `half_const(2) = 1.678347` | 0.50 | **0.2135** | |
| `half_const(3) = 2.3302562` | 0.50 | **0.0833** | |
| `range_const(1) = 3` | 0.0498 | 0.049787 | correct (documented rounding) |
| `half_const(1) = 0.6931472` | 0.50 | 0.500000 | correct |

Correct values are the same roots divided by the sqrt factor:

| kernel | quantity | inherited | correct | inflation |
|---|---|---|---|---|
| Matern 3/2 | practical range | 4.743865 | **2.738871** | 1.73x |
| Matern 5/2 | practical range | 5.918649 | **2.646900** | 2.24x |
| Matern 3/2 | half distance | 1.678347 | **0.968994** | 1.73x |
| Matern 5/2 | half distance | 2.330256 | **1.042122** | 2.24x |

## Consequence here: none yet, and that is luck rather than care

Every fit in this project so far uses `kernel = 1`, the exponential, where the
constants are right. No number of ours is affected. But `kernel` is a data
switch on `composition_gp.stan`, deliberately, so the first sensitivity run
that flips it to Matern 3/2 or 5/2 would have inflated every reported range and
half-distance by 73 or 124 percent, and nothing in the pipeline would have
objected.

## Consequence for ../mataa: potentially real, and theirs to assess

`../mataa`'s own header calls d50 "the PRIMARY reported scale" and states that
"every kernel inflates" the 0.05 tail "while d50 is recovered by the right
kernel". Their `docs/superpowers/specs/2026-08-28-kernel-choice-by-process-simulation.md`
makes smoothness a measured choice, so Matern 3/2 and 5/2 fits are not
hypothetical there. Any d50 or practical range they report under those two
kernels is inflated by 1.73x or 2.24x. Their spatial findings
(`output/02_spatial/findings/ahu_range_recovery.md`, `ahu_kernel_comparison.md`,
`03_process/findings/radius_to_d50_map.md`) are where this would land.

This is flagged, not fixed: that repository is not ours to edit, and the
assessment of which reported numbers are affected belongs to whoever owns it.

## Why the inherited test did not catch it

The header says the constants are "pinned by test against
`R/simulate_gp.R::range_const()`, which must agree to 1e-6". If that R function
computes the constant the same way, the test compares an implementation with a
copy of itself. That is the exact failure `../mataa`'s own
`findings/numerical_floor.md` names: "an oracle that shares a
[computation] with the thing it audits is not an oracle".

`tests/testthat/test-kernels.R` avoids it in two ways. It solves the roots from
the kernel definitions with `uniroot`, sharing no code with either project. And
it asserts the **consequence** rather than the value: for each kernel it
evaluates the correlation at `range_const * rho` and requires it to be 0.05, and
at `half_const * rho` and requires 0.5. That assertion is what a constant is
*for*, and it fails loudly on the inherited values.
