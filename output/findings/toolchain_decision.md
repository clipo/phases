# Toolchain decision: rstan builds on aarch64

Date: 2026-08-31. Produced by `analyses/00_setup/00_install_deps.R` and
verified by `analyses/00_setup/01_verify_toolchain.R`.
Plan: `docs/superpowers/plans/2026-08-31-stan-rebuild-program.md`, item 0 task 1.

**Status: settled. rstan, no fallback needed.**

## The risk, and why it was a real one

The program plan named the rstan build as its one platform risk. This machine is
**aarch64 (ARM64) Linux**, which is not the platform `../mataa` builds on, so
rstan and StanHeaders compile from source rather than installing from a binary.
The recorded fallback ladder was: retry at `Ncpus=1`, then cmdstanr, which
downloads a prebuilt CmdStan and avoids compiling the rstan R/C++ interface
entirely. Falling back would have cost `bridgesampling`, which is the stated
reason for choosing rstan over cmdstanr in the first place.

Neither rung was needed.

## What was measured

| | |
|---|---|
| Platform | `aarch64-unknown-linux-gnu` |
| R | 4.3.3 (2024-02-29) |
| Compiler | g++ 13.3.0 |
| rstan | **2.32.7**, built from source, `Ncpus = 4` |
| StanHeaders | **2.32.10** |
| bridgesampling | 1.2.1, installed |
| Smoke model compile | 32.7 s |

Verification against a conjugate target computed in R from the prior and the
likelihood, never read off the fit (rule 3). Prior `mu ~ Normal(0, 10)`,
likelihood `y ~ Normal(mu, 1)` on three points with `ybar` exactly 1.0:

| quantity | analytic | NUTS | difference |
|---|---|---|---|
| posterior mean | 0.996678 | 0.992636 | 4.04e-03 |
| posterior sd | 0.576390 | 0.570572 | 5.82e-03 |

Diagnostics, the full rule-16 set even for a smoke test: R-hat 1.0008, ESS 3307,
divergences 0 of 8000 (0.00%), treedepth saturation 0. The pass tolerance on the
mean is 10 times the Monte Carlo standard error, derived rather than chosen.

## The consequence that matters downstream

StanHeaders is on the **2.32** line, the same line `../mataa` reports
(they record Stan 2.32.2). That predates the built-in
`dirichlet_multinomial_lpmf`, so their reason for hand-writing the density and
naming it `dm_lpmf` rather than the built-in name holds here for the same
reason, not merely by inheritance. Task 2 carries their corrected sum-of-logs
form forward on that basis.

The two projects are not on identical toolchains (aarch64 against x86, 2.32.10
against 2.32.2). Any claim that they share a toolchain has to be stated at that
resolution.
