# R and Stan environment

Recorded 2026-09-22 by `analyses/00_setup/02_record_environment.R`.

| | |
|---|---|
| platform | `aarch64-unknown-linux-gnu` |
| R | R version 4.3.3 (2024-02-29) |
| C++ | g++ 13.3.0 |
| rstan | 2.32.7 |
| StanHeaders | 2.32.10 |
| bridgesampling | 1.2.1 |
| posterior | 1.7.0 |
| loo | 2.10.1 |
| bayesplot | 1.16.0 |
| here | 1.0.2 |
| testthat | 3.3.2 |
| renv | 1.2.4 |

## Comparison with ../mataa

mataa builds on Linux Mint / x86 and records Stan 2.32.2. This box is
aarch64 with StanHeaders on the same 2.32 line. Same Stan line, different
platform. The consequence that matters: 2.32 predates the built-in
`dirichlet_multinomial_lpmf`, so the hand-written `dm_lpmf` is required here
for the same reason it is required there.

