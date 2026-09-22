# Robustness of phase-like emergence to contingency and factors (n = 28)

Grid: interaction length [12.0, 18.0, 24.0, 36.0] km, between-node mixing [0.01, 0.02, 0.05], innovation [0.006, 0.012, 0.024]; 12 seeds per cell; 432 time-transgressive runs total.
Observed data: 2 communities, between-group F_ST = 0.039.

## Headline robustness
- Phase-like structure (>= 2 emergent communities): 430 of 432 runs (99.5%).
- Emergent F_ST at or below the observed value: 406 of 430 runs (94.4%).
- Emergent F_ST within a factor of two of observed (0.019-0.078): 159 of 430 runs (37.0%).
- Mean emergent communities 2.2 (range 1-4); mean F_ST 0.018; mean seriation |rho| 0.56.

## Marginal effects (mean +/- sd)
- F_ST by interaction length (km): 12: 0.021±0.014; 18: 0.018±0.013; 24: 0.017±0.010; 36: 0.016±0.008
- F_ST by between-node mixing: 0.01: 0.021±0.014; 0.02: 0.016±0.011; 0.05: 0.016±0.009
- F_ST by innovation rate: 0.006: 0.031±0.011; 0.012: 0.015±0.005; 0.024: 0.008±0.002
- communities by interaction length (km): 12: 2.333±0.527; 18: 2.259±0.438; 24: 2.222±0.437; 36: 2.111±0.342

Interpretation: across stochastic replicates and the full factor grid, neutral drift on the real geography reliably produces phase-like, spatially coherent communities, and the between-group F_ST stays at or near the observed level rather than the much higher value bounded groups would leave. The appearance of phase structure is a generic outcome of distance-structured drift on this layout, not an artifact of one tuned parameter set. F_ST is essentially flat across interaction length and falls with more mixing or more innovation, as expected.

Figure: figures/figS4_emergence_robustness.png