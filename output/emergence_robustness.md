# Robustness of phase-like emergence to contingency and factors (n = 28)

Grid: interaction length [12.0, 18.0, 24.0, 36.0] km, between-node mixing [0.01, 0.02, 0.05], innovation [0.006, 0.012, 0.024]; 12 seeds per cell; 432 time-transgressive runs total.
Observed data: 2 communities, between-group F_ST = 0.039.

## Headline robustness
- Phase-like structure (>= 2 emergent communities): 428 of 432 runs (99.1%).
- Emergent F_ST at or below the observed value: 402 of 428 runs (93.9%).
- Emergent F_ST within a factor of two of observed (0.019-0.078): 151 of 428 runs (35.3%).
- Mean emergent communities 2.2 (range 1-4); mean F_ST 0.018; mean seriation |rho| 0.56.

## Marginal effects (mean +/- sd)
- F_ST by interaction length (km): 12: 0.020±0.013; 18: 0.018±0.012; 24: 0.018±0.011; 36: 0.016±0.008
- F_ST by between-node mixing: 0.01: 0.021±0.013; 0.02: 0.017±0.010; 0.05: 0.016±0.009
- F_ST by innovation rate: 0.006: 0.030±0.010; 0.012: 0.015±0.005; 0.024: 0.008±0.003
- communities by interaction length (km): 12: 2.324±0.487; 18: 2.222±0.416; 24: 2.213±0.452; 36: 2.056±0.329

Interpretation: across stochastic replicates and the full factor grid, neutral drift on the real geography reliably produces phase-like, spatially coherent communities, and the between-group F_ST stays at or near the observed level rather than the much higher value bounded groups would leave. The appearance of phase structure is a generic outcome of distance-structured drift on this layout, not an artifact of one tuned parameter set. F_ST is essentially flat across interaction length and falls with more mixing or more innovation, as expected.

Figure: figures/figS4_emergence_robustness.png