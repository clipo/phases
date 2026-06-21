# Robustness of phase-like emergence to contingency and factors (n = 29)

Grid: interaction length [12.0, 18.0, 24.0, 36.0] km, between-node mixing [0.01, 0.02, 0.05], innovation [0.006, 0.012, 0.024]; 12 seeds per cell; 432 time-transgressive runs total.
Observed data: 2 communities, between-group F_ST = 0.043.

## Headline robustness
- Phase-like structure (>= 2 emergent communities): 99.8% of runs.
- Emergent F_ST at or below the observed value: 85% of runs.
- Emergent F_ST within a factor of two of observed (0.022-0.086): 53% of runs.
- Mean emergent communities 2.2 (range 1-3); mean F_ST 0.028; mean seriation |rho| 0.44.

## Marginal effects (mean +/- sd)
- F_ST by interaction length (km): 12: 0.028±0.017; 18: 0.029±0.019; 24: 0.028±0.017; 36: 0.026±0.015
- F_ST by between-node mixing: 0.01: 0.034±0.022; 0.02: 0.026±0.015; 0.05: 0.022±0.010
- F_ST by innovation rate: 0.006: 0.043±0.018; 0.012: 0.025±0.009; 0.024: 0.014±0.005
- communities by interaction length (km): 12: 2.269±0.443; 18: 2.213±0.409; 24: 2.120±0.353; 36: 2.130±0.336

Interpretation: across stochastic replicates and the full factor grid, neutral drift on the real geography reliably produces phase-like, spatially coherent communities, and the between-group F_ST stays at or near the observed level rather than the much higher value bounded groups would leave. The appearance of phase structure is a generic outcome of distance-structured drift on this layout, not an artifact of one tuned parameter set. F_ST is essentially flat across interaction length and falls with more mixing or more innovation, as expected.

Figure: figures/figS5_emergence_robustness.png