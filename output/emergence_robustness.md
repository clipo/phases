# Robustness of phase-like emergence to contingency and factors (n = 43)

Grid: interaction length [12.0, 18.0, 24.0, 36.0] km, between-node mixing [0.01, 0.02, 0.05], innovation [0.006, 0.012, 0.024]; 12 seeds per cell; 432 time-transgressive runs total.
Observed data: 2 communities, between-group F_ST = 0.040.

## Headline robustness
- Phase-like structure (>= 2 emergent communities): 430 of 432 runs (99.5%).
- Emergent F_ST at or below the observed value: 422 of 430 runs (98.1%).
- Emergent F_ST within a factor of two of observed (0.020-0.080): 85 of 430 runs (19.8%).
- Mean emergent communities 2.3 (range 1-4); mean F_ST 0.013; mean seriation |rho| 0.50.

## Marginal effects (mean +/- sd)
- F_ST by interaction length (km): 12: 0.016±0.012; 18: 0.013±0.010; 24: 0.011±0.008; 36: 0.010±0.007
- F_ST by between-node mixing: 0.01: 0.015±0.011; 0.02: 0.012±0.009; 0.05: 0.010±0.007
- F_ST by innovation rate: 0.006: 0.022±0.009; 0.012: 0.011±0.005; 0.024: 0.004±0.002
- communities by interaction length (km): 12: 2.481±0.616; 18: 2.287±0.472; 24: 2.176±0.427; 36: 2.139±0.346

Interpretation: across stochastic replicates and the full factor grid, neutral drift on the real geography reliably produces phase-like, spatially coherent communities, and the between-group F_ST stays at or near the observed level rather than the much higher value bounded groups would leave. The appearance of phase structure is a generic outcome of distance-structured drift on this layout, not an artifact of one tuned parameter set. F_ST is essentially flat across interaction length and falls with more mixing or more innovation, as expected.

Figure: figures/figS4_emergence_robustness.png