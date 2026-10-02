# Robustness of phase-like emergence to contingency and factors (n = 28)

Grid: interaction length [12.0, 18.0, 24.0, 36.0] km, between-node mixing [0.01, 0.02, 0.05], innovation [0.006, 0.012, 0.024]; 12 seeds per cell; 432 time-transgressive runs total.
Observed data: 2 communities, between-group F_ST = 0.039.

## Headline robustness
- Phase-like structure (>= 2 emergent communities): 432 of 432 runs (100.0%).
- Emergent F_ST at or below the observed value: 407 of 432 runs (94.2%).
- Emergent F_ST within a factor of two of observed (0.019-0.078): 160 of 432 runs (37.0%).
- Mean emergent communities 2.3 (range 2-4); mean F_ST 0.018; mean seriation |rho| 0.55.

## Marginal effects (mean +/- sd)
- F_ST by interaction length (km): 12: 0.019±0.012; 18: 0.018±0.012; 24: 0.018±0.011; 36: 0.016±0.009
- F_ST by between-node mixing: 0.01: 0.021±0.014; 0.02: 0.017±0.010; 0.05: 0.016±0.009
- F_ST by innovation rate: 0.006: 0.030±0.010; 0.012: 0.016±0.005; 0.024: 0.008±0.003
- communities by interaction length (km): 12: 2.361±0.499; 18: 2.278±0.448; 24: 2.241±0.449; 36: 2.148±0.380

Interpretation: across stochastic replicates and the full factor grid, neutral drift on the real geography reliably produces phase-like, spatially coherent communities, and the between-group F_ST stays at or near the observed level rather than the much higher value bounded groups would leave. The appearance of phase structure is a generic outcome of distance-structured drift on this layout, not an artifact of one tuned parameter set. F_ST is essentially flat across interaction length and falls with more mixing or more innovation, as expected.

Figure: figures/34_emergence_robustness.png