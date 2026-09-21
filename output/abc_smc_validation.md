# ABC-SMC validation battery

Config: full, 20 worker process(es). SBC n=200 (particles 300, rounds 5); coverage grid [np.float64(-0.4), np.float64(-0.3), np.float64(-0.2), np.float64(-0.1), np.float64(0.0), np.float64(0.1), np.float64(0.2), np.float64(0.3), np.float64(0.4)] x 20 reps.

## 1. Simulation-based calibration

- b rank KS vs uniform: D = 0.088, p = 0.086 (p > 0.05 is consistent with calibration).
- Rank histograms: figures/abc_smc_sbc (should be flat).

## 2. Coverage / power on b

- true b = -0.40: 95% coverage 0.95, mean width 0.264.
- true b = -0.30: 95% coverage 0.90, mean width 0.307.
- true b = -0.20: 95% coverage 0.75, mean width 0.326.
- true b = -0.10: 95% coverage 0.80, mean width 0.347.
- true b = +0.00: 95% coverage 0.90, mean width 0.219.
- true b = +0.10: 95% coverage 0.65, mean width 0.246.
- true b = +0.20: 95% coverage 1.00, mean width 0.334.
- true b = +0.30: 95% coverage 0.90, mean width 0.349.
- true b = +0.40: 95% coverage 0.90, mean width 0.357.

## 3. Cross-check vs rejection ABC (analysis 19)

- ABC-SMC b = +0.001 [-0.026, +0.029].
- Rejection ABC b = +0.029 [-0.028, +0.065].
- Posterior-mean gap = 0.029 (should be within Monte-Carlo error).
- Overlay: figures/abc_smc_crosscheck.