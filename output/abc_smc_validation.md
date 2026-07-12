# ABC-SMC validation battery

Config: full, 12 worker process(es). SBC n=200 (particles 300, rounds 5); coverage grid [np.float64(-0.4), np.float64(-0.3), np.float64(-0.2), np.float64(-0.1), np.float64(0.0), np.float64(0.1), np.float64(0.2), np.float64(0.3), np.float64(0.4)] x 20 reps.

## 1. Simulation-based calibration

- b rank KS vs uniform: D = 0.068, p = 0.305 (p > 0.05 is consistent with calibration).
- Rank histograms: figures/abc_smc_sbc (should be flat).

## 2. Coverage / power on b

- true b = -0.40: 95% coverage 0.90, mean width 0.257.
- true b = -0.30: 95% coverage 0.85, mean width 0.308.
- true b = -0.20: 95% coverage 0.85, mean width 0.321.
- true b = -0.10: 95% coverage 0.95, mean width 0.361.
- true b = +0.00: 95% coverage 0.90, mean width 0.241.
- true b = +0.10: 95% coverage 0.65, mean width 0.245.
- true b = +0.20: 95% coverage 1.00, mean width 0.332.
- true b = +0.30: 95% coverage 1.00, mean width 0.345.
- true b = +0.40: 95% coverage 0.90, mean width 0.359.

## 3. Cross-check vs rejection ABC (analysis 19)

- ABC-SMC b = +0.006 [-0.020, +0.035].
- Rejection ABC b = +0.028 [-0.021, +0.063].
- Posterior-mean gap = 0.022 (should be within Monte-Carlo error).
- Overlay: figures/abc_smc_crosscheck.