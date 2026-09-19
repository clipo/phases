# ABC-SMC validation battery

Config: full, 20 worker process(es). SBC n=200 (particles 300, rounds 5); coverage grid [np.float64(-0.4), np.float64(-0.3), np.float64(-0.2), np.float64(-0.1), np.float64(0.0), np.float64(0.1), np.float64(0.2), np.float64(0.3), np.float64(0.4)] x 20 reps.

## 1. Simulation-based calibration

- b rank KS vs uniform: D = 0.070, p = 0.264 (p > 0.05 is consistent with calibration).
- Rank histograms: figures/abc_smc_sbc (should be flat).

## 2. Coverage / power on b

- true b = -0.40: 95% coverage 0.95, mean width 0.264.
- true b = -0.30: 95% coverage 0.95, mean width 0.308.
- true b = -0.20: 95% coverage 0.75, mean width 0.323.
- true b = -0.10: 95% coverage 0.80, mean width 0.365.
- true b = +0.00: 95% coverage 0.85, mean width 0.246.
- true b = +0.10: 95% coverage 0.80, mean width 0.238.
- true b = +0.20: 95% coverage 1.00, mean width 0.331.
- true b = +0.30: 95% coverage 0.95, mean width 0.345.
- true b = +0.40: 95% coverage 0.85, mean width 0.359.

## 3. Cross-check vs rejection ABC (analysis 19)

- SKIPPED: run analyses/19 and analyses/38 first to produce both posterior files.