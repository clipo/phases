# Would a culture historian have found these phases in pottery made by drift?

Basin phase set, 28 assemblages, 3 phases; agreement scored on the 26 assemblages inside a Phillips (1970) phase area.
300 realizations of the calibrated pooled drift model (2000 learners, innovation 0.001, mixing 0.02), which contains no groups,
sampled at the observed sherd totals and clustered exactly as the real record is in `76_phase_recovery.py`
(chi-square composition), by two methods: k-means (median of 10 initializations) and Ward linkage.

## k-means

| weight on composition | real pottery | drift pottery, median (95 percent) | share of drift records at or above the real one |
|---|---|---|---|
| 0.0 | **0.480** | 0.480 (0.480 to 0.480) | 100% |
| 0.1 | **0.480** | 0.480 (0.480 to 0.863) | 100% |
| 0.25 | **0.740** | 0.480 (0.480 to 0.863) | 13% |
| 0.5 | **0.740** | 0.480 (0.480 to 0.863) | 13% |
| 1.0 | **0.740** | 0.510 (0.376 to 0.771) | 6% |
| 2.0 | **0.539** | 0.454 (0.308 to 0.682) | 21% |
| pottery alone | **0.358** | 0.123 (-0.042 to 0.386) | 5% |

## Ward linkage

| weight on composition | real pottery | drift pottery, median (95 percent) | share of drift records at or above the real one |
|---|---|---|---|
| 0.0 | **0.740** | 0.740 (0.740 to 0.740) | 100% |
| 0.1 | **0.740** | 0.740 (0.480 to 0.863) | 95% |
| 0.25 | **0.740** | 0.740 (0.454 to 0.863) | 68% |
| 0.5 | **0.740** | 0.628 (0.445 to 0.863) | 36% |
| 1.0 | **0.740** | 0.487 (0.348 to 0.771) | 6% |
| 2.0 | **0.454** | 0.439 (0.247 to 0.670) | 46% |
| pottery alone | **0.358** | 0.098 (-0.045 to 0.416) | 6% |

The weight-zero row is the site map alone and is the same in both columns by construction.

Read the rows where pottery enters. Where the real pottery sits inside what drift pottery gives,
the published phases are as recoverable from a record with no social groups in it as from the
real one, and the appearance of phases needs nothing beyond distance decay and where the sites are.
Where the real pottery sits above it, the record carries phase-aligned structure that drift on
this geography does not produce.

Figure written to fig15_phases_under_drift.png and its siblings.
