# The NAD27 conversion: ledger and measured effect

Produced by `scripts/measure_datum.py`. Tolerance 10 m (PFG UTMs are rounded to 10 m). Basin set of 28 assemblages from `basin_members_curated.txt`.

## Assemblage coordinates (`mainfort-pfg-cplXY.txt`, after the published corrections)

| class | all assemblages | basin members | meaning |
|---|---|---|---|
| converted | 20 | 16 | file coordinate is PFG's UTM read as NAD83; moved to the NAD27 reading |
| near | 11 | 8 | within 100 m of that reading but not within tolerance; derivation not demonstrated; left as recorded |
| other | 3 | 2 | more than 100 m from PFG's UTM either way; left as recorded |
| no_pfg_utm | 21 | 2 | no PFG UTM to compare against; left as recorded |

### Basin members, row by row

| assemblage | site | class | m from NAD83 reading | m from NAD27 reading | m from Mainfort's own XY |
|---|---|---|---|---|---|
| Cramor_Place | 12-O-5 | converted | 0 | 210 | - |
| Walls | 13-P-1 | converted | 0 | 210 | 0 |
| Barton_Ranch | 11-O-10 | converted | 0 | 210 | - |
| Hollywood | 13-O-10 | converted | 0 | 210 | 0 |
| Castile_Landing | 13-N-21 | converted | 0 | 210 | 0 |
| Cummins | 11-O-4 | converted | 0 | 210 | 0 |
| Fortune | 11-N-15 | converted | 0 | 210 | 0 |
| Irby | 13-P-10 | converted | 0 | 210 | 0 |
| Neeleys_Ferry | 11-N-4 | converted | 0 | 210 | 0 |
| Nickel | 13-N-15/B | converted | 0 | 210 | 0 |
| Commerce | 13-O-11 | converted | 0 | 210 | 0 |
| Lake_Cormorant | 13-P-8 | converted | 0 | 210 | 0 |
| Turnbow | 11-N-12 | converted | 0 | 210 | 0 |
| Rose_Mound | 12-N-3/A&B | converted | 0 | 210 | 0 |
| Davis | 13-N-5 | converted | 0 | 210 | 1926 |
| Starkley | 13-N-16 | converted | 0 | 210 | 2459 |
| Clay_Hill | 13-N-7 | near | 14 | 224 | 0 |
| Mound_Place | 12-P-1 | near | 19 | 220 | 0 |
| Woodlyn | 13-P-11 | near | 24 | 186 | 0 |
| Vernon_Paul | 11-N-9 | near | 38 | 209 | - |
| Parkin | 11-N-1 | near | 41 | 226 | 0 |
| Kent_Place | 13-N-4/B,C,D,E | near | 48 | 258 | - |
| Williamson | 11-N-13 | near | 53 | 158 | - |
| Big_Eddy | 12-N-4 | near | 54 | 181 | 0 |
| Grant | nan | no_pfg_utm | - | - | 0 |
| Holden_Lake | nan | no_pfg_utm | - | - | - |
| Beck | 13-O-7 | other | 639 | 781 | 21719 |
| Belle_Meade | 13-O-5 | other | 5433 | 5546 | 4747 |

Of the 11 near cases, 8 coincide with Mainfort's own coordinate for the site (`mainfortXY` sheet) within tolerance; those are his points, whose datum is not recorded, and they stand as recorded under the 2026-09-19 ruling on non-PFG coordinates.

## Settlement table (`LMVData.xlsx`)

367 rows whose Source begins with "PFG", or that the settlement correction table replaced with PFG UTMs, converted from EPSG:267zz to EPSG:269zz; all others as recorded. For the basin's assemblage sites the compilation's UTMs are identical to PFG's site table digit for digit, which is why the Source label is trusted.

## Effect on the basin set

Assemblages moved: 16 of 28, by 209 to 209 m.

| k | same partition? | plug-in F_ST before | after |
|---|---|---|---|
| 2 | yes | 0.006240 | 0.006240 |
| 3 | yes | 0.015590 | 0.015590 |
| 4 | NO | 0.019438 | 0.018650 |
| 5 | NO | 0.021422 | 0.020634 |
| 6 | yes | 0.032704 | 0.032704 |

**The partition or the F_ST changed at some k; every position-dependent result must be rerun and re-read.**
