# The primary phase tests under Phillips's (1970) phase areas

Produced by `analyses/95_phillips_assignment_impact.py`. Same 28 assemblages, same counts, same
coordinates, same seeds; the two columns differ only in the phase labels. The superseded column
is Mainfort (1996: Figure 1) with unmapped assemblages given the territory they fall in, which the
pipeline used until 2026-10-01; the Phillips column, the pipeline's since then, is the outline each
assemblage lies in (`data/raw/phillips1970_phase_outlines.csv`, Lipo 2001: Figure 2.6), with an
assemblage inside no outline given the nearest.

## The assignments

| scheme | Kent | Parkin | Walls | assigned directly | by territory or nearness |
|---|---|---|---|---|---|
| Mainfort 1996 (superseded) | 8 | 11 | 9 | 25 | 3 |
| Phillips 1970 (pipeline) | 11 | 12 | 5 | 26 | 2 |

5 assemblages change phase: Beck (Walls to Kent); Belle_Meade (Walls to Kent); Castile_Landing (Kent to Parkin); Commerce (Walls to Kent); Hollywood (Walls to Kent).
Assigned by nearness under Phillips: Cramor_Place (Kent, 1.0 km outside the line); Holden_Lake (Parkin, 1.0 km outside the line).

## Do the phase lines separate the pottery better than other divisions? (2000 posterior draws, 300 alternatives of each kind)

| three phases | Mainfort 1996 (superseded) | Phillips 1970 (pipeline) |
|---|---|---|
| plug-in F_ST | 0.0110 | 0.0123 |
| F_ST under the scheme, posterior median (95% CrI) | 0.0109 (0.0091 to 0.0128) | 0.0121 (0.0101 to 0.0144) |
| median F_ST, same sizes around random centers | 0.0144 (0.0115 to 0.0179) | 0.0112 (0.0091 to 0.0144) |
| median F_ST, same sizes made compact | 0.0156 (0.0112 to 0.0218) | 0.0103 (0.0079 to 0.0127) |
| P(scheme separates better than random centers) | 0.21 | 0.57 |
| P(scheme separates better than compact) | 0.05 | 0.72 |
| boundary excess at the lines, posterior median (95% CrI) | +1.0 (-1.6 to +3.1) | +7.3 (+4.6 to +9.8) |
| P(larger boundary excess than random centers) | 0.15 | 0.58 |
| P(larger boundary excess than compact) | 0.06 | 0.65 |

| Parkin against the rest | Mainfort 1996 (superseded) | Phillips 1970 (pipeline) |
|---|---|---|
| plug-in F_ST | 0.0062 | 0.0067 |
| F_ST, posterior median (95% CrI) | 0.0061 (0.0051 to 0.0074) | 0.0066 (0.0052 to 0.0083) |
| median F_ST, same sizes around random centers | 0.0064 (0.0052 to 0.0093) | 0.0062 (0.0049 to 0.0081) |
| median F_ST, same sizes made compact | 0.0076 (0.0053 to 0.0111) | 0.0057 (0.0042 to 0.0071) |
| P(separates better than random centers) | 0.28 | 0.52 |
| P(separates better than compact) | 0.00 | 0.94 |
| boundary excess, posterior median (95% CrI) | +6.1 (+3.3 to +8.9) | +5.9 (+3.5 to +8.3) |
| P(larger boundary excess than random centers) | 0.28 | 0.43 |
| P(larger boundary excess than compact) | 0.00 | 0.49 |

## What recovers the scheme? (adjusted Rand index over the directly assigned assemblages, k-means, 40 starts, chi-square composition; median, with the range over starts)

| weight on composition | Mainfort 1996 (superseded) | Phillips 1970 (pipeline) |
|---|---|---|
| site map alone | 1.000 (1.000 to 1.000) | 0.480 (0.480 to 0.480) |
| 0.1 | 1.000 (1.000 to 1.000) | 0.480 (0.480 to 0.480) |
| 0.25 | 0.607 (0.607 to 0.607) | 0.740 (0.740 to 0.740) |
| 0.5 | 0.607 (0.607 to 0.607) | 0.740 (0.740 to 0.740) |
| 1 | 0.607 (0.607 to 0.607) | 0.740 (0.740 to 0.740) |
| composition alone | 0.204 (0.204 to 0.204) | 0.358 (0.358 to 0.358) |
| same-size divisions around random centers (median of 400) | 0.558 | 0.485 |
