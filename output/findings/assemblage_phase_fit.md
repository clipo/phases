# Which phase does each assemblage's pottery fit best?

Produced by `analyses/96_assemblage_phase_fit.py`. 28 assemblages, phases Kent 11, Parkin 12, Walls 5 (Phillips's areas). Ceramic distance is the chi-square distance from the assemblage's decorated-class profile to each phase's pooled profile, the assemblage left out of its own phase's pool. P(own) is the share of 4,000 draws of Dirichlet(counts + 1/2), for the assemblage and for every pool, in which its own phase is the closest. km is the straight-line distance to the nearest other member of each phase.

- By direct distance, 5 of 28 assemblages are closer to another phase's profile than to their own: Belle_Meade (Kent, closer to Parkin), Big_Eddy (Parkin, closer to Kent), Holden_Lake (Parkin, closer to Kent), Mound_Place (Walls, closer to Parkin), Turnbow (Parkin, closer to Kent).
- With sampling carried through, 22 fit their own phase best with probability above 0.75, 4 fit another phase best (probability of their own below 0.25), and 2 are about even between two phases (0.25 to 0.75): Castile_Landing, Mound_Place.
- 2 assemblages lie nearer a member of another phase than any other member of their own: Castile_Landing (Parkin, nearest Kent), Clay_Hill (Kent, nearest Parkin).

| assemblage | sherds | phase | ceramic distance to Kent / Parkin / Walls | closest | P(own) | km to nearest member of Kent / Parkin / Walls |
|---|---|---|---|---|---|---|
| Big_Eddy | 229 | Parkin | 0.929 / 1.177 / 1.330 | Kent | 0.00 | 22 / 15 / 42 |
| Holden_Lake | 360 | Parkin (nearest area) | 1.242 / 1.372 / 1.772 | Kent | 0.00 | 27 / 13 / 25 |
| Turnbow | 128 | Parkin | 0.940 / 1.024 / 1.490 | Kent | 0.00 | 48 / 3 / 48 |
| Belle_Meade | 1709 | Kent | 0.397 / 0.302 / 0.807 | Parkin | 0.01 | 1 / 27 / 16 |
| Mound_Place | 168 | Walls | 0.496 / 0.496 / 0.538 | Parkin | 0.40 | 22 / 25 / 11 |
| Castile_Landing | 614 | Parkin | 0.571 / 0.555 / 0.794 | Parkin | 0.60 | 12 / 16 / 42 |
| Cramor_Place | 424 | Kent (nearest area) | 0.185 / 0.259 / 0.834 | Kent | 0.87 | 7 / 18 / 24 |
| Nickel | 1300 | Kent | 0.275 / 0.344 / 0.837 | Kent | 0.95 | 7 / 12 / 30 |
| Rose_Mound | 983 | Parkin | 0.271 / 0.139 / 0.885 | Parkin | 0.97 | 28 / 7 / 38 |
| Beck | 1233 | Kent | 0.284 / 0.364 / 0.858 | Kent | 0.98 | 1 / 27 / 17 |
| Davis | 79 | Kent | 0.879 / 1.010 / 1.147 | Kent | 0.99 | 6 / 18 / 31 |
| Woodlyn | 107 | Walls | 0.854 / 0.829 / 0.305 | Walls | 0.99 | 16 / 37 / 3 |
| Fortune | 200 | Parkin | 0.529 / 0.290 / 0.895 | Parkin | 1.00 | 46 / 3 / 47 |
| Grant | 93 | Kent | 1.259 / 1.457 / 1.543 | Kent | 1.00 | 3 / 25 / 42 |
| Hollywood | 226 | Kent | 0.726 / 0.936 / 1.022 | Kent | 1.00 | 5 / 36 / 22 |
| Walls | 360 | Walls | 0.855 / 0.832 / 0.451 | Walls | 1.00 | 18 / 35 / 3 |
| Cummins | 261 | Parkin | 0.663 / 0.457 / 0.874 | Parkin | 1.00 | 50 / 5 / 47 |
| Lake_Cormorant | 192 | Walls | 0.994 / 1.124 / 0.635 | Walls | 1.00 | 18 / 41 / 4 |
| Commerce | 128 | Kent | 0.532 / 0.720 / 1.059 | Kent | 1.00 | 5 / 32 / 20 |
| Barton_Ranch | 638 | Parkin | 0.477 / 0.190 / 0.921 | Parkin | 1.00 | 42 / 12 / 36 |
| Clay_Hill | 181 | Kent | 0.607 / 0.775 / 1.014 | Kent | 1.00 | 17 / 12 / 47 |
| Irby | 91 | Walls | 1.570 / 1.694 / 1.101 | Walls | 1.00 | 22 / 44 / 4 |
| Kent_Place | 241 | Kent | 0.814 / 0.992 / 1.121 | Kent | 1.00 | 6 / 23 / 33 |
| Neeleys_Ferry | 1266 | Parkin | 0.410 / 0.223 / 0.917 | Parkin | 1.00 | 38 / 2 / 40 |
| Parkin | 1303 | Parkin | 0.454 / 0.248 / 0.780 | Parkin | 1.00 | 31 / 4 / 36 |
| Starkley | 171 | Kent | 0.836 / 1.035 / 1.192 | Kent | 1.00 | 3 / 25 / 39 |
| Vernon_Paul | 758 | Parkin | 0.512 / 0.309 / 0.810 | Parkin | 1.00 | 37 / 2 / 39 |
| Williamson | 658 | Parkin | 0.503 / 0.226 / 0.949 | Parkin | 1.00 | 34 / 4 / 40 |

## If Big Eddy and Castile Landing are counted as Kent

Recovery is the adjusted Rand index over all 28 assemblages (chi-square composition); the phase-line comparison is analysis 86's, with its seeds (2,000 posterior draws, 300 alternative divisions of each kind).

| | Phillips's areas | the two counted as Kent |
|---|---|---|
| phase sizes (Kent / Parkin / Walls) | 11 / 12 / 5 | 13 / 10 / 5 |
| site map alone, k-means | 0.509 | 0.666 |
| site map alone, Ward | 0.756 | 1.000 |
| map plus composition (weight 0.5), k-means | 0.756 | 1.000 |
| composition alone, k-means | 0.326 | 0.327 |
| F_ST, posterior median (95% CrI) | 0.0121 (0.0101 to 0.0144) | 0.0142 (0.0124 to 0.0162) |
| median F_ST, random-center / compact divisions | 0.0112 / 0.0103 | 0.0123 / 0.0121 |
| P(lines separate better than random-center / compact) | 0.57 / 0.72 | 0.61 / 0.57 |
| boundary excess, posterior median | +7.3 | +12.4 |
| P(larger boundary excess than random-center / compact) | 0.58 / 0.65 | 0.79 / 0.71 |
