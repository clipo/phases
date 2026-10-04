# The assemblages and Phillips's (1970) phase areas

Written by `analyses/93_phillips1970_phase_map.py`. The areas are those of
`data/raw/phillips1970_phase_outlines.csv` (Lipo 2001: Figure 2.6, after Phillips 1970).
The analysis assigns each assemblage the area it lies in, or the nearest when it lies in
none (`36_canonical_phase_map.assign_primary_phases`). Containment is point-in-polygon in
UTM 15N; the areas carry the source's drafting error (see the CSV header), so an assemblage
within 1.25 km of a line is not securely on either side of it. The last two columns
record the Mainfort (1996: Figure 1) labels the pipeline used until 2026-10-01.

- Drawn in Figure 1: 28 analyzed assemblages and 8 other assemblages of the
  38-assemblage matrix; 2 fall outside the frame.
- Of the 28 analyzed assemblages, 26 lie inside a phase area and 2 take the
  nearest: Cramor_Place (Kent, 1.0 km outside the line); Holden_Lake (Parkin, 1.0 km outside the line).
- By phase: Nodena 0, Parkin 12, Walls 5, Kent 11.
- 5 analyzed assemblages carry a different phase than under the 1996 labels: Beck (Walls to Kent); Belle_Meade (Walls to Kent); Castile_Landing (Kent to Parkin); Commerce (Walls to Kent); Hollywood (Walls to Kent).
- Area sizes (km2): Nodena 649, Parkin 1199, Walls 640, Kent 1722. No two areas overlap.

| assemblage | set | Phillips phase | assigned by | km outside the line | Mainfort 1996 label (superseded) | changed |
|---|---|---|---|---|---|---|
| Barton_Ranch | analyzed | Parkin | inside the area | 0 | Parkin |  |
| Beck | analyzed | Kent | inside the area | 0 | Walls | changed |
| Belle_Meade | analyzed | Kent | inside the area | 0 | Walls | changed |
| Big_Eddy | analyzed | Parkin | inside the area | 0 | Parkin |  |
| Carson_Lake | not analyzed | Nodena | inside the area | 0 | Nodena |  |
| Castile_Landing | analyzed | Parkin | inside the area | 0 | Kent | changed |
| Cheatham | not analyzed | Walls | inside the area | 0 | Walls |  |
| Clay_Hill | analyzed | Kent | inside the area | 0 | Kent |  |
| Commerce | analyzed | Kent | inside the area | 0 | Walls | changed |
| Connor | not analyzed | Kent | inside the area | 0 | Kent |  |
| Cramor_Place | analyzed | Kent | nearest area | 1.0 | Kent |  |
| Cummins | analyzed | Parkin | inside the area | 0 | Parkin |  |
| Davis | analyzed | Kent | inside the area | 0 | Kent |  |
| Dundee | not analyzed | Kent | nearest area | 3.9 | Parchman | changed |
| Fortune | analyzed | Parkin | inside the area | 0 | Parkin |  |
| Grant | analyzed | Kent | inside the area | 0 | Kent |  |
| Holden_Lake | analyzed | Parkin | nearest area | 1.0 | Parkin |  |
| Hollywood | analyzed | Kent | inside the area | 0 | Walls | changed |
| Irby | analyzed | Walls | inside the area | 0 | Walls |  |
| Kent_Place | analyzed | Kent | inside the area | 0 | Kent |  |
| Lake_Cormorant | analyzed | Walls | inside the area | 0 | Walls |  |
| Mound_Place | analyzed | Walls | inside the area | 0 | Walls |  |
| Neeleys_Ferry | analyzed | Parkin | inside the area | 0 | Parkin |  |
| Nickel | analyzed | Kent | inside the area | 0 | Kent |  |
| Norfolk | not analyzed | Walls | inside the area | 0 | Walls |  |
| Notgrass | not analyzed | Nodena | inside the area | 0 | Nodena |  |
| Parkin | analyzed | Parkin | inside the area | 0 | Parkin |  |
| Pouncey | not analyzed | Walls | inside the area | 0 | Walls |  |
| Rose_Mound | analyzed | Parkin | inside the area | 0 | Parkin |  |
| Starkley | analyzed | Kent | inside the area | 0 | Kent |  |
| Turnbow | analyzed | Parkin | inside the area | 0 | Parkin |  |
| Upper_Nodena | not analyzed | Nodena | inside the area | 0 | Nodena |  |
| Vernon_Paul | analyzed | Parkin | inside the area | 0 | Parkin |  |
| Walls | analyzed | Walls | inside the area | 0 | Walls |  |
| Williamson | analyzed | Parkin | inside the area | 0 | Parkin |  |
| Woodlyn | analyzed | Walls | inside the area | 0 | Walls |  |
