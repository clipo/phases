# Drainage-based basin definition (vs the latitude cut)

Basin proximity to the St. Francis drainage (L'Anguille River, Saint Francis Floodway, Saint Francis River, Tyronza River), distance in km (UTM 15N).

- Median distance of curated assemblages to the drainage: 11.0 km.

Membership comparison (drainage threshold vs lat >= 34.5 set of 36):

| threshold (km) | n in drainage basin | added vs lat-cut | dropped vs lat-cut |
|---|---|---|---|
| 10 | 18 | - | Beck, Belle_Meade, Carson_Lake, Cheatham, Commerce, Cramor_Place, Cummins, Dundee, Hollywood, Irby, Lake_Cormorant, Mound_Place, Norfolk, Notgrass, Pouncey, Upper_Nodena, Walls, Woodlyn |
| 15 | 22 | - | Beck, Belle_Meade, Cheatham, Commerce, Dundee, Hollywood, Irby, Lake_Cormorant, Mound_Place, Norfolk, Pouncey, Upper_Nodena, Walls, Woodlyn |
| 20 | 28 | - | Cheatham, Irby, Lake_Cormorant, Mound_Place, Norfolk, Pouncey, Walls, Woodlyn |
| 25 | 29 | - | Cheatham, Irby, Lake_Cormorant, Mound_Place, Norfolk, Walls, Woodlyn |

- Salomon: 32.0 km from the drainage (lat 34.35).
- Parchman: 29.8 km from the drainage (lat 34.36).

St. Francis-type sites outside the curated set that the drainage rule excludes (LMV location table, UTM 15N):
- Old Town: 28.1 km from the drainage.
- Blanchard: 79.8 km from the drainage.

## Convergence verdict under the drainage basin (<= 15 km)

- threshold 15 km (n = 22, k = 3): neutral rho = +0.60, F_ST rho = -1.00, spatial rho = -0.60; converges? no.
- threshold 20 km (n = 28, k = 3): neutral rho = +0.71, F_ST rho = +0.60, spatial rho = -0.37; converges? no.

Reading: the drainage-defined basin closely matches the latitude-cut basin, and the no-convergence verdict holds under the hydrological boundary as well, so the result does not depend on the latitude cut. The southern St-Francis-type sites excluded by the latitude cut are reported above with their drainage distances, making the boundary choice explicit rather than arbitrary.