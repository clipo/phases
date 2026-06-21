# Drainage-based basin definition (vs the latitude cut)

Basin proximity to the St. Francis drainage (L'Anguille River, Saint Francis Floodway, Saint Francis River, Tyronza River), distance in km (UTM 15N).

- Median distance of curated assemblages to the drainage: 18.5 km.

Membership comparison (drainage threshold vs lat >= 34.5 set of 53):

| threshold (km) | n in drainage basin | added vs lat-cut | dropped vs lat-cut |
|---|---|---|---|
| 10 | 19 | - | 40LA007, 40TP026, Beck, Belle_Meade, Bishop, Carson_Lake, Cheatham, Chuccalissa, Commerce, Cramor_Place, Cummins, Dundee, Fullen, Graves_Lake, Hatchie, Hollywood, Irby, Jeter, Jones_Bayou, Lake_Cormorant, Mound_Place, Norfolk, Notgrass, Porter, Pouncey, Rast, Richardsons_Landing, Upper_Nodena, Wall, Walls, West_Mounds, Wilder, Woodlyn, Young |
| 15 | 25 | - | 40LA007, 40TP026, Beck, Bishop, Cheatham, Chuccalissa, Commerce, Dundee, Fullen, Graves_Lake, Hatchie, Hollywood, Irby, Jeter, Jones_Bayou, Lake_Cormorant, Mound_Place, Norfolk, Porter, Pouncey, Rast, Richardsons_Landing, Upper_Nodena, Wall, Walls, Wilder, Woodlyn, Young |
| 20 | 29 | - | 40LA007, 40TP026, Beck, Bishop, Cheatham, Chuccalissa, Fullen, Graves_Lake, Hatchie, Irby, Jeter, Jones_Bayou, Lake_Cormorant, Mound_Place, Norfolk, Porter, Pouncey, Rast, Richardsons_Landing, Wall, Walls, Wilder, Woodlyn, Young |
| 25 | 32 | - | 40LA007, Beck, Bishop, Cheatham, Chuccalissa, Fullen, Graves_Lake, Hatchie, Irby, Jeter, Jones_Bayou, Lake_Cormorant, Mound_Place, Norfolk, Porter, Rast, Wall, Walls, Wilder, Woodlyn, Young |

- Salomon: 32.0 km from the drainage (lat 34.35).
- Parchman: 29.8 km from the drainage (lat 34.36).

## Convergence verdict under the drainage basin (<= 15 km)

- threshold 15 km (n = 25, k = 3): neutral rho = +1.00, F_ST rho = -0.60, spatial rho = +0.03; converges? no.
- threshold 20 km (n = 29, k = 3): neutral rho = +0.94, F_ST rho = -0.60, spatial rho = -0.03; converges? no.

Reading: the drainage-defined basin closely matches the latitude-cut basin, and the no-convergence verdict holds under the hydrological boundary as well, so the result does not depend on the latitude cut. The southern St-Francis-type sites excluded by the latitude cut are reported above with their drainage distances, making the boundary choice explicit rather than arbitrary.