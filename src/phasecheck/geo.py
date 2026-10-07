"""Distances and a local map projection from latitude and longitude.

Two products, used for different things:

- `great_circle_km` gives the pairwise distances (km) that the boundary
  excess (question 2 and 3) and the copying model (questions 6 to 9) use
  when no distance matrix is supplied.
- `planar_km` gives flat east/north coordinates (km) from which the
  alternative divisions of the map, their spread, the clustering on location
  (question 1) and the nearest-neighbor distances of question 4 are built.
  These always come from the coordinates, even when a distance matrix along
  rivers or trails is supplied.

The flat map is refused for a study area over `MAX_EXTENT_KM` across or
centered poleward of 70 degrees, rather than silently distorted.
"""
from __future__ import annotations

import numpy as np

EARTH_RADIUS_KM = 6371.0088         # mean Earth radius (IUGG), for the haversine
KM_PER_DEG = 111.32                 # km per degree of latitude, as the paper's analyses use
MAX_EXTENT_KM = 1000.0              # beyond this the flat-map approximation is refused


def great_circle_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Pairwise great-circle distance in kilometers (haversine).

    Parameters
    ----------
    lat, lon : array_like, shape (n,)
        Decimal degrees; latitude in -90..90, longitude in -180..180.

    Returns
    -------
    numpy.ndarray, shape (n, n)
        Symmetric distances in km, zero on the diagonal.
    """
    la, lo = np.deg2rad(np.asarray(lat, float)), np.deg2rad(np.asarray(lon, float))
    dla = la[:, None] - la[None, :]
    dlo = lo[:, None] - lo[None, :]
    # Clipping guards arcsin against rounding just above 1 for antipodal pairs.
    h = np.sin(dla / 2) ** 2 + np.cos(la)[:, None] * np.cos(la)[None, :] * np.sin(dlo / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(h, 0.0, 1.0)))


def planar_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Coordinates as centered kilometers, equirectangular about the mean latitude.

    Returns an (n, 2) array of east and north offsets in km from the
    centroid. The approximation is good to well under 1 percent over a few
    hundred kilometers at middle latitudes.

    Parameters
    ----------
    lat, lon : array_like, shape (n,)
        Decimal degrees.

    Returns
    -------
    numpy.ndarray, shape (n, 2)
        Column 0 is km east, column 1 km north, both centered on the mean.

    Raises
    ------
    ValueError
        If the mean latitude is poleward of 70 degrees (a common sign of
        swapped columns), or the points span more than `MAX_EXTENT_KM`
        (the diagonal of their bounding box on the flat map).
    """
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    # Longitudes are taken relative to their circular mean, so a region that
    # straddles 180 degrees is not read as spanning the globe.
    r = np.deg2rad(lon)
    center = np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean()))
    lon = (lon - center + 180.0) % 360.0 - 180.0
    mlat = float(lat.mean())
    if abs(mlat) > 70:
        raise ValueError(f"mean latitude {mlat:.1f} is too far poleward for the flat-map approximation; "
                         "check that the latitude and longitude columns are not swapped")
    p = np.column_stack([lon * np.cos(np.deg2rad(mlat)) * KM_PER_DEG, lat * KM_PER_DEG])
    p = p - p.mean(0)
    # Extent is the diagonal of the bounding box, a slight overstatement of the
    # largest pairwise distance, so the refusal errs toward caution.
    extent = float(np.sqrt(((p.max(0) - p.min(0)) ** 2).sum()))
    if extent > MAX_EXTENT_KM:
        raise ValueError(f"the assemblages span about {extent:.0f} km; the flat-map approximation used for "
                         f"dividing the map is limited to {MAX_EXTENT_KM:.0f} km. Analyze a smaller region.")
    return p
