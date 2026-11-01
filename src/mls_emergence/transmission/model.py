"""Couple the ceramic-copying simulator to the monument-mls emergence engine.

The bridge has three stages:

1.  ``phi_trajectory`` integrates the parent replicator dynamics to obtain phi(t),
    the fraction of monument-signaling (cooperation-relevant) groups over time, at
    a bistable operating point so trajectories are non-trivial.
2.  ``simulate_copying`` maps phi(t) to a ceramic-style assortment level
    ``a = coupling * phi(t)`` and generates a group-by-type count slice at each
    time via the shared transmission simulator (``simulate_slice``). The single
    parameter ``coupling`` in [0,1] is the load-bearing assumption: how strongly
    ceramic-style assortment tracks the cooperation-relevant assortment carried by
    monuments. coupling=1 means perfect tracking; coupling=0 means decoupled.
3.  ``emit_signatures`` / ``coupling_robustness`` push the resulting slices through
    the same four-signature pipeline used for the blind validation, so the
    convergence criterion is applied identically to the coupled model. The
    robustness sweep characterizes the coupling range over which the model yields a
    detectable convergent signature.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from mls_emergence.signatures import convergence
from mls_emergence.validation.harness import SIGNATURE_COLUMNS, signatures_over_axis
from mls_emergence.validation.mechanisms import simulate_slice

try:
    from signaling.emergence import replicator_dynamics
except ImportError:  # monument-mls optional; only the (paper-cut) coupling demo uses it
    replicator_dynamics = None

# Chosen bistable operating point. phi_star(0.5, 0.5) ~= 0.4299 is a finite
# interior saddle, so phi=0 and phi=1 are both attracting and phi_star +/- 0.1 are
# valid initial conditions on opposite sides. Scanned over sigma in [0.3,0.8] and
# lambda_W in [0.5,2.0]; this point keeps the saddle comfortably away from the
# absorbing boundaries (others, e.g. sigma=0.8,lambda_W=0.5, push phi_star to
# ~0.11 where phi_star-0.1 is nearly at the boundary).
SIGMA = 0.5
LAMBDA_W = 0.5

# Number of ordinal slices generated from a (subsampled) phi trajectory.
N_SLICES_DEFAULT = 8


def phi_trajectory(
    sigma: float, lambda_W: float, phi_0: float, t_max: float = 200.0
) -> np.ndarray:
    """Return phi over time from the parent replicator dynamics at (sigma, lambda_W).

    phi is the fraction of monument-signaling groups; phi=0 and phi=1 are
    absorbing and the interior saddle is phi_star(sigma, lambda_W).
    """
    if replicator_dynamics is None:
        raise ImportError(
            "phi_trajectory requires the monument-mls package (signaling); "
            "install with: pip install -e ../monument-mls"
        )
    out = replicator_dynamics(sigma, lambda_W, phi_0, t_max=t_max)
    return np.asarray(out["phi"], dtype=float)


def simulate_copying(
    phi_t: np.ndarray,
    coupling: float,
    n_groups: int,
    n_per_group: int,
    n_types: int,
    coords: np.ndarray,
    seed: int,
    n_slices: int = N_SLICES_DEFAULT,
) -> list[np.ndarray]:
    """Generate ceramic-copying slices driven by the monument-emergence phi.

    For each (subsampled) time step the latent ceramic assortment level is
    ``a = clip(coupling * phi_t[t], 0, 1)`` and a G x K count slice is drawn via
    ``simulate_slice`` with ``between_divergence = within_conformity = a`` on the
    bounded spatial rule. coupling scales how much phi propagates into ceramic
    assortment: coupling=0 holds a=0 (no emergence in the ceramic record regardless
    of monument dynamics); coupling=1 makes ceramic assortment track phi exactly.
    """
    phi_t = np.asarray(phi_t, dtype=float)
    if phi_t.size == 0:
        raise ValueError("phi_t is empty")
    # Subsample to a manageable, evenly spaced set of slices.
    idx = np.linspace(0, phi_t.size - 1, num=min(n_slices, phi_t.size)).round().astype(int)
    phi_sub = phi_t[idx]

    rng = np.random.default_rng(seed)
    coords = np.asarray(coords, dtype=float)
    slices: list[np.ndarray] = []
    for phi in phi_sub:
        a = float(np.clip(coupling * phi, 0.0, 1.0))
        slices.append(
            simulate_slice(
                n_groups=n_groups,
                n_per_group=n_per_group,
                n_types=n_types,
                between_divergence=a,
                within_conformity=a,
                spatial_rule="bounded",
                coords=coords,
                rng=rng,
            )
        )
    return slices


def emit_signatures(slices: list[np.ndarray], coords: np.ndarray) -> pd.DataFrame:
    """Compute the four cultural-transmission signatures over the slice sequence."""
    return signatures_over_axis(slices, coords)


def coupling_robustness(
    sigma: float,
    lambda_W: float,
    phi_0: float,
    couplings,
    n_groups: int,
    n_per_group: int,
    n_types: int,
    coords: np.ndarray,
    seed: int,
    t_max: float = 200.0,
    n_slices: int = N_SLICES_DEFAULT,
) -> pd.DataFrame:
    """Sweep coupling and report the convergence trend produced by each value.

    For every coupling: simulate the phi-driven ceramic record, emit the signature
    panel, and compute (a) the OLS ordinal slope of the combined convergence score
    and (b) the per-signature ordinal slopes. ``all_trend_up`` flags whether all
    four signatures rise (positive slope), the convergence criterion in raw form.

    The resulting table characterizes the coupling range over which the
    monument-phi-driven ceramic assortment yields a detectable convergent
    signature: the sensitivity of the empirical method to the shared-latent-
    assortment assumption.
    """
    phi_t = phi_trajectory(sigma, lambda_W, phi_0, t_max=t_max)
    rows = []
    for coupling in couplings:
        slices = simulate_copying(
            phi_t,
            coupling=coupling,
            n_groups=n_groups,
            n_per_group=n_per_group,
            n_types=n_types,
            coords=coords,
            seed=seed,
            n_slices=n_slices,
        )
        panel = emit_signatures(slices, coords)
        conv_trend = convergence.time_derivative(convergence.convergence_score(panel))
        row = {
            "coupling": float(coupling),
            "convergence_trend": float(conv_trend),
        }
        per_sig_up = []
        for col in SIGNATURE_COLUMNS:
            slope = convergence.time_derivative(panel[col])
            row[f"trend_{col}"] = float(slope)
            per_sig_up.append(slope > 0)
        row["all_trend_up"] = bool(all(per_sig_up))
        rows.append(row)
    return pd.DataFrame(rows)
