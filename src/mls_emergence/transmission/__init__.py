"""Transmission-layer coupling between the ceramic-copying simulator and the
parent monument-mls emergence engine.

This package makes explicit the paper's "shared latent assortment" bridge: the
monument-emergence dynamics (signaling.emergence.replicator_dynamics) produce a
fraction phi of cooperation-signaling groups over time, and a coupling parameter
controls how strongly that phi drives the latent assortment level of the ceramic
style-copying process. The criterion's ability to detect emergence therefore
depends on this coupling, which is the load-bearing assumption tested here.
"""
from mls_emergence.transmission.model import (
    LAMBDA_W,
    SIGMA,
    coupling_robustness,
    emit_signatures,
    phi_trajectory,
    simulate_copying,
)

__all__ = [
    "SIGMA",
    "LAMBDA_W",
    "phi_trajectory",
    "simulate_copying",
    "emit_signatures",
    "coupling_robustness",
]
