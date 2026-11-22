import pytest


def test_monument_mls_emergence_importable():
    # The parent model is used only by the earlier version's coupling code and is installed separately.
    emergence = pytest.importorskip("signaling.emergence", reason="monument-mls (module signaling) is not installed")
    assert hasattr(emergence, "replicator_dynamics")
    assert hasattr(emergence, "phi_star")

def test_local_package_importable():
    import mls_emergence
    assert mls_emergence is not None
