def test_monument_mls_emergence_importable():
    from signaling import emergence
    assert hasattr(emergence, "replicator_dynamics")
    assert hasattr(emergence, "phi_star")

def test_local_package_importable():
    import mls_emergence
    assert mls_emergence is not None
