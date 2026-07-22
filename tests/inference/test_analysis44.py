import importlib
import numpy as np


def _mod():
    import sys
    from pathlib import Path
    ROOT = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "analyses"))
    return importlib.import_module("44_bayesian_fst_validation")


def test_simulate_bn_higher_f_gives_higher_plugin():
    from mls_emergence.signatures.variance import cultural_fst
    a44 = _mod()
    p_anc = np.array([0.4, 0.35, 0.25])
    sizes = np.array([200, 200, 200])
    rng = np.random.default_rng(0)
    lo = np.mean([cultural_fst(a44.simulate_bn(0.02, p_anc, sizes, rng)[0].astype(float))
                  for _ in range(15)])
    hi = np.mean([cultural_fst(a44.simulate_bn(0.30, p_anc, sizes, rng)[0].astype(float))
                  for _ in range(15)])
    assert hi > lo
