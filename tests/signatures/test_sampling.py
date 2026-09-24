import numpy as np
import pytest
from mls_emergence.signatures.sampling import rarefy


def test_full_sample_is_unchanged():
    a = np.array([[25, 25], [1, 0], [0, 0]])
    np.testing.assert_array_equal(rarefy(a, 50, np.random.default_rng(0)), a)


def test_subsample_is_a_subset_with_requested_total():
    a = np.array([[1, 40, 80], [30, 50, 0]])
    for seed in range(20):
        out = rarefy(a, 50, np.random.default_rng(seed))
        assert (out <= a).all()
        np.testing.assert_array_equal(out.sum(1), [50, 50])


def test_fractional_counts_are_rejected():
    with pytest.raises(ValueError, match="integer"):
        rarefy([[1.5, 2]], 1, np.random.default_rng(0))
