"""Unit tests for the global-normalisation bounds in ``bin/compute_cellpose_norm.py``.

Cellpose normalises every patch to that patch's own 1-99 percentile, so an
identical cell lands on a different scale in a dense versus a sparse patch -- the
residual tiling grid once sopa's CLAHE/gaussian are off. This script computes one
[low, high] pair per channel from the whole image so every patch shares a scale.
The bit that must not regress is the numeric core: correct percentiles, the
memory-guard subsample, and the flat-channel guard (Cellpose divides by
high - low, so equal bounds would give inf/NaN).

Pure python -- no pipeline run, and no xarray/spatialdata import (those are lazy
in the script, so only its pure-numpy helper is exercised here).
"""

import numpy as np

from compute_cellpose_norm import _percentile_bounds


def test_bounds_are_the_requested_percentiles():
    """A ramp's 1/99 percentiles are its near-extremes; bounds match np.percentile."""
    values = np.arange(0, 10000, dtype=np.float32).reshape(100, 100)
    lo, hi = _percentile_bounds(values, 1.0, 99.0, max_pixels=10_000_000)
    exp_lo, exp_hi = np.percentile(values.ravel(), [1.0, 99.0])
    assert lo == float(exp_lo)
    assert hi == float(exp_hi)


def test_flat_channel_does_not_produce_a_zero_division_range():
    """A constant channel has lo == hi; the guard nudges hi so hi - lo > 0."""
    values = np.full((50, 50), 1234.0, dtype=np.float32)
    lo, hi = _percentile_bounds(values, 1.0, 99.0, max_pixels=10_000_000)
    assert lo == 1234.0
    assert hi > lo  # Cellpose divides by (hi - lo); must never be zero


def test_returns_plain_python_floats_for_literal_output():
    """The result is JSON/ast.literal_eval-safe: plain floats, not numpy scalars."""
    values = np.arange(400, dtype=np.float32).reshape(20, 20)
    bounds = _percentile_bounds(values, 1.0, 99.0, max_pixels=10_000_000)
    assert all(type(b) is float for b in bounds)


def test_subsample_kicks_in_above_the_cap_and_stays_close():
    """Above max_pixels the array is strided; 1/99 percentiles barely move."""
    rng = np.random.default_rng(0)
    values = rng.integers(0, 65535, size=(2000, 2000)).astype(np.float32)
    full = _percentile_bounds(values, 1.0, 99.0, max_pixels=10_000_000)  # no subsample
    capped = _percentile_bounds(values, 1.0, 99.0, max_pixels=100_000)  # forces stride
    # Same order of magnitude and within a few percent -- subsampling is only a
    # memory guard, not a change of statistic.
    assert abs(capped[0] - full[0]) <= max(5.0, 0.05 * abs(full[0]) + 1.0)
    assert abs(capped[1] - full[1]) <= 0.05 * full[1]
