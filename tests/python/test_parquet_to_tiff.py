"""Unit tests for label-dtype selection in ``bin/parquet_to_tiff.py``.

Label IDs are stored as pixel values, so a mask dtype too narrow for the largest
ID silently truncates the segmentation. This is worth pinning because there is
no error to notice: ``rasterio.features.rasterize`` does not raise on overflow,
it caps at the dtype maximum and every cell past that point vanishes. It only
bites on whole-slide images, which is exactly where nobody is watching.

Pure python -- no pipeline run required.
"""

import numpy as np
import pytest
import rasterio.features
from shapely.geometry import box

import tifffile

from parquet_to_tiff import (
    LABEL_DTYPES,
    MASK_COMPRESSION,
    MASK_COMPRESSION_ARGS,
    select_label_dtype,
)

UINT16_MAX = np.iinfo(np.uint16).max  # 65535
UINT32_MAX = np.iinfo(np.uint32).max


def _label_shapes(n):
    """``n`` disjoint 1px cells with 1-based IDs, laid out in a 100-wide grid."""
    return ((box(i % 100, i // 100, i % 100 + 0.9, i // 100 + 0.9), i) for i in range(1, n + 1))


def _rasterize(n, dtype):
    height = (n // 100) + 2
    return rasterio.features.rasterize(
        _label_shapes(n), out_shape=(height, 110), fill=0, dtype=dtype, all_touched=True
    )


# ----------------------------------------------------------------------------
# dtype selection
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "n_labels, expected",
    [
        (0, "uint16"),
        (1, "uint16"),
        (UINT16_MAX - 1, "uint16"),
        (UINT16_MAX, "uint16"),  # boundary: still fits
        (UINT16_MAX + 1, "uint32"),  # boundary: must widen
        (200_000, "uint32"),
        (UINT32_MAX, "uint32"),
    ],
)
def test_select_label_dtype_picks_narrowest_that_fits(n_labels, expected):
    assert select_label_dtype(n_labels).name == expected


def test_select_label_dtype_rejects_counts_beyond_widest_dtype():
    with pytest.raises(ValueError, match="exceeds the"):
        select_label_dtype(UINT32_MAX + 1)


def test_label_dtypes_are_ordered_narrowest_first():
    sizes = [np.dtype(dt).itemsize for dt in LABEL_DTYPES]
    assert sizes == sorted(sizes)


# ----------------------------------------------------------------------------
# the truncation this guards against
# ----------------------------------------------------------------------------


def test_uint16_silently_truncates_beyond_its_maximum():
    """Documents the failure mode: no exception, just missing cells."""
    mask = _rasterize(UINT16_MAX + 2, np.uint16)

    assert mask.max() == UINT16_MAX
    assert len(np.unique(mask)) - 1 == UINT16_MAX  # two labels lost, silently


def test_selected_dtype_preserves_every_label():
    n = UINT16_MAX + 2
    mask = _rasterize(n, select_label_dtype(n))

    assert mask.dtype == np.uint32
    assert mask.max() == n
    assert len(np.unique(mask)) - 1 == n


# ----------------------------------------------------------------------------
# compression
# ----------------------------------------------------------------------------


def test_mask_compression_is_lossless_and_smaller(tmp_path):
    """These masks are written to and re-read from shared scratch, so raw byte
    count is the cost. Compression must not change a single label."""
    mask = _rasterize(5000, np.uint32)
    path = tmp_path / "mask.tiff"

    tifffile.imwrite(
        path, mask,
        compression=MASK_COMPRESSION, compressionargs=MASK_COMPRESSION_ARGS,
    )

    assert path.stat().st_size < mask.nbytes
    assert np.array_equal(tifffile.imread(path), mask)
