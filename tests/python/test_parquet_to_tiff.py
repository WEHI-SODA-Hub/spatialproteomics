"""Unit tests for label-dtype selection in ``bin/parquet_to_tiff.py``.

Label IDs are stored as pixel values, so a mask dtype too narrow for the largest
ID silently truncates the segmentation. This is worth pinning because there is
no error to notice: ``rasterio.features.rasterize`` does not raise on overflow,
it caps at the dtype maximum and every cell past that point vanishes. It only
bites on whole-slide images, which is exactly where nobody is watching.

Pure python -- no pipeline run required.
"""

import numpy as np
import rasterio.features
from shapely.geometry import box

import tifffile

from parquet_to_tiff import (
    LABEL_DTYPE,
    MASK_COMPRESSION,
    MASK_COMPRESSION_ARGS,
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
# the label dtype
# ----------------------------------------------------------------------------


def test_label_dtype_indexes_far_more_cells_than_any_real_slide():
    """uint32 is the wide fixed choice: 4.29e9 IDs, well past any real mask."""
    assert np.dtype(LABEL_DTYPE).name == "uint32"
    assert np.iinfo(LABEL_DTYPE).max == UINT32_MAX


# ----------------------------------------------------------------------------
# the truncation this guards against
# ----------------------------------------------------------------------------


def test_uint16_silently_truncates_beyond_its_maximum():
    """Documents the failure mode uint32 avoids: no exception, missing cells."""
    mask = _rasterize(UINT16_MAX + 2, np.uint16)

    assert mask.max() == UINT16_MAX
    assert len(np.unique(mask)) - 1 == UINT16_MAX  # two labels lost, silently


def test_label_dtype_preserves_every_label_past_uint16():
    n = UINT16_MAX + 2
    mask = _rasterize(n, LABEL_DTYPE)

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
