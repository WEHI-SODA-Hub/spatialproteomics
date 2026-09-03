"""Unit tests for ``bin/aviti_wholecell_segment.py``'s pure-numpy mask
post-processing (small-object removal).

These do not import cellpose/torch (those imports live inside ``main()``),
so they run without a GPU or a Cellpose environment.
"""

import numpy as np

from aviti_wholecell_segment import remove_small_cells


def test_remove_small_cells_zeroes_labels_below_min_area():
    mask = np.array([[1, 1, 2], [1, 1, 2], [3, 3, 3]], dtype=np.uint16)
    # label 1: area 4, label 2: area 2, label 3: area 3
    out = remove_small_cells(mask, min_area=3)
    expected = np.array([[1, 1, 0], [1, 1, 0], [3, 3, 3]], dtype=np.uint16)
    np.testing.assert_array_equal(out, expected)


def test_remove_small_cells_noop_when_min_area_is_zero():
    mask = np.array([[1, 2], [0, 2]], dtype=np.uint16)
    out = remove_small_cells(mask, min_area=0)
    np.testing.assert_array_equal(out, mask)
