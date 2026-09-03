"""Unit tests for ``bin/aviti_nuclear_segment.py``'s pure-numpy mask
post-processing (small-object removal, binarization to Elembio's
Nuclear.tif convention).

These do not import cellpose/torch (those imports live inside ``main()``),
so they run without a GPU or a Cellpose environment.
"""

import numpy as np

from aviti_nuclear_segment import binarize, remove_small_cells


def test_remove_small_cells_zeroes_labels_below_min_area():
    mask = np.array([[1, 1, 2], [1, 1, 2], [3, 3, 3]], dtype=np.int32)
    # label 1: area 4, label 2: area 2, label 3: area 3
    out = remove_small_cells(mask, min_area=3)
    expected = np.array([[1, 1, 0], [1, 1, 0], [3, 3, 3]], dtype=np.int32)
    np.testing.assert_array_equal(out, expected)


def test_remove_small_cells_noop_when_min_area_is_zero():
    mask = np.array([[1, 2], [0, 2]], dtype=np.int32)
    out = remove_small_cells(mask, min_area=0)
    np.testing.assert_array_equal(out, mask)


def test_binarize_collapses_instance_labels_to_zero_one():
    # Per Elembio's notebook, Nuclear.tif is a uint8 0/1 presence mask, not
    # instance-labeled -- distinct nucleus IDs (1, 2, 3, ...) must all
    # collapse to 1.
    mask = np.array([[0, 1, 2], [3, 0, 7]], dtype=np.int32)
    out = binarize(mask)
    expected = np.array([[0, 1, 1], [1, 0, 1]], dtype=np.uint8)
    np.testing.assert_array_equal(out, expected)
    assert out.dtype == np.uint8


def test_binarize_handles_a_high_instance_count_without_wraparound():
    # A real AVITI tile can contain thousands of nuclei -- confirming that
    # binarizing (rather than narrowing instance IDs to uint8 directly)
    # avoids the label-collision-via-modulo-wraparound risk that motivated
    # this fix: label values here deliberately exceed uint8's range.
    mask = np.array([[0, 300, 1000, 65000]], dtype=np.int32)
    out = binarize(mask)
    np.testing.assert_array_equal(out, np.array([[0, 1, 1, 1]], dtype=np.uint8))


def test_remove_small_cells_then_binarize_preserves_presence_after_filtering():
    mask = np.array([[1, 1, 2, 3]], dtype=np.int32)  # label areas: 1->2, 2->1, 3->1
    filtered = remove_small_cells(mask, min_area=2)  # drops labels 2 and 3
    out = binarize(filtered)
    np.testing.assert_array_equal(out, np.array([[1, 1, 0, 0]], dtype=np.uint8))
