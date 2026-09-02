"""Unit tests for per-well stitching in ``bin/aviti_stitch.py``.

These exercise the pure-numpy stitching logic (offset computation, occupancy-
based overlap resolution, label relabeling) with tiny synthetic tiles -- no
real AVITI data, no GPU, no cellpose.
"""

import numpy as np
import pytest
import tifffile

from aviti_stitch import (
    compute_offsets,
    get_channel_names,
    mm_to_px,
    read_tile_rows,
    stitch_image,
    stitch_masks,
)


# ----------------------------------------------------------------------------
# coordinate conversion
# ----------------------------------------------------------------------------


def test_mm_to_px_uses_the_given_pixel_size():
    # 1 mm = 1000 um; at 0.25 um/px that's 4000 px.
    assert mm_to_px(1.0, pixel_size_microns=0.25) == 4000


def test_mm_to_px_rounds_to_nearest_pixel():
    assert mm_to_px(0.00013, pixel_size_microns=0.25) == round(0.00013 * 1000 / 0.25)


def test_compute_offsets_are_zero_based_from_the_minimum_tile():
    rows = [
        {"x_mm": "0.0", "y_mm": "0.0"},
        {"x_mm": "0.001", "y_mm": "0.0"},  # 4 px to the right at 0.25 um/px
    ]
    offsets = compute_offsets(rows, pixel_size_microns=0.25)
    assert offsets == [(0, 0), (4, 0)]


# ----------------------------------------------------------------------------
# manifest reading
# ----------------------------------------------------------------------------


def test_read_tile_rows_orders_deterministically_by_tile_name(tmp_path):
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(
        "tile,x_mm,y_mm,cell_mask,nuclear_mask,image_tif\n"
        "tileB,1,0,b_cell.tif,b_nuc.tif,b_img.tif\n"
        "tileA,0,0,a_cell.tif,a_nuc.tif,a_img.tif\n"
    )
    rows = read_tile_rows(manifest)
    assert [r["tile"] for r in rows] == ["tileA", "tileB"]


def test_read_tile_rows_fails_on_an_empty_manifest(tmp_path):
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("tile,x_mm,y_mm,cell_mask,nuclear_mask,image_tif\n")
    with pytest.raises(ValueError, match="No rows"):
        read_tile_rows(manifest)


# ----------------------------------------------------------------------------
# mask stitching: placement and overlap resolution
# ----------------------------------------------------------------------------


def _write_mask(path, mask):
    tifffile.imwrite(path, mask.astype(np.uint32))


def test_non_overlapping_tiles_are_placed_at_their_offsets(tmp_path):
    tile_a = np.array([[1, 1], [0, 0]], dtype=np.uint32)
    tile_b = np.array([[2, 2], [0, 0]], dtype=np.uint32)
    _write_mask(tmp_path / "a.tif", tile_a)
    _write_mask(tmp_path / "b.tif", tile_b)

    rows = [{"cell_mask": str(tmp_path / "a.tif")}, {"cell_mask": str(tmp_path / "b.tif")}]
    offsets = [(0, 0), (2, 0)]

    canvas = stitch_masks(rows, offsets, "cell_mask", canvas_h=2, canvas_w=4)

    assert np.array_equal(canvas[:, 0:2], tile_a)
    # Tile B's label 2 is offset by tile A's running max label (1) -> becomes 3.
    assert np.array_equal(canvas[:, 2:4], np.array([[3, 3], [0, 0]], dtype=np.uint32))


def test_first_tile_wins_on_overlap(tmp_path):
    tile_a = np.array([[1, 1]], dtype=np.uint32)
    tile_b = np.array([[2, 2]], dtype=np.uint32)
    _write_mask(tmp_path / "a.tif", tile_a)
    _write_mask(tmp_path / "b.tif", tile_b)

    rows = [{"cell_mask": str(tmp_path / "a.tif")}, {"cell_mask": str(tmp_path / "b.tif")}]
    # Fully overlapping placement.
    offsets = [(0, 0), (0, 0)]

    canvas = stitch_masks(rows, offsets, "cell_mask", canvas_h=1, canvas_w=2)

    # Tile A (placed first) keeps every pixel; tile B contributes nothing.
    assert np.array_equal(canvas, tile_a)


def test_background_pixels_do_not_block_a_later_tile(tmp_path):
    # Tile A has a background (0) hole; tile B should be free to fill it.
    tile_a = np.array([[1, 0]], dtype=np.uint32)
    tile_b = np.array([[0, 2]], dtype=np.uint32)
    _write_mask(tmp_path / "a.tif", tile_a)
    _write_mask(tmp_path / "b.tif", tile_b)

    rows = [{"cell_mask": str(tmp_path / "a.tif")}, {"cell_mask": str(tmp_path / "b.tif")}]
    offsets = [(0, 0), (0, 0)]

    canvas = stitch_masks(rows, offsets, "cell_mask", canvas_h=1, canvas_w=2)

    # Tile A's real label 1 stays; tile A's background at column 1 is filled
    # by tile B's label (offset by tile A's running max, 1 -> 3), not
    # silently zeroed by the "background wins" bug this guards against.
    assert canvas[0, 0] == 1
    assert canvas[0, 1] == 3


def test_relabeling_keeps_ids_unique_across_many_tiles(tmp_path):
    tiles = [np.array([[i + 1]], dtype=np.uint32) for i in range(3)]
    paths = []
    for i, tile in enumerate(tiles):
        p = tmp_path / f"tile{i}.tif"
        _write_mask(p, tile)
        paths.append(p)

    rows = [{"cell_mask": str(p)} for p in paths]
    offsets = [(i, 0) for i in range(3)]

    canvas = stitch_masks(rows, offsets, "cell_mask", canvas_h=1, canvas_w=3)

    # Every placed label must be unique: no two tiles' "label 1" collide.
    labels = canvas[canvas > 0]
    assert len(labels) == len(set(labels.tolist()))


# ----------------------------------------------------------------------------
# image stitching
# ----------------------------------------------------------------------------


def test_stitch_image_places_multichannel_tiles_and_reads_channel_names(tmp_path):
    img_a = np.zeros((2, 1, 2), dtype=np.uint16)
    img_a[0] = [[10, 20]]
    img_a[1] = [[30, 40]]
    img_b = np.zeros((2, 1, 2), dtype=np.uint16)
    img_b[0] = [[50, 60]]
    img_b[1] = [[70, 80]]

    path_a = tmp_path / "a.ome.tif"
    path_b = tmp_path / "b.ome.tif"
    tifffile.imwrite(
        path_a, img_a, ome=True,
        metadata={"axes": "CYX", "Channel": {"Name": ["Nucleus", "Cell-Membrane"]}},
    )
    tifffile.imwrite(path_b, img_b)

    rows = [{"image_tif": str(path_a)}, {"image_tif": str(path_b)}]
    offsets = [(0, 0), (2, 0)]

    canvas, channel_names = stitch_image(rows, offsets, canvas_h=1, canvas_w=4)

    assert channel_names == ["Nucleus", "Cell-Membrane"]
    assert canvas.shape == (2, 1, 4)
    assert np.array_equal(canvas[:, :, 0:2], img_a)
    assert np.array_equal(canvas[:, :, 2:4], img_b)


def test_get_channel_names_falls_back_to_generic_names_without_ome_metadata(tmp_path):
    path = tmp_path / "plain.tif"
    tifffile.imwrite(path, np.zeros((3, 2, 2), dtype=np.uint16), photometric="minisblack")
    assert get_channel_names(path, 3) == ["Channel_0", "Channel_1", "Channel_2"]
