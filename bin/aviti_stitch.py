#!/usr/bin/env python
'''
Module      : aviti_stitch
Description : Stitches per-tile whole-cell masks, nuclear masks, and merged
              intensity images for one AVITI well into single per-well
              TIFFs.

              Tiles are placed on a *grid*, not at a literal pixel conversion
              of their stage coordinates (XMillimeters/YMillimeters). AVITI's
              stage frame is not aligned with its own image row/column frame:
              treating stage millimetres as image pixels directly produces
              both wildly oversized, non-uniform gaps between tiles (stage
              spacing does not match the exported field-of-view pixel size)
              and a well that is rotated 90 degrees relative to how AVITI's
              own Cytocanvas viewer lays it out.

              Each tile's stage X/Y is instead used only to *rank* it among
              the other tiles in its well, and those ranks are mapped onto
              a Cytocanvas-orientated pixel grid with a small fixed visual
              gap (--tile-gap-microns) between adjacent tiles:

                  output_row    = reverse_rank(x_mm)       descending
                  output_column = rank(y_mm)               ascending

              i.e. tiles with a larger stage X land in earlier (higher, closer
              to the top of the image) output rows, and tiles with a smaller
              stage Y land in earlier (further left) output columns. This
              90-degree rotation was determined empirically by comparing this
              pipeline's naive stage-coordinate placement against AVITI's own
              Cytocanvas layout for the same well, and is not otherwise
              documented by Elembio.

              Label IDs are made unique across tiles by adding a running
              offset per tile for the whole-cell mask, so a stitched well
              mask never silently merges two different tiles' cell 1 into
              one label. The nuclear mask is Elembio's binary Nuclear.tif
              convention (0/1 presence, not instance-labeled), so it is
              composited as-is with no offsetting.

              Where tiles overlap, whichever tile is placed first
              (deterministic tile-name order) wins the disputed pixels
              outright, for both masks and the intensity image -- there is
              no seam blending and no attempt to reconcile a cell split
              across two tiles' claimed regions.
Copyright   : (c) WEHI SODA Hub, 2026
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import csv
import sys
from pathlib import Path
from typing import Annotated, List

import numpy as np
import tifffile
import typer

app = typer.Typer(add_completion=False)

MASK_COMPRESSION = "zlib"
MASK_COMPRESSION_ARGS = {"level": 1}
# As we need to stitch several images together, use larger uint to avoid overflow
LABEL_DTYPE = np.uint32


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def read_tile_rows(manifest: Path) -> List[dict]:
    with open(manifest, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"No rows found in stitch manifest {manifest}")
    # Deterministic placement order, independent of CSV row order.
    return sorted(rows, key=lambda r: r["tile"])


def microns_to_px(value_microns: float, pixel_size_microns: float) -> int:
    return round(float(value_microns) / pixel_size_microns)


def compute_grid_offsets(rows: List[dict], tile_h: int, tile_w: int, gap_px: int):
    '''
    Rank each tile's stage X/Y among its well's other tiles and map those
    ranks onto a Cytocanvas-orientated pixel grid, spaced by the tile size
    plus a fixed visual gap (rather than by stage distance).

    See the module docstring for why the axes are rotated/reversed rather
    than mapped straight through: output_row = reverse_rank(x_mm) descending,
    output_column = rank(y_mm) ascending.

    Returns (offsets, n_rows, n_cols), where offsets are (x0, y0) pixel
    positions in the same order as ``rows``, and n_rows/n_cols are the
    resulting grid dimensions (for logging/sanity-checking a well that is
    missing tiles).
    '''
    xs = [float(r["x_mm"]) for r in rows]
    ys = [float(r["y_mm"]) for r in rows]

    x_values = sorted(set(xs))
    y_values = sorted(set(ys))
    x_rank = {x: i for i, x in enumerate(x_values)}
    y_rank = {y: i for i, y in enumerate(y_values)}
    n_rows, n_cols = len(x_values), len(y_values)

    offsets = []
    for x, y in zip(xs, ys):
        out_row = (n_rows - 1) - x_rank[x]
        out_col = y_rank[y]
        offsets.append((out_col * (tile_w + gap_px), out_row * (tile_h + gap_px)))

    return offsets, n_rows, n_cols


def get_channel_names(tiff_path: Path, n_channels: int) -> List[str]:
    try:
        with tifffile.TiffFile(tiff_path) as tif:
            if tif.is_ome:
                from xml.etree import ElementTree as ET
                root = ET.fromstring(tif.ome_metadata)
                names = [
                    ch.get("Name") for ch in root.iter()
                    if ch.tag.split("}", 1)[-1] == "Channel"
                ]
                if len(names) == n_channels and all(names):
                    return names
    except Exception:
        pass
    return [f"Channel_{i}" for i in range(n_channels)]


def stitch_masks(rows, offsets, path_key, canvas_h, canvas_w, binary: bool = False):
    '''
    Composite one mask type across all tiles in a well.

    ``binary=True`` is for Elembio's Nuclear.tif convention (a 0/1 presence
    mask, not instance-labeled -- see aviti_nuclear_segment.py): pixels are
    copied through as-is with no running-label offset, since there are no
    per-tile instance IDs to keep unique. Cell masks (``binary=False``) are
    instance-labeled, so a running offset is added per tile to keep every
    tile's "cell 1" etc. distinct in the stitched well.
    '''
    dtype = np.uint8 if binary else LABEL_DTYPE
    canvas = np.zeros((canvas_h, canvas_w), dtype=dtype)
    occupied = np.zeros((canvas_h, canvas_w), dtype=bool)
    running_max = 0

    for row, (x0, y0) in zip(rows, offsets):
        tile_mask = tifffile.imread(row[path_key]).astype(dtype)
        h, w = tile_mask.shape
        region_occupied = occupied[y0:y0 + h, x0:x0 + w]
        new_pixels = ~region_occupied

        if binary:
            offset_mask = tile_mask
        else:
            offset_mask = tile_mask.copy()
            offset_mask[offset_mask > 0] += running_max
            running_max = int(max(running_max, offset_mask.max()))

        region_canvas = canvas[y0:y0 + h, x0:x0 + w]
        region_canvas[new_pixels] = offset_mask[new_pixels]
        canvas[y0:y0 + h, x0:x0 + w] = region_canvas
        occupied[y0:y0 + h, x0:x0 + w] = region_occupied | (tile_mask > 0)

    return canvas


def stitch_image(rows, offsets, canvas_h, canvas_w):
    first = tifffile.imread(rows[0]["image_tif"])
    if first.ndim == 2:
        first = first[np.newaxis, ...]
    n_channels = first.shape[0]
    dtype = first.dtype
    channel_names = get_channel_names(Path(rows[0]["image_tif"]), n_channels)

    canvas = np.zeros((n_channels, canvas_h, canvas_w), dtype=dtype)
    occupied = np.zeros((canvas_h, canvas_w), dtype=bool)

    for row, (x0, y0) in zip(rows, offsets):
        img = tifffile.imread(row["image_tif"])
        if img.ndim == 2:
            img = img[np.newaxis, ...]
        _, h, w = img.shape
        region_occupied = occupied[y0:y0 + h, x0:x0 + w]
        new_pixels = ~region_occupied

        region_canvas = canvas[:, y0:y0 + h, x0:x0 + w]
        for c in range(n_channels):
            channel_region = region_canvas[c]
            channel_region[new_pixels] = img[c][new_pixels]
            region_canvas[c] = channel_region
        canvas[:, y0:y0 + h, x0:x0 + w] = region_canvas
        occupied[y0:y0 + h, x0:x0 + w] = region_occupied | new_pixels

    return canvas, channel_names


@app.command()
def main(
    manifest: Annotated[Path, typer.Argument(
        exists=True,
        help="Per-well stitch manifest CSV with columns tile,x_mm,y_mm,"
             "cell_mask,nuclear_mask,image_tif."
    )],
    output_cell: Annotated[Path, typer.Option(help="Path to write the stitched whole-cell mask.")],
    output_nuclear: Annotated[Path, typer.Option(help="Path to write the stitched nuclear mask.")],
    output_image: Annotated[Path, typer.Option(help="Path to write the stitched intensity image.")],
    pixel_size_microns: Annotated[float, typer.Option(
        help="AVITI pixel size in microns/pixel, used only to convert "
             "--tile-gap-microns into a pixel gap. Tile placement is grid-"
             "based (see module docstring), not a literal stage-coordinate "
             "pixel conversion."
    )] = 0.25,
    tile_gap_microns: Annotated[float, typer.Option(
        help="Visual gap in microns to insert between adjacent tiles in the "
             "stitched output. This is purely for viewer clarity, not a "
             "measurement of the true physical stage gap between tiles."
    )] = 32.0,
):
    '''
    Stitch one well's per-tile masks and intensity image into per-well TIFFs.
    '''
    rows = read_tile_rows(manifest)

    tile_shape = tifffile.imread(rows[0]["cell_mask"]).shape
    tile_h, tile_w = tile_shape

    gap_px = microns_to_px(tile_gap_microns, pixel_size_microns)
    offsets, n_rows, n_cols = compute_grid_offsets(rows, tile_h, tile_w, gap_px)
    if n_rows * n_cols != len(rows):
        log(
            f"WARNING: {len(rows)} tile(s) span a {n_rows}x{n_cols} grid "
            f"({n_rows * n_cols} cells) -- this well may be missing tile(s)."
        )

    canvas_w = max(x for x, _ in offsets) + tile_w
    canvas_h = max(y for _, y in offsets) + tile_h
    log(
        f"Stitching {len(rows)} tile(s) into a {n_rows}x{n_cols} grid "
        f"({canvas_w}x{canvas_h}px canvas, {gap_px}px/{tile_gap_microns}um gap)"
    )

    cell_canvas = stitch_masks(rows, offsets, "cell_mask", canvas_h, canvas_w)
    tifffile.imwrite(
        output_cell, cell_canvas,
        compression=MASK_COMPRESSION, compressionargs=MASK_COMPRESSION_ARGS,
    )

    nuclear_canvas = stitch_masks(rows, offsets, "nuclear_mask", canvas_h, canvas_w, binary=True)
    tifffile.imwrite(
        output_nuclear, nuclear_canvas,
        compression=MASK_COMPRESSION, compressionargs=MASK_COMPRESSION_ARGS,
    )

    image_canvas, channel_names = stitch_image(rows, offsets, canvas_h, canvas_w)
    tifffile.imwrite(
        output_image, image_canvas,
        ome=True,
        metadata={"axes": "CYX", "Channel": {"Name": channel_names}},
    )

    log(
        f"Stitched cell mask ({int(cell_canvas.max())} max label), "
        f"nuclear mask ({int(nuclear_canvas.sum())} nucleus px), and "
        f"{len(channel_names)}-channel image written."
    )


if __name__ == "__main__":
    app()
