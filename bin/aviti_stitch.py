#!/usr/bin/env python
'''
Module      : aviti_stitch
Description : Stitches per-tile whole-cell masks, nuclear masks, and merged
              intensity images for one AVITI well into single per-well
              TIFFs, using the tile stage coordinates (XMillimeters,
              YMillimeters) discovered from RunParameters.json.

              Tiles are placed on a shared canvas at their nominal stage
              position converted to pixels via --pixel-size-microns. Where
              tiles overlap, whichever tile is placed first (deterministic
              tile-name order) wins the disputed pixels outright, for both
              masks and the intensity image -- there is no seam blending and
              no attempt to reconcile a cell split across two tiles' claimed
              regions. This matches AVITI's own tiling, which is a
              near-contiguous grid rather than the heavily overlapping
              patches this pipeline uses elsewhere for whole-slide COMET
              images.

              Label IDs are made unique across tiles by adding a running
              offset per mask type (cell, nuclear) before compositing, so a
              stitched well mask never silently merges two different tiles'
              cell 1 into one label.
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


def mm_to_px(value_mm: float, pixel_size_microns: float) -> int:
    return round(float(value_mm) * 1000.0 / pixel_size_microns)


def compute_offsets(rows: List[dict], pixel_size_microns: float):
    xs = [mm_to_px(r["x_mm"], pixel_size_microns) for r in rows]
    ys = [mm_to_px(r["y_mm"], pixel_size_microns) for r in rows]
    min_x, min_y = min(xs), min(ys)
    return [(x - min_x, y - min_y) for x, y in zip(xs, ys)]


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


def stitch_masks(rows, offsets, path_key, canvas_h, canvas_w):
    canvas = np.zeros((canvas_h, canvas_w), dtype=LABEL_DTYPE)
    occupied = np.zeros((canvas_h, canvas_w), dtype=bool)
    running_max = 0

    for row, (x0, y0) in zip(rows, offsets):
        tile_mask = tifffile.imread(row[path_key]).astype(LABEL_DTYPE)
        h, w = tile_mask.shape
        region_occupied = occupied[y0:y0 + h, x0:x0 + w]
        new_pixels = ~region_occupied

        offset_mask = tile_mask.copy()
        offset_mask[offset_mask > 0] += running_max

        region_canvas = canvas[y0:y0 + h, x0:x0 + w]
        region_canvas[new_pixels] = offset_mask[new_pixels]
        canvas[y0:y0 + h, x0:x0 + w] = region_canvas
        occupied[y0:y0 + h, x0:x0 + w] = region_occupied | (tile_mask > 0)

        running_max = int(max(running_max, offset_mask.max()))

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
        help="AVITI pixel size in microns/pixel, used to convert stage "
             "XMillimeters/YMillimeters into pixel offsets."
    )] = 0.25,
):
    '''
    Stitch one well's per-tile masks and intensity image into per-well TIFFs.
    '''
    rows = read_tile_rows(manifest)
    offsets = compute_offsets(rows, pixel_size_microns)

    tile_shape = tifffile.imread(rows[0]["cell_mask"]).shape
    tile_h, tile_w = tile_shape
    canvas_w = max(x for x, _ in offsets) + tile_w
    canvas_h = max(y for _, y in offsets) + tile_h
    log(f"Stitching {len(rows)} tile(s) into a {canvas_w}x{canvas_h}px canvas")

    cell_canvas = stitch_masks(rows, offsets, "cell_mask", canvas_h, canvas_w)
    tifffile.imwrite(
        output_cell, cell_canvas,
        compression=MASK_COMPRESSION, compressionargs=MASK_COMPRESSION_ARGS,
    )

    nuclear_canvas = stitch_masks(rows, offsets, "nuclear_mask", canvas_h, canvas_w)
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
        f"nuclear mask ({int(nuclear_canvas.max())} max label), and "
        f"{len(channel_names)}-channel image written."
    )


if __name__ == "__main__":
    app()
