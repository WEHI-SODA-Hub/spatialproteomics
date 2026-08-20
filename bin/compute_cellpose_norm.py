#!/usr/bin/env python
'''
Module      : compute_cellpose_norm
Description : Compute per-channel global intensity bounds (low, high) for a
              SpatialData .zarr image, to be handed to Cellpose as its
              `normalize={"lowhigh": ...}` argument.

              Cellpose normalises every patch to that patch's own 1-99
              percentile. Identical cells therefore land on different scales in
              a dense versus a sparse patch, which is the source of the residual
              tiling grid once sopa's CLAHE/gaussian stages are disabled.
              Computing one (low, high) pair per channel from the WHOLE image
              and applying it to every patch removes that patch dependence.

              Bounds are emitted in the SAME channel order they are passed on
              the command line, which must match the order in which the caller
              passes `--channels` to `sopa segmentation cellpose` (membrane then
              nucleus for whole-cell, nucleus alone for nuclear-only). The
              output is a JSON/Python-literal list of [low, high] pairs, ready to
              splice into `--method-kwargs`.
Copyright   : (c) WEHI SODA Hub, 2025
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import json
from pathlib import Path
from typing import Annotated, List, Optional

import numpy as np
import typer

app = typer.Typer()


def _get_image(zarr_path: str):
    """Return the single image element of a SpatialData store as a DataArray or
    DataTree (multiscale). Uses sopa's helper when available, else falls back to
    the sole image in the store."""
    import spatialdata as sd

    sdata = sd.read_zarr(zarr_path)
    try:
        from sopa.utils import get_spatial_image

        return get_spatial_image(sdata)
    except Exception:
        keys = list(sdata.images.keys())
        if len(keys) != 1:
            raise ValueError(
                f"Expected exactly one image in {zarr_path}, found {keys}"
            )
        return sdata.images[keys[0]]


def _pick_level(image, max_pixels: int):
    """From a (possibly multiscale) image, return a single DataArray to compute
    percentiles on. For a multiscale pyramid, pick the finest level whose Y*X is
    at or below `max_pixels` (falling back to the coarsest if even that is
    larger). p1/p99 are robust to this downsampling, and it bounds memory on
    large slides."""
    # Multiscale images are DataTree-like: they expose scale groups.
    scales = None
    try:
        scales = [k for k in image.keys() if str(k).startswith("scale")]
    except Exception:
        scales = None

    if scales:
        # scale0 is full resolution; scales are ordered coarser as the index grows.
        def _yx(scale):
            arr = image[scale]["image"]
            return int(arr.sizes["y"]) * int(arr.sizes["x"])

        scales_sorted = sorted(scales, key=lambda s: _yx(s), reverse=True)
        chosen = scales_sorted[-1]  # coarsest as a safe default
        for s in scales_sorted:  # finest first
            if _yx(s) <= max_pixels:
                chosen = s
                break
        return image[chosen]["image"]

    # Single-scale DataArray.
    return image


def _percentile_bounds(values: np.ndarray, plow: float, phigh: float,
                       max_pixels: int) -> List[float]:
    """Compute [low, high] percentile bounds from a 2D intensity array.

    Pure numpy so it is unit-testable without xarray/spatialdata. Stride-
    subsamples a single-scale array far above the cap to bound memory, and
    guards a flat channel: Cellpose divides by (high - low), so equal bounds
    would give inf/NaN."""
    values = np.asarray(values)
    n = values.size
    if n > max_pixels:
        stride = int(np.ceil(np.sqrt(n / max_pixels)))
        values = values[::stride, ::stride]

    lo, hi = np.percentile(values.ravel(), [plow, phigh])
    lo, hi = float(lo), float(hi)
    if hi <= lo:
        hi = lo + 1.0
    return [lo, hi]


def _bounds_for_channel(arr, channel: str, plow: float, phigh: float,
                        max_pixels: int) -> List[float]:
    """Compute [low, high] percentile bounds for one named channel."""
    a = arr.sel(c=channel).data  # dask or numpy, shape (Y, X)
    return _percentile_bounds(np.asarray(a), plow, phigh, max_pixels)


@app.command()
def main(
    zarr: Annotated[str, typer.Argument(help="Path to the SpatialData .zarr store")],
    channel: Annotated[
        List[str],
        typer.Option(
            "--channel",
            help="Channel name to compute bounds for. Repeat in the SAME order "
            "the channels are passed to `sopa segmentation cellpose` "
            "(membrane then nucleus for whole-cell).",
        ),
    ],
    percentiles: Annotated[
        str,
        typer.Option(help="Low,high percentiles for the bounds."),
    ] = "1,99",
    max_pixels: Annotated[
        int,
        typer.Option(
            help="Compute on the finest multiscale level with at most this many "
            "pixels (single-scale images above it are stride-subsampled).",
        ),
    ] = 8_000_000,
    as_method_kwargs: Annotated[
        bool,
        typer.Option(
            "--as-method-kwargs",
            help="Emit the full Cellpose method-kwargs dict "
            '`{"normalize": {"lowhigh": ...}}` instead of the bare list of '
            "pairs, ready to pass to `sopa segmentation cellpose --method-kwargs`.",
        ),
    ] = False,
    output: Annotated[
        Optional[str],
        typer.Option("--output", "-o", help="Write JSON here instead of stdout."),
    ] = None,
):
    """Emit per-channel [low, high] bounds as a JSON list of pairs, ordered to
    match the given `--channel` order."""
    plow, phigh = (float(x) for x in percentiles.split(","))

    image = _get_image(zarr)
    arr = _pick_level(image, max_pixels)

    available = {str(c) for c in arr.coords["c"].values}
    missing = [c for c in channel if c not in available]
    if missing:
        raise ValueError(
            f"Channel(s) {missing} not found in image. Available: {sorted(available)}"
        )

    lowhigh = [
        _bounds_for_channel(arr, c, plow, phigh, max_pixels) for c in channel
    ]

    payload = {"normalize": {"lowhigh": lowhigh}} if as_method_kwargs else lowhigh
    text = json.dumps(payload)
    if output:
        Path(output).write_text(text)
    else:
        typer.echo(text)


if __name__ == "__main__":
    app()
