#!/usr/bin/env python
'''
Module      : aviti_merge_tile_channels
Description : Stacks the raw (non-normalised) AVITI cell-paint tile channels
              (Nucleus, Cell-Membrane, optionally Actin) into a single
              multi-channel OME-TIFF with named channels.

              This is the per-tile intensity image later stitched into a
              per-well image and fed to CELLMEASUREMENT/KRONOS2EMBEDDINGS --
              it deliberately carries raw pixel values, not the
              percentile-normalised stack used only as Cellpose input, so
              that downstream per-cell intensity measurements are computed
              from the instrument's own values.
Copyright   : (c) WEHI SODA Hub, 2026
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import sys
from pathlib import Path
from typing import Annotated, Optional

import numpy as np
import tifffile
import typer

app = typer.Typer(add_completion=False)


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def read_channel(path: Path) -> np.ndarray:
    img = tifffile.imread(path)
    if img.ndim != 2:
        raise ValueError(f"Expected a single-page 2D channel image, got shape {img.shape} from {path}")
    return img


def build_stack(nucleus_tif: Path, membrane_tif: Path, actin_tif: Optional[Path]):
    '''
    Returns (stack, channel_names) with stack shape (C, Y, X).
    '''
    nucleus = read_channel(nucleus_tif)
    membrane = read_channel(membrane_tif)

    if membrane.shape != nucleus.shape:
        raise ValueError(
            f"Channel shape mismatch: nucleus {nucleus.shape} vs membrane {membrane.shape}"
        )

    channels = [nucleus, membrane]
    names = ["Nucleus", "Cell-Membrane"]

    if actin_tif is not None:
        actin = read_channel(actin_tif)
        if actin.shape != nucleus.shape:
            raise ValueError(
                f"Channel shape mismatch: nucleus {nucleus.shape} vs actin {actin.shape}"
            )
        channels.append(actin)
        names.append("Actin")

    stack = np.stack(channels, axis=0)
    return stack, names


@app.command()
def main(
    nucleus_tif: Annotated[Path, typer.Option(exists=True, help="Raw Nucleus channel tile TIFF.")],
    membrane_tif: Annotated[Path, typer.Option(exists=True, help="Raw Cell-Membrane channel tile TIFF.")],
    output: Annotated[Path, typer.Option(help="Path to write the merged multi-channel OME-TIFF.")],
    actin_tif: Annotated[Optional[Path], typer.Option(
        exists=True, help="Raw Actin channel tile TIFF (3-channel mode only)."
    )] = None,
):
    '''
    Merge raw AVITI tile channels into one named multi-channel OME-TIFF.
    '''
    stack, names = build_stack(nucleus_tif, membrane_tif, actin_tif)

    tifffile.imwrite(
        output,
        stack,
        ome=True,
        metadata={"axes": "CYX", "Channel": {"Name": names}},
    )
    log(f"Wrote {len(names)}-channel image {names} with shape {stack.shape} to {output}")


if __name__ == "__main__":
    app()
