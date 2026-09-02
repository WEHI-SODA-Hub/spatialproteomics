#!/usr/bin/env python
'''
Module      : aviti_wholecell_segment
Description : Whole-cell/membrane segmentation of a single AVITI cell-paint
              tile using a Cellpose v4 Segment-Anything ("cpsam") model,
              fed the Nucleus + Cell-Membrane (+ Actin, 3-channel mode)
              channels stacked into one multi-channel image.

              This is an independent implementation of the general approach
              described in Elembio's own (BSD-licensed) analysis notebook --
              stack the cell-paint channels and run a Cellpose-SAM model per
              tile -- not a port of that notebook's code. Normalisation is
              left to Cellpose's own built-in percentile normalisation
              (channel-independent 1st/99th percentile rescale), matching
              this pipeline's existing Cellpose preprocessing contract
              (see docs/usage.md) rather than reimplementing a bespoke
              normalisation step.
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

# Label masks are long runs of a repeated ID, so deflate shrinks them ~90x.
MASK_COMPRESSION = "zlib"
MASK_COMPRESSION_ARGS = {"level": 1}


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def read_channel(path: Path) -> np.ndarray:
    img = tifffile.imread(path)
    if img.ndim != 2:
        raise ValueError(f"Expected a single-page 2D channel image, got shape {img.shape} from {path}")
    return img


def build_stack(nucleus_tif: Path, membrane_tif: Path, actin_tif: Optional[Path]) -> np.ndarray:
    '''
    Returns a channel-last (Y, X, C) stack for Cellpose's channel_axis=-1.
    '''
    nucleus = read_channel(nucleus_tif)
    membrane = read_channel(membrane_tif)
    if membrane.shape != nucleus.shape:
        raise ValueError(
            f"Channel shape mismatch: nucleus {nucleus.shape} vs membrane {membrane.shape}"
        )

    channels = [nucleus, membrane]
    if actin_tif is not None:
        actin = read_channel(actin_tif)
        if actin.shape != nucleus.shape:
            raise ValueError(
                f"Channel shape mismatch: nucleus {nucleus.shape} vs actin {actin.shape}"
            )
        channels.append(actin)

    return np.stack(channels, axis=-1)


def remove_small_cells(mask: np.ndarray, min_area: int) -> np.ndarray:
    if min_area <= 0:
        return mask
    labels, counts = np.unique(mask, return_counts=True)
    small = labels[(labels != 0) & (counts < min_area)]
    if small.size == 0:
        return mask
    out = mask.copy()
    out[np.isin(out, small)] = 0
    return out


@app.command()
def main(
    nucleus_tif: Annotated[Path, typer.Option(exists=True, help="Nucleus channel tile TIFF.")],
    membrane_tif: Annotated[Path, typer.Option(exists=True, help="Cell-Membrane channel tile TIFF.")],
    output: Annotated[Path, typer.Option(help="Path to write the whole-cell label mask TIFF.")],
    actin_tif: Annotated[Optional[Path], typer.Option(
        exists=True, help="Actin channel tile TIFF (3-channel mode only)."
    )] = None,
    pretrained_model: Annotated[str, typer.Option(
        help="Built-in Cellpose model name (e.g. cpsam_v2) or path to a custom model."
    )] = "cpsam_v2",
    diameter: Annotated[float, typer.Option(help="Expected cell diameter in pixels.")] = 30.0,
    flow_threshold: Annotated[float, typer.Option()] = 0.4,
    cellprob_threshold: Annotated[float, typer.Option()] = 0.0,
    min_area: Annotated[int, typer.Option(help="Discard cells smaller than this many px^2. 0 disables.")] = 0,
    gpu: Annotated[bool, typer.Option(help="Run on GPU.")] = True,
):
    '''
    Segment whole cells in one AVITI tile with a Cellpose-SAM model.
    '''
    import torch
    from cellpose import models

    log(f"PyTorch version: {torch.__version__}")
    log(f"CUDA available: {torch.cuda.is_available()}")
    use_gpu = gpu and torch.cuda.is_available()
    if gpu and not use_gpu:
        log("WARNING: --gpu requested but no CUDA device is available; running on CPU.")

    stack = build_stack(nucleus_tif, membrane_tif, actin_tif)
    log(f"Built {stack.shape[-1]}-channel stack with shape {stack.shape} for whole-cell segmentation")

    model = models.CellposeModel(gpu=use_gpu, pretrained_model=pretrained_model)
    masks, _flows, _styles = model.eval(
        stack,
        channel_axis=-1,
        diameter=diameter,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )

    masks = remove_small_cells(masks.astype(np.uint32), min_area)

    tifffile.imwrite(
        output, masks,
        compression=MASK_COMPRESSION, compressionargs=MASK_COMPRESSION_ARGS,
    )
    n_cells = len(np.unique(masks)) - 1
    log(f"Whole-cell segmentation: {n_cells} cell(s) written to {output}")


if __name__ == "__main__":
    app()
