#!/usr/bin/env python
'''
Module      : aviti_nuclear_segment
Description : Nuclear segmentation of a single AVITI cell-paint tile using a
              Cellpose 3.x custom model trained specifically for AVITI nuclei
              (e.g. 20250212_cellpose_nuc_8diam), run in a separate
              environment from the Cellpose v4 SAM whole-cell path because the
              two major Cellpose versions are not installable side by side.

              Independent implementation -- not derived from Elembio's own
              (BSD-licensed) analysis notebook.
Copyright   : (c) WEHI SODA Hub, 2026
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import sys
from pathlib import Path
from typing import Annotated

import numpy as np
import tifffile
import typer

app = typer.Typer(add_completion=False)

# Label masks are long runs of a repeated ID, so deflate shrinks them ~90x.
MASK_COMPRESSION = "zlib"
MASK_COMPRESSION_ARGS = {"level": 1}


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


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
    nucleus_tif: Annotated[Path, typer.Argument(exists=True, help="Nucleus channel tile TIFF.")],
    output: Annotated[Path, typer.Option(help="Path to write the nuclear label mask TIFF.")],
    model_path: Annotated[Path, typer.Option(
        exists=True,
        help="Path to the custom Cellpose 3.x nuclear model "
             "(e.g. 20250212_cellpose_nuc_8diam)."
    )],
    diameter: Annotated[float, typer.Option(
        help="Expected nuclear diameter in pixels. Matches the model's own "
             "training diameter unless you have a reason to override it."
    )] = 8.0,
    flow_threshold: Annotated[float, typer.Option()] = 0.4,
    cellprob_threshold: Annotated[float, typer.Option()] = 0.0,
    min_area: Annotated[int, typer.Option(help="Discard nuclei smaller than this many px^2. 0 disables.")] = 0,
    gpu: Annotated[bool, typer.Option(help="Run on GPU.")] = True,
):
    '''
    Segment nuclei in one AVITI tile with a custom Cellpose 3.x model.
    '''
    import torch
    from cellpose import models

    log(f"PyTorch version: {torch.__version__}")
    log(f"CUDA available: {torch.cuda.is_available()}")
    use_gpu = gpu and torch.cuda.is_available()
    if gpu and not use_gpu:
        log("WARNING: --gpu requested but no CUDA device is available; running on CPU.")

    img = tifffile.imread(nucleus_tif)
    if img.ndim != 2:
        raise ValueError(f"Expected a single-page 2D nucleus image, got shape {img.shape}")

    model = models.CellposeModel(gpu=use_gpu, pretrained_model=str(model_path))
    # channels=[0, 0]: grayscale input with no separate nuclear channel index --
    # the whole image *is* the nuclear channel, which is Cellpose's convention
    # for single-channel input on the pre-v4 API this model was trained under.
    masks, _flows, _styles = model.eval(
        img,
        channels=[0, 0],
        diameter=diameter,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )

    masks = remove_small_cells(masks.astype(np.uint32), min_area)

    tifffile.imwrite(
        output, masks,
        compression=MASK_COMPRESSION, compressionargs=MASK_COMPRESSION_ARGS,
    )
    n_nuclei = len(np.unique(masks)) - 1
    log(f"Nuclear segmentation: {n_nuclei} nucleus/nuclei written to {output}")


if __name__ == "__main__":
    app()
