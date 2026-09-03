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

# Matches Elembio's own onboard segmentation output tag-for-tag (confirmed by
# comparing tifffile tags against a real onboard-segmented Cell.tif):
# ImageJ-format ImageDescription/hyperstack metadata with AdobeDeflate
# compression, not tifffile's default "shaped" JSON metadata.
MASK_COMPRESSION = "zlib"


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


def binarize(mask: np.ndarray) -> np.ndarray:
    '''
    Collapse an instance-labeled mask to Elembio's own Nuclear.tif
    convention: a uint8 mask where 0 means no nucleus and 1 means a nucleus
    is present, per their notebook's documented Nuclear.tif spec (unlike
    Cell.tif, which is a uint16 instance-labeled mask). Nuclear.tif is never
    instance-labeled, so there is no label-count-vs-dtype overflow concern
    here the way there is for the whole-cell mask.
    '''
    return (mask > 0).astype(np.uint8)


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

    masks = remove_small_cells(masks, min_area)
    n_nuclei = len(np.unique(masks)) - 1

    # Per Elembio's own notebook, Nuclear.tif is a uint8 *binary* mask (0 =
    # no nucleus, 1 = nucleus present), unlike Cell.tif which is uint16
    # instance-labeled -- so collapse instances to presence here rather than
    # keeping (and potentially overflowing) per-nucleus label IDs.
    masks = binarize(masks)

    # imagej=True matches Elembio's own onboard segmentation TIFFs (ImageJ-
    # format ImageDescription/hyperstack metadata); see the module-level
    # MASK_COMPRESSION comment for why this matters more than compression
    # itself.
    tifffile.imwrite(output, masks, imagej=True, compression=MASK_COMPRESSION)
    log(f"Nuclear segmentation: {n_nuclei} nucleus/nuclei written to {output}")


if __name__ == "__main__":
    app()
