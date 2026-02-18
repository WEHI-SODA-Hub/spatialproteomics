#!/usr/bin/env python3
"""
Smooth label segmentation masks to reduce polygon complexity.

Supports two methods:
  - morphological: per-label binary close + open with a disk structuring element
  - gaussian: per-label Gaussian blur then re-threshold at 0.5

This prevents StackOverflowError in QuPath's GeoJSON export when cell ROIs have
overly complex/jagged polygon boundaries.
"""

import argparse
import sys

import numpy as np
import tifffile
from scipy import ndimage
from skimage.morphology import disk, binary_closing, binary_opening
from skimage.filters import gaussian


def smooth_label_morphological(label_mask, kernel_size=2):
    """
    Smooth each label in a label mask using morphological close then open.

    Parameters
    ----------
    label_mask : np.ndarray
        Integer label image where 0 is background.
    kernel_size : int
        Radius of the disk structuring element (default: 2).

    Returns
    -------
    np.ndarray
        Smoothed label mask with the same dtype.
    """
    selem = disk(kernel_size)
    labels = np.unique(label_mask)
    labels = labels[labels != 0]

    smoothed = np.zeros_like(label_mask)

    for label_id in labels:
        binary = label_mask == label_id
        # Close first (fills small holes, connects nearby regions)
        binary = binary_closing(binary, selem)
        # Then open (removes small protrusions, smooths boundary)
        binary = binary_opening(binary, selem)
        smoothed[binary & (smoothed == 0)] = label_id

    return smoothed


def smooth_label_gaussian(label_mask, sigma=1.0):
    """
    Smooth each label in a label mask using Gaussian blur + re-threshold.

    Parameters
    ----------
    label_mask : np.ndarray
        Integer label image where 0 is background.
    sigma : float
        Standard deviation of Gaussian kernel (default: 1.0).

    Returns
    -------
    np.ndarray
        Smoothed label mask with the same dtype.
    """
    labels = np.unique(label_mask)
    labels = labels[labels != 0]

    smoothed = np.zeros_like(label_mask)
    # Track confidence for overlap resolution
    confidence = np.zeros(label_mask.shape, dtype=np.float32)

    for label_id in labels:
        binary = (label_mask == label_id).astype(np.float32)
        blurred = gaussian(binary, sigma=sigma, preserve_range=True)
        mask = blurred > 0.5

        # Resolve overlaps: keep label with highest blur confidence
        overlap = mask & (smoothed != 0)
        if np.any(overlap):
            better = blurred > confidence
            smoothed[mask & ((smoothed == 0) | better)] = label_id
            confidence[mask & ((confidence < blurred))] = blurred[mask & ((confidence < blurred))]
        else:
            smoothed[mask] = label_id
            confidence[mask] = blurred[mask]

    return smoothed


def main():
    parser = argparse.ArgumentParser(
        description="Smooth label segmentation masks to reduce polygon complexity."
    )
    parser.add_argument(
        "input_mask",
        help="Path to input label mask TIFF file."
    )
    parser.add_argument(
        "output_mask",
        help="Path to output smoothed label mask TIFF file."
    )
    parser.add_argument(
        "--method",
        choices=["morphological", "gaussian"],
        default="morphological",
        help="Smoothing method: 'morphological' (close+open with disk) or "
             "'gaussian' (blur+threshold). Default: morphological."
    )
    parser.add_argument(
        "--kernel-size",
        type=float,
        default=2,
        help="For morphological: disk radius (int, default: 2). "
             "For gaussian: sigma value (float, default: 2)."
    )
    args = parser.parse_args()

    print(f"Reading mask: {args.input_mask}")
    mask = tifffile.imread(args.input_mask)
    original_dtype = mask.dtype
    print(f"Mask shape: {mask.shape}, dtype: {original_dtype}")

    n_labels = len(np.unique(mask)) - 1  # exclude background
    print(f"Number of labels: {n_labels}")

    if n_labels == 0:
        print("No labels found, copying input to output unchanged.")
        tifffile.imwrite(args.output_mask, mask)
        sys.exit(0)

    print(f"Smoothing with method={args.method}, kernel_size={args.kernel_size}")

    if args.method == "morphological":
        smoothed = smooth_label_morphological(mask, kernel_size=int(args.kernel_size))
    elif args.method == "gaussian":
        smoothed = smooth_label_gaussian(mask, sigma=args.kernel_size)
    else:
        print(f"Unknown method: {args.method}", file=sys.stderr)
        sys.exit(1)

    # Preserve original dtype
    smoothed = smoothed.astype(original_dtype)

    n_labels_after = len(np.unique(smoothed)) - 1
    print(f"Labels after smoothing: {n_labels_after} (lost {n_labels - n_labels_after})")

    print(f"Writing smoothed mask: {args.output_mask}")
    tifffile.imwrite(args.output_mask, smoothed)
    print("Done.")


if __name__ == "__main__":
    main()
