#!/usr/bin/env python3
import tifffile
import numpy as np
import argparse


def postprocess_mask(input_tiff, mask_path, force_transpose=False):
    """
    Post-process mesmer output mask to fix transposed dimensions.

    Mesmer sometimes outputs masks with swapped X,Y dimensions.
    This script detects and corrects that issue.

    Use --force-transpose to always transpose regardless of shape check,
    which is useful when the mask spatial orientation is visibly flipped.
    """
    # Read input image to get expected dimensions
    input_img = tifffile.imread(input_tiff)
    if input_img.ndim == 3:
        expected_shape = input_img.shape[1:]  # (Y, X) from (C, Y, X)
    else:
        expected_shape = input_img.shape

    # Read mesmer output mask
    mask = tifffile.imread(mask_path)

    if force_transpose:
        print(f"Force-transposing mask from {mask.shape} to {mask.shape[::-1]}")
        mask = mask.T
        tifffile.imwrite(mask_path, mask, compression='deflate')
    # Check if dimensions are swapped (mask is X,Y instead of Y,X)
    elif mask.shape != expected_shape and mask.shape == expected_shape[::-1]:
        print(f"Transposing mask from {mask.shape} to {expected_shape}")
        mask = mask.T
        tifffile.imwrite(mask_path, mask, compression='deflate')
    else:
        print(f"Mask dimensions correct: {mask.shape}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Post-process mesmer segmentation masks')
    parser.add_argument('input_tiff', help='Original input TIFF file')
    parser.add_argument('mask_path', help='Mesmer output mask to post-process')
    parser.add_argument('--force-transpose', action='store_true',
                        help='Always transpose the mask, regardless of shape check')

    args = parser.parse_args()
    postprocess_mask(args.input_tiff, args.mask_path, force_transpose=args.force_transpose)
