#!/usr/bin/env python3
"""
Smooth label segmentation masks to reduce polygon complexity.

Uses morphological close + open per label, parallelised via shared memory
for scalability on large masks (e.g. 40k x 40k with 100k+ labels).

This prevents StackOverflowError in QuPath's GeoJSON export when cell ROIs
have overly complex/jagged polygon boundaries.
"""

import argparse
import sys
import numpy as np
import tifffile
from multiprocessing import shared_memory
from concurrent.futures import ProcessPoolExecutor
from scipy.ndimage import find_objects
from skimage.morphology import disk, binary_closing, binary_opening


# ── Worker (must be module-level for pickling) ────────────────────────────────

def _process_label(args):
    """
    Process a single label: extract crop → morphological close+open → return pixels.
    Operates in crop-local space then offsets to global coordinates.
    """
    (label_id, slc_starts, slc_stops, shm_name,
     mask_shape, mask_dtype, kernel_size) = args

    shm = shared_memory.SharedMemory(name=shm_name)
    label_mask = np.ndarray(mask_shape, dtype=mask_dtype, buffer=shm.buf)

    try:
        slc = tuple(slice(s, e) for s, e in zip(slc_starts, slc_stops))
        crop = label_mask[slc] == label_id

        if not np.any(crop):
            return label_id, None, None

        selem = disk(kernel_size)

        # Close first (fills small holes, connects nearby regions)
        binary = binary_closing(crop, selem)
        # Then open (removes small protrusions, smooths boundary)
        binary = binary_opening(binary, selem)

        if not np.any(binary):
            return label_id, None, None

        # Get local pixel coordinates and offset to global
        local_rr, local_cc = np.where(binary)
        minr, minc = slc_starts
        rr = local_rr + minr
        cc = local_cc + minc

        # Clip to mask bounds — guards against kernel expanding beyond crop edge
        valid = (
            (rr >= 0) & (rr < mask_shape[0]) &
            (cc >= 0) & (cc < mask_shape[1])
        )

        return label_id, rr[valid], cc[valid]

    except Exception as e:
        print(f"  [warn] label {label_id} failed: {e}", flush=True)
        # Fallback: return original pixels unchanged
        orig_rr, orig_cc = np.where(label_mask == label_id)
        return label_id, orig_rr, orig_cc

    finally:
        shm.close()  # Never unlink in workers — only the creator does that


# ── Main smoothing function ───────────────────────────────────────────────────

def smooth_label_morphological_parallel(
    label_mask, kernel_size=2, n_workers=8, chunksize=200
):
    """
    Parallelised morphological label smoothing using shared memory.

    Parameters
    ----------
    label_mask : np.ndarray
        Integer label image where 0 is background.
    kernel_size : int
        Radius of the disk structuring element (default: 2).
    n_workers : int
        Number of parallel worker processes.
    chunksize : int
        Number of labels per executor dispatch chunk.
    """
    print(f"Setting up shared memory for mask {label_mask.shape} {label_mask.dtype}...")
    shm = shared_memory.SharedMemory(create=True, size=label_mask.nbytes)
    shared_arr = np.ndarray(label_mask.shape, dtype=label_mask.dtype, buffer=shm.buf)
    np.copyto(shared_arr, label_mask)

    print("Finding label bounding boxes...")
    slices = find_objects(label_mask)
    n_labels = sum(1 for s in slices if s is not None)
    print(f"Processing {n_labels} labels with {n_workers} workers, chunksize={chunksize}...")

    # Serialise slices as (starts, stops) tuples — slice objects aren't picklable
    work = []
    for label_id, slc in enumerate(slices, start=1):
        if slc is None:
            continue
        slc_starts = tuple(s.start for s in slc)
        slc_stops  = tuple(s.stop  for s in slc)
        work.append((
            label_id, slc_starts, slc_stops,
            shm.name, label_mask.shape, label_mask.dtype,
            kernel_size
        ))

    smoothed = np.zeros_like(label_mask)
    done = 0
    warned = 0

    try:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            for label_id, rr, cc in executor.map(
                _process_label, work, chunksize=chunksize
            ):
                if rr is not None and len(rr) > 0:
                    # First-write wins for any overlap between adjacent labels
                    free = smoothed[rr, cc] == 0
                    smoothed[rr[free], cc[free]] = label_id
                else:
                    warned += 1
                done += 1
                if done % 5000 == 0:
                    print(f"  {done}/{n_labels} labels done...", flush=True)
    finally:
        shm.close()
        shm.unlink()  # Only unlink in the creator process

    if warned > 0:
        print(f"  [warn] {warned} labels produced no output pixels (removed by morphological op or too small)")

    return smoothed


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Smooth label segmentation masks to reduce polygon complexity."
    )
    parser.add_argument("input_mask",  help="Input label mask TIFF.")
    parser.add_argument("output_mask", help="Output smoothed mask TIFF.")
    parser.add_argument(
        "--kernel-size", type=int, default=2,
        help="Radius of the disk structuring element for morphological "
             "close+open (default: 2)."
    )
    parser.add_argument(
        "--n-workers", type=int, default=8,
        help="Number of parallel worker processes (default: 8). "
             "Match to --cpus-per-task in your SLURM script."
    )
    parser.add_argument(
        "--chunksize", type=int, default=200,
        help="Labels per worker chunk (default: 200)."
    )
    args = parser.parse_args()

    print(f"Reading mask: {args.input_mask}")
    mask = tifffile.imread(args.input_mask)
    original_dtype = mask.dtype
    print(f"Mask shape: {mask.shape}, dtype: {original_dtype}")

    n_labels = len(np.unique(mask)) - 1  # exclude background
    print(f"Number of labels: {n_labels}")

    if n_labels == 0:
        print("No labels found — writing input unchanged.")
        tifffile.imwrite(args.output_mask, mask)
        sys.exit(0)

    smoothed = smooth_label_morphological_parallel(
        mask,
        kernel_size=args.kernel_size,
        n_workers=args.n_workers,
        chunksize=args.chunksize,
    )

    # Preserve original dtype
    smoothed = smoothed.astype(original_dtype)

    n_after = len(np.unique(smoothed)) - 1
    print(f"Labels after smoothing: {n_after} (lost {n_labels - n_after})")

    print(f"Writing smoothed mask: {args.output_mask}")
    tifffile.imwrite(args.output_mask, smoothed)
    print("Done.")


if __name__ == "__main__":
    main()