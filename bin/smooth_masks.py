#!/usr/bin/env python3
"""
Smooth label segmentation masks using Shapely polygon simplification.
Parallelised for large masks (e.g. 40k x 40k with 100k+ labels).

Rasterizes in crop-local space to minimise memory overhead per worker.
"""

import argparse
import sys
import numpy as np
import tifffile
from multiprocessing import shared_memory
from concurrent.futures import ProcessPoolExecutor
from scipy.ndimage import find_objects
from skimage.measure import find_contours
from skimage.draw import polygon as draw_polygon
from shapely.geometry import Polygon
from shapely.validation import make_valid


# ── Worker (must be module-level for pickling) ────────────────────────────────

def _process_label(args):
    """
    Process a single label: contour → simplify → rasterize.
    Rasterizes in crop-local space then offsets to global coordinates.
    """
    (label_id, slc_starts, slc_stops, shm_name,
     mask_shape, mask_dtype, tolerance, min_area) = args

    shm = shared_memory.SharedMemory(name=shm_name)
    label_mask = np.ndarray(mask_shape, dtype=mask_dtype, buffer=shm.buf)

    try:
        slc = tuple(slice(s, e) for s, e in zip(slc_starts, slc_stops))
        crop = label_mask[slc] == label_id

        if not np.any(crop):
            return label_id, None, None

        contours = find_contours(crop, level=0.5)
        if not contours:
            return label_id, None, None

        largest = max(contours, key=lambda c: c.shape[0])
        if len(largest) < 4:
            return label_id, None, None

        minr, minc = slc_starts
        crop_h = slc_stops[0] - slc_starts[0]
        crop_w = slc_stops[1] - slc_starts[1]
        crop_shape = (crop_h, crop_w)

        # Contour coords are crop-local (row, col) → convert to (x, y) for Shapely
        coords = [(c[1], c[0]) for c in largest]

        poly = Polygon(coords)
        if not poly.is_valid:
            poly = make_valid(poly)
        if poly.is_empty or poly.area < min_area:
            return label_id, None, None

        simplified = poly.simplify(tolerance=tolerance, preserve_topology=True)
        if simplified.is_empty:
            return label_id, None, None

        polys = (
            [simplified] if simplified.geom_type == "Polygon"
            else list(simplified.geoms)
        )

        all_rr, all_cc = [], []
        for p in polys:
            if p.is_empty or p.area < min_area:
                continue

            x, y = p.exterior.coords.xy

            # Rasterize in crop-local space
            local_y = np.array(y)  # already crop-local
            local_x = np.array(x)  # already crop-local
            rr, cc = draw_polygon(local_y, local_x, shape=crop_shape)

            # Offset to global coordinates
            rr = rr + minr
            cc = cc + minc

            # Clip to mask bounds — guards against sub-pixel edge cases
            valid = (
                (rr >= 0) & (rr < mask_shape[0]) &
                (cc >= 0) & (cc < mask_shape[1])
            )
            all_rr.append(rr[valid])
            all_cc.append(cc[valid])

        if not all_rr:
            return label_id, None, None

        return label_id, np.concatenate(all_rr), np.concatenate(all_cc)

    except Exception as e:
        print(f"  [warn] label {label_id} failed: {e}", flush=True)
        # Fallback: return original pixels for this label
        orig_rr, orig_cc = np.where(label_mask == label_id)
        return label_id, orig_rr, orig_cc

    finally:
        shm.close()  # Never unlink in workers — only the creator does that


# ── Main smoothing function ───────────────────────────────────────────────────

def smooth_label_shapely_parallel(
    label_mask, tolerance=1.0, min_area=10, n_workers=8, chunksize=200
):
    """
    Parallelised Shapely-based label smoothing using shared memory.

    Parameters
    ----------
    label_mask : np.ndarray
        Integer label image where 0 is background.
    tolerance : float
        Douglas-Peucker tolerance in pixels. Higher = more simplification.
        Typical values: 0.5–2.0.
    min_area : int
        Minimum polygon area (px²) to retain after simplification.
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
            tolerance, min_area
        ))

    smoothed = np.zeros_like(label_mask)
    done = 0

    try:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            for label_id, rr, cc in executor.map(
                _process_label, work, chunksize=chunksize
            ):
                if rr is not None and len(rr) > 0:
                    # First-write wins for any overlap between adjacent labels
                    free = smoothed[rr, cc] == 0
                    smoothed[rr[free], cc[free]] = label_id
                done += 1
                if done % 5000 == 0:
                    print(f"  {done}/{n_labels} labels done...", flush=True)
    finally:
        shm.close()
        shm.unlink()  # Only unlink in the creator process

    return smoothed


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Parallelised Shapely smoothing of label masks for QuPath export."
    )
    parser.add_argument("input_mask",  help="Input label mask TIFF.")
    parser.add_argument("output_mask", help="Output smoothed mask TIFF.")
    parser.add_argument(
        "--tolerance", type=float, default=1.0,
        help="Douglas-Peucker tolerance in pixels (default: 1.0)."
    )
    parser.add_argument(
        "--min-area", type=int, default=10,
        help="Minimum polygon area in px² to retain (default: 10)."
    )
    parser.add_argument(
        "--n-workers", type=int, default=8,
        help="Number of parallel worker processes (default: 8). "
             "Match to --cpus-per-task in your SLURM script."
    )
    parser.add_argument(
        "--chunksize", type=int, default=200,
        help="Labels per worker chunk (default: 200). "
             "Increase for many small cells, decrease for fewer large cells."
    )
    args = parser.parse_args()

    print(f"Reading: {args.input_mask}")
    mask = tifffile.imread(args.input_mask)
    original_dtype = mask.dtype
    print(f"Shape: {mask.shape}, dtype: {original_dtype}")

    n_labels = len(np.unique(mask)) - 1  # exclude background
    print(f"Labels found: {n_labels}")

    if n_labels == 0:
        print("No labels found — writing input unchanged.")
        tifffile.imwrite(args.output_mask, mask)
        sys.exit(0)

    smoothed = smooth_label_shapely_parallel(
        mask,
        tolerance=args.tolerance,
        min_area=args.min_area,
        n_workers=args.n_workers,
        chunksize=args.chunksize,
    )

    # Preserve original dtype
    smoothed = smoothed.astype(original_dtype)

    n_after = len(np.unique(smoothed)) - 1
    print(f"Labels after smoothing: {n_after} (lost {n_labels - n_after})")

    print(f"Writing: {args.output_mask}")
    tifffile.imwrite(args.output_mask, smoothed)
    print("Done.")


if __name__ == "__main__":
    main()