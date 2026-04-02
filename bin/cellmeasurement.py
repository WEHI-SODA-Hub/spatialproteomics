#!/usr/bin/env python3
"""
Python implementation of the cellmeasurement Groovy app.

This script measures cell and nucleus compartments from labeled masks and a
multi-channel TIFF image, then exports a GeoJSON FeatureCollection.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, mapping, shape
from shapely.ops import unary_union
from skimage.measure import block_reduce, find_contours, regionprops
from skimage.morphology import binary_erosion, disk
from skimage.segmentation import watershed

_DISK_1 = disk(1)
_DISK_1.flags.writeable = False


@dataclass
class CellRecord:
    cell_id: int
    cell_label: Optional[int]
    nucleus_label: Optional[int]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract cell measurements from masks and image TIFF.")
    p.add_argument("-n", "--nuclear-mask", required=True, help="Nuclear segmentation mask TIFF")
    p.add_argument("-w", "--whole-cell-mask", required=True, help="Whole-cell segmentation mask TIFF")
    p.add_argument("-f", "--tiff-file", required=True, help="Multi-channel TIFF image")
    p.add_argument("-o", "--output-file", required=True, help="Output GeoJSON path")
    p.add_argument("-d", "--downsample-factor", type=float, default=1.0)
    p.add_argument("-p", "--pixel-size-microns", type=float, default=0.5)
    p.add_argument("--skip-measurements", action="store_true")
    p.add_argument("--simplify-rois", action="store_true",
                   help="Simplify ROI geometry with Douglas-Peucker. Enabled by default; use --no-simplify-rois to disable.")
    p.add_argument("--no-simplify-rois", dest="simplify_rois", action="store_false")
    p.set_defaults(simplify_rois=True)
    p.add_argument("--tolerance", type=float, default=1.4,
                   help="Simplification tolerance in pixels (default 1.4, matching Cellpose extension VW distance)")
    p.add_argument("--percentiles", default="")
    p.add_argument("--erosion-steps", default="")
    p.add_argument("--expansion-steps", default="",
                   help="CSV of positive ints: measure intensity in annular rings dilated outward from cell boundary")
    p.add_argument("-i", "--dist-threshold", type=float, default=10.0)
    p.add_argument("-e", "--estimate-cell-boundary-dist", type=float, default=3.0)
    p.add_argument("-t", "--threads", type=int, default=1)
    p.add_argument("--tile-size", type=int, default=2048)
    p.add_argument("--tile-overlap", type=int, default=200)
    p.add_argument("--pretty-json", action="store_true", help="Write indented GeoJSON output")
    p.add_argument("--output-mask", default="",
                   help="Write a rasterized label mask TIFF from the final cell geometries")
    return p.parse_args()


def tile_flags_explicit(argv: Sequence[str]) -> bool:
    return any(
        a == "--tile-size"
        or a.startswith("--tile-size=")
        or a == "--tile-overlap"
        or a.startswith("--tile-overlap=")
        for a in argv
    )


def load_label_mask(path: str) -> np.ndarray:
    arr = tifffile.imread(path)
    if arr.ndim != 2:
        raise ValueError(f"Mask must be 2D label image, got shape={arr.shape} for {path}")
    # Normalize to signed integer labels to keep downstream ops consistent.
    return arr.astype(np.int64, copy=False)


def load_image(path: str) -> Tuple[np.ndarray, List[str]]:
    img = tifffile.imread(path)
    ch_names: List[str] = []

    with tifffile.TiffFile(path) as tf:
        # Strategy 1: OME-XML Channel/@Name (covers OME-TIFF, OPAL QPTIFF, COMET)
        try:
            import xml.etree.ElementTree as ET

            ome = tf.ome_metadata
            # Fallback: OPAL QPTIFF stores OME-XML in first page ImageDescription
            if not ome and tf.pages:
                first_desc = tf.pages[0].description
                if isinstance(first_desc, str) and first_desc.strip().startswith("<"):
                    ome = first_desc
            if ome:
                root = ET.fromstring(ome)
                ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
                channels = root.findall(".//ome:Channel", ns)
                if channels:
                    ch_names = [
                        ch.get("Name") or ch.get("ID") or ""
                        for ch in channels
                    ]
        except Exception as e:
            print(f"Warning: failed to parse OME metadata for channel names ({e})")
            ch_names = []

        # Strategy 2: MIBI JSON metadata in per-page ImageDescription
        if not ch_names:
            try:
                import json as _json
                first_desc = tf.pages[0].description if tf.pages else ""
                _json.loads(first_desc)  # probe for JSON
                ch_names = [
                    str(_json.loads(p.description).get("channel.target", ""))
                    for p in tf.pages
                ]
            except (ValueError, TypeError, KeyError, AttributeError):
                pass

        # Strategy 3: ImageJ metadata Labels
        if not ch_names:
            try:
                ij = tf.imagej_metadata
                if ij and "Labels" in ij:
                    ch_names = [str(lbl) for lbl in ij["Labels"]]
            except Exception:
                pass

    if img.ndim == 2:
        img = img[np.newaxis, ...]
    elif img.ndim == 3:
        # Accept C,Y,X or Y,X,C. Heuristic: if first dim is small, treat as channels.
        if img.shape[0] <= 64:
            pass
        elif img.shape[2] <= 64:
            img = np.transpose(img, (2, 0, 1))
        else:
            raise ValueError(f"Unsupported 3D image layout: {img.shape}")
    else:
        raise ValueError(f"Unsupported image dimensions: {img.shape}")

    if not ch_names or len(ch_names) != img.shape[0]:
        print(f"Warning: found {len(ch_names)} channel names but image has {img.shape[0]} channels; using fallback names")
        ch_names = [f"Channel {i + 1}" for i in range(img.shape[0])]
    else:
        print(f"Detected channel names: {ch_names}")

    return img.astype(np.float32, copy=False), ch_names


def maybe_downsample(image_cyx: np.ndarray, nuc: np.ndarray, whole: np.ndarray, ds: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if ds <= 1.0:
        return image_cyx, nuc, whole
    step = int(round(ds))
    if step < 2:
        return image_cyx, nuc, whole
    # Crop to divisible shape first to avoid block_reduce zero-padding bias at edges.
    c, h, w = image_cyx.shape
    h2 = (h // step) * step
    w2 = (w // step) * step
    if h2 == 0 or w2 == 0:
        return image_cyx, nuc, whole

    image_crop = image_cyx[:, :h2, :w2]
    nuc_crop = nuc[:h2, :w2]
    whole_crop = whole[:h2, :w2]

    # Use reshape+mean for image intensities (faster for this regular block case)
    # and block max for label masks.
    image_ds = image_crop.reshape(c, h2 // step, step, w2 // step, step).mean(axis=(2, 4))
    nuc_ds = block_reduce(nuc_crop, block_size=(step, step), func=np.max)
    whole_ds = block_reduce(whole_crop, block_size=(step, step), func=np.max)
    return image_ds.astype(np.float32, copy=False), nuc_ds.astype(np.int64, copy=False), whole_ds.astype(np.int64, copy=False)


def label_props_dict(label_img: np.ndarray) -> Dict[int, Dict[str, Any]]:
    out: Dict[int, Dict[str, Any]] = {}
    for r in regionprops(label_img):
        out[int(r.label)] = {
            "centroid": (float(r.centroid[0]), float(r.centroid[1])),
            "bbox": (int(r.bbox[0]), int(r.bbox[1]), int(r.bbox[2]), int(r.bbox[3])),
        }
    return out


def match_cells(
    nuc: np.ndarray, whole: np.ndarray, dist_threshold: float, estimate_dist: float
) -> Tuple[np.ndarray, np.ndarray, List[CellRecord], Dict[str, int], Dict[int, Tuple[slice, slice]]]:
    whole_props = label_props_dict(whole)
    nuc_props = label_props_dict(nuc)

    records: List[CellRecord] = []
    bbox_map: Dict[int, Tuple[slice, slice]] = {}
    out_cell = np.zeros_like(whole, dtype=np.int64)
    out_nuc = np.zeros_like(nuc, dtype=np.int64)

    whole_labels = sorted(whole_props.keys())
    whole_pts = np.array([whole_props[l]["centroid"] for l in whole_labels], dtype=np.float64)
    tree = cKDTree(whole_pts) if len(whole_pts) else None

    next_id = 1
    used_whole: set[int] = set()
    dropped_synth_cells = 0

    # Track which nuclei need synthesized boundaries for the watershed pass.
    unmatched_nuclei: List[Tuple[int, int]] = []  # (nuc_label, assigned_cell_id)

    # --- Pass 1: match nuclei to whole-cell masks; defer unmatched nuclei ---
    for nlab in sorted(nuc_props.keys()):
        nr, nc = nuc_props[nlab]["centroid"]
        nminr, nminc, nmaxr, nmaxc = nuc_props[nlab]["bbox"]
        matched_whole = None

        if tree is not None:
            dist, idx = tree.query([nr, nc], k=1, distance_upper_bound=dist_threshold)
            if np.isfinite(dist) and idx < len(whole_labels):
                candidate = whole_labels[int(idx)]
                if candidate not in used_whole:
                    matched_whole = candidate

        if matched_whole is not None:
            cminr, cminc, cmaxr, cmaxc = whole_props[matched_whole]["bbox"]
            minr = min(nminr, cminr)
            minc = min(nminc, cminc)
            maxr = max(nmaxr, cmaxr)
            maxc = max(nmaxc, cmaxc)
            used_whole.add(matched_whole)

            rs = slice(minr, maxr)
            cs = slice(minc, maxc)
            npatch = nuc[rs, cs] == nlab
            cpatch = whole[rs, cs] == matched_whole

            # Avoid overlaps in synthesized output labels.
            available = out_cell[rs, cs] == 0
            cpatch = cpatch & available
            npatch = npatch & cpatch

            if not np.any(cpatch):
                dropped_synth_cells += 1
                continue

            out_cell[rs, cs][cpatch] = next_id
            out_nuc[rs, cs][npatch] = next_id

            cell_rows, cell_cols = np.nonzero(cpatch)
            cell_minr = minr + int(cell_rows.min())
            cell_minc = minc + int(cell_cols.min())
            cell_maxr = minr + int(cell_rows.max()) + 1
            cell_maxc = minc + int(cell_cols.max()) + 1
            bbox_map[next_id] = (
                slice(max(0, cell_minr - 1), min(nuc.shape[0], cell_maxr + 1)),
                slice(max(0, cell_minc - 1), min(nuc.shape[1], cell_maxc + 1)),
            )
            records.append(CellRecord(next_id, matched_whole, int(nlab)))
            next_id += 1
        else:
            # Reserve an id and record it; actual pixels assigned in pass 2.
            unmatched_nuclei.append((nlab, next_id))
            next_id += 1

    # --- Pass 2: watershed-partition unmatched nuclei so they don't overlap ---
    if unmatched_nuclei:
        # Build a seed image: each unmatched nucleus gets its reserved cell id.
        seeds = np.zeros_like(out_cell)
        for nlab, cid in unmatched_nuclei:
            seeds[nuc == nlab] = cid

        # Restrict growth to pixels within estimate_dist of any unmatched nucleus
        # and not already claimed by pass-1 cells.
        any_unmatched_nuc = seeds > 0
        growth_zone = ndi.binary_dilation(
            any_unmatched_nuc,
            structure=disk(max(1, int(round(estimate_dist)))),
        )
        growth_zone = growth_zone & (out_cell == 0)

        # Euclidean distance from each background pixel to the nearest seed;
        # watershed floods lowest-distance pixels first → Voronoi-like partition.
        dist_map = ndi.distance_transform_edt(seeds == 0)
        ws = watershed(dist_map, markers=seeds, mask=growth_zone)

        for nlab, cid in unmatched_nuclei:
            cpatch_full = ws == cid
            npatch_full = (nuc == nlab) & cpatch_full

            if not np.any(cpatch_full):
                dropped_synth_cells += 1
                continue

            out_cell[cpatch_full] = cid
            out_nuc[npatch_full] = cid

            cell_rows, cell_cols = np.nonzero(cpatch_full)
            cell_minr = int(cell_rows.min())
            cell_minc = int(cell_cols.min())
            cell_maxr = int(cell_rows.max()) + 1
            cell_maxc = int(cell_cols.max()) + 1
            bbox_map[cid] = (
                slice(max(0, cell_minr - 1), min(nuc.shape[0], cell_maxr + 1)),
                slice(max(0, cell_minc - 1), min(nuc.shape[1], cell_maxc + 1)),
            )
            records.append(CellRecord(cid, None, int(nlab)))

    unmatched_whole = len(set(whole_labels) - used_whole)
    stats = {
        "nucleus_count": len(nuc_props),
        "whole_cell_count": len(whole_props),
        "matched_cells": len(records),
        "unmatched_whole_cells": unmatched_whole,
        "dropped_synth_cells": dropped_synth_cells,
    }
    return out_cell, out_nuc, records, stats, bbox_map


def mask_to_geometry(
    mask: np.ndarray,
    simplify: bool,
    tolerance: float,
    row_offset: int = 0,
    col_offset: int = 0,
):
    if not np.any(mask):
        return None

    contours = find_contours(mask.astype(np.uint8), level=0.5)
    if not contours:
        return None

    polys = []
    for c in contours:
        if len(c) < 3:
            continue
        xy = [(float(col_offset + p[1]), float(row_offset + p[0])) for p in c]
        poly = Polygon(xy)
        if poly.is_empty:
            continue
        if not poly.is_valid:
            poly = poly.buffer(0)
        if not poly.is_empty:
            polys.append(poly)

    if not polys:
        return None

    g = unary_union(polys)
    if g.geom_type == "GeometryCollection":
        keep = [geom for geom in g.geoms if geom.geom_type in ("Polygon", "MultiPolygon") and not geom.is_empty]
        if not keep:
            return None
        g = unary_union(keep)

    if simplify and tolerance > 0:
        g = g.simplify(tolerance, preserve_topology=True)
    if not g.is_valid:
        g = g.buffer(0)
    if g.is_empty:
        return None

    return g


def basic_shape_metrics(cell_mask: np.ndarray, nuc_mask: Optional[np.ndarray], px_um: float) -> Dict[str, float]:
    rp = regionprops(cell_mask.astype(np.uint8))
    if not rp:
        return {}
    r = rp[0]
    area_um2 = float(r.area) * (px_um ** 2)
    perimeter_um = float(r.perimeter) * px_um if r.perimeter > 0 else 0.0
    circularity = float(4 * math.pi * r.area / (r.perimeter ** 2)) if r.perimeter > 0 else 0.0
    maj = float(r.major_axis_length) * px_um
    mino = float(r.minor_axis_length) * px_um
    solidity = float(r.solidity) if r.solidity is not None else 0.0

    out = {
        "Cell: Area µm^2": area_um2,
        "Cell: Circularity": circularity,
        "Cell: Length µm": perimeter_um,
        "Cell: Max diameter µm": maj,
        "Cell: Min diameter µm": mino,
        "Cell: Solidity": solidity,
    }

    if nuc_mask is not None and np.any(nuc_mask):
        nrp = regionprops(nuc_mask.astype(np.uint8))
        if nrp:
            nr = nrp[0]
            n_area = float(nr.area) * (px_um ** 2)
            out["Nucleus: Area µm^2"] = n_area
            out["Nucleus: Circularity"] = float(4 * math.pi * nr.area / (nr.perimeter ** 2)) if nr.perimeter > 0 else 0.0
            out["Nucleus: Length µm"] = float(nr.perimeter) * px_um if nr.perimeter > 0 else 0.0
            out["Nucleus: Max diameter µm"] = float(nr.major_axis_length) * px_um
            out["Nucleus: Min diameter µm"] = float(nr.minor_axis_length) * px_um
            out["Nucleus: Solidity"] = float(nr.solidity) if nr.solidity is not None else 0.0
            out["Nucleus/Cell area ratio"] = float(n_area / area_um2) if area_um2 > 0 else 0.0

    return out


def compartment_masks(cell_mask: np.ndarray, nuc_mask: np.ndarray) -> Dict[str, np.ndarray]:
    cm = cell_mask.astype(bool)
    nm = nuc_mask.astype(bool)
    cyto = cm & ~nm
    mem = cm & ~binary_erosion(cm, _DISK_1)
    return {
        "CELL": cm,
        "NUCLEUS": nm,
        "CYTOPLASM": cyto,
        "MEMBRANE": mem,
    }


def stat_values(vals: np.ndarray) -> Dict[str, float]:
    if vals.size == 0:
        return {}
    return {
        "Mean": float(np.mean(vals)),
        "Median": float(np.median(vals)),
        "Min": float(np.min(vals)),
        "Max": float(np.max(vals)),
        "Std.Dev.": float(np.std(vals)),
    }


def add_intensity_measurements(props: Dict[str, Any], image_cyx: np.ndarray, ch_names: Sequence[str], comp_masks: Dict[str, np.ndarray]):
    labels = {"CELL": "Cell", "NUCLEUS": "Nucleus", "CYTOPLASM": "Cytoplasm", "MEMBRANE": "Membrane"}
    for ci, ch in enumerate(ch_names):
        ch_img = image_cyx[ci]
        for comp, m in comp_masks.items():
            vals = ch_img[m]
            if vals.size == 0:
                continue
            for k, v in stat_values(vals).items():
                props[f"{ch}: {labels[comp]}: {k}"] = v


def add_percentiles(props: Dict[str, Any], image_cyx: np.ndarray, ch_names: Sequence[str], comp_masks: Dict[str, np.ndarray], percentiles: Sequence[float]):
    labels = {"CELL": "Cell", "NUCLEUS": "Nucleus", "CYTOPLASM": "Cytoplasm", "MEMBRANE": "Membrane"}
    for ci, ch in enumerate(ch_names):
        ch_img = image_cyx[ci]
        for comp, m in comp_masks.items():
            vals = ch_img[m]
            if vals.size == 0:
                continue
            for p in percentiles:
                props[f"{ch}: {labels[comp]}: Percentile: {p}"] = float(np.percentile(vals, p))


def erode_steps(mask: np.ndarray, step: int) -> np.ndarray:
    if step <= 0:
        return mask
    return ndi.binary_erosion(mask, structure=_DISK_1, iterations=step)


def add_erosion_measurements(props: Dict[str, Any], image_cyx: np.ndarray, ch_names: Sequence[str], comp_masks: Dict[str, np.ndarray], steps: Sequence[int]):
    ordered_steps = sorted(set(int(s) for s in steps if int(s) > 0))
    if not ordered_steps:
        return

    for comp in ("CELL", "NUCLEUS"):
        base = comp_masks[comp]
        base_area = int(np.count_nonzero(base))
        if base_area == 0:
            continue
        prev = 0
        cur = base.copy()
        comp_name = comp.capitalize()
        for s in ordered_steps:
            extra = s - prev
            cur = erode_steps(cur, extra)
            prev = s
            area = int(np.count_nonzero(cur))
            props[f"{comp_name}: Eroded_{s}px: Area_Fraction"] = float(area / base_area)
            if area == 0:
                continue
            for ci, ch in enumerate(ch_names):
                vals = image_cyx[ci][cur]
                if vals.size == 0:
                    continue
                props[f"{ch}: {comp_name}: Eroded_{s}px: Mean"] = float(np.mean(vals))
                props[f"{ch}: {comp_name}: Eroded_{s}px: Median"] = float(np.median(vals))


def add_expansion_measurements(
    props: Dict[str, Any],
    image_cyx: np.ndarray,
    ch_names: Sequence[str],
    cell_mask: np.ndarray,
    steps: Sequence[int],
):
    """Measure intensity in annular rings dilated outward from the cell boundary.

    Each ring at step *s* covers the zone [dilated_{s-1}, dilated_s) so that
    rings are mutually exclusive and together form a radial profile of the
    pericellular neighbourhood.
    """
    ordered_steps = sorted(set(int(s) for s in steps if int(s) > 0))
    if not ordered_steps:
        return

    base_area = int(np.count_nonzero(cell_mask))
    if base_area == 0:
        return

    prev_dilated = cell_mask.astype(bool)
    prev_step = 0
    for s in ordered_steps:
        extra = s - prev_step
        cur_dilated = ndi.binary_dilation(prev_dilated, structure=_DISK_1, iterations=extra)
        ring = cur_dilated & ~prev_dilated
        ring_area = int(np.count_nonzero(ring))
        props[f"Cell: Expanded_{s}px: Area_Fraction"] = float(ring_area / base_area)

        if ring_area > 0:
            for ci, ch in enumerate(ch_names):
                vals = image_cyx[ci][ring]
                if vals.size == 0:
                    continue
                props[f"{ch}: Cell: Expanded_{s}px: Mean"] = float(np.mean(vals))
                props[f"{ch}: Cell: Expanded_{s}px: Median"] = float(np.median(vals))

        prev_dilated = cur_dilated
        prev_step = s


def parse_csv_numbers(s: str, cast=float, positive_only=False) -> List:
    if not s or not s.strip():
        return []
    out = []
    for x in s.split(","):
        x = x.strip()
        if not x:
            continue
        try:
            v = cast(x)
        except (TypeError, ValueError) as e:
            raise ValueError(f"Invalid value '{x}' in list '{s}'") from e
        if positive_only and v <= 0:
            continue
        out.append(v)
    return out


def feature_for_cell(
    cell_id: int,
    cell_crop: np.ndarray,
    nuc_crop: np.ndarray,
    image_crop: np.ndarray,
    row_offset: int,
    col_offset: int,
    rec_cell_label: Optional[int],
    rec_nucleus_label: Optional[int],
    simplify_rois: bool,
    tolerance: float,
    pixel_size_microns: float,
    skip_measurements: bool,
    ch_names: Sequence[str],
    percentiles: Sequence[float],
    erosion_steps: Sequence[int],
    expansion_steps: Sequence[int] = (),
):
    cmask = cell_crop == cell_id
    nmask = nuc_crop == cell_id
    geom = mask_to_geometry(cmask, simplify_rois, tolerance, row_offset=row_offset, col_offset=col_offset)
    if geom is None:
        return None
    nuc_geom = mask_to_geometry(nmask, simplify_rois, tolerance, row_offset=row_offset, col_offset=col_offset)

    measurements: Dict[str, Any] = {
        "id": int(cell_id),
        "cell_label": int(rec_cell_label) if rec_cell_label is not None else None,
        "nucleus_label": int(rec_nucleus_label) if rec_nucleus_label is not None else None,
    }
    measurements.update(basic_shape_metrics(cmask, nmask, pixel_size_microns))

    if not skip_measurements:
        comps = compartment_masks(cmask, nmask)
        add_intensity_measurements(measurements, image_crop, ch_names, comps)
        if percentiles:
            add_percentiles(measurements, image_crop, ch_names, comps, percentiles)
        if erosion_steps:
            add_erosion_measurements(measurements, image_crop, ch_names, comps, erosion_steps)
        if expansion_steps:
            add_expansion_measurements(measurements, image_crop, ch_names, cmask, expansion_steps)

    feature: Dict[str, Any] = {
        "type": "Feature",
        "id": f"cell-{cell_id}",
        "geometry": mapping(geom),
        "properties": {
            "objectType": "cell",
            "id": int(cell_id),
            "cell_label": int(rec_cell_label) if rec_cell_label is not None else None,
            "nucleus_label": int(rec_nucleus_label) if rec_nucleus_label is not None else None,
            "measurements": measurements,
        },
    }
    if nuc_geom is not None:
        feature["nucleusGeometry"] = mapping(nuc_geom)

    return feature


def iter_tasks(
    unique_cells: Sequence[int],
    bbox_map: Dict[int, Tuple[slice, slice]],
    records_by_id: Dict[int, CellRecord],
    cell_labels: np.ndarray,
    nuc_labels: np.ndarray,
    img_cyx: np.ndarray,
    args: argparse.Namespace,
    ch_names: Sequence[str],
    percentiles: Sequence[float],
    erosion_steps: Sequence[int],
    expansion_steps: Sequence[int] = (),
):
    max_expand = max(expansion_steps) if expansion_steps else 0
    h, w = cell_labels.shape[:2]
    for cid in unique_cells:
        rs, cs = bbox_map[cid]
        if max_expand > 0:
            r0 = max(rs.start - max_expand, 0)
            r1 = min(rs.stop + max_expand, h)
            c0 = max(cs.start - max_expand, 0)
            c1 = min(cs.stop + max_expand, w)
            rs_pad = slice(r0, r1)
            cs_pad = slice(c0, c1)
        else:
            rs_pad, cs_pad = rs, cs
        rec = records_by_id.get(cid)
        # Copy in parent process to keep each task payload self-contained for worker pickling.
        yield (
            cid,
            cell_labels[rs_pad, cs_pad].copy(),
            nuc_labels[rs_pad, cs_pad].copy(),
            img_cyx[:, rs_pad, cs_pad].copy(),
            rs_pad.start,
            cs_pad.start,
            rec.cell_label if rec else None,
            rec.nucleus_label if rec else None,
            args.simplify_rois,
            args.tolerance,
            args.pixel_size_microns,
            args.skip_measurements,
            tuple(ch_names),
            tuple(percentiles),
            tuple(erosion_steps),
            tuple(expansion_steps),
        )


def constrain_cell_overlaps(features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Clip overlapping cell geometries so no two cells share area.

    Equivalent to QuPath's ``CellTools.constrainCellOverlaps()``.
    For each pair of overlapping cells the overlap is assigned to the cell
    with the **smaller** area (the larger cell is trimmed), which preserves
    small cells and matches QuPath's heuristic.

    Uses a grid-based spatial hash for broad-phase collision detection,
    avoiding any Shapely STRtree version-dependent behaviour.
    """
    if not features:
        return features

    n = len(features)
    geoms = [shape(f["geometry"]) for f in features]
    areas = [g.area for g in geoms]

    # -- broad-phase: grid spatial hash based on geometry bounding boxes --
    bounds = [g.bounds for g in geoms]  # (minx, miny, maxx, maxy)
    all_minx = min(b[0] for b in bounds)
    all_miny = min(b[1] for b in bounds)
    all_maxx = max(b[2] for b in bounds)
    all_maxy = max(b[3] for b in bounds)
    span = max(all_maxx - all_minx, all_maxy - all_miny, 1.0)
    # Target ~sqrt(n) cells so each cell has ~1 geometry on average
    grid_size = max(span / max(int(n ** 0.5), 1), 1.0)

    grid: Dict[Tuple[int, int], List[int]] = {}
    for i, (minx, miny, maxx, maxy) in enumerate(bounds):
        gx0 = int((minx - all_minx) / grid_size)
        gy0 = int((miny - all_miny) / grid_size)
        gx1 = int((maxx - all_minx) / grid_size)
        gy1 = int((maxy - all_miny) / grid_size)
        for gx in range(gx0, gx1 + 1):
            for gy in range(gy0, gy1 + 1):
                grid.setdefault((gx, gy), []).append(i)

    # -- narrow-phase: pairwise intersection test and clipping --
    checked: set = set()
    clipped = 0
    for cell_list in grid.values():
        for ii in range(len(cell_list)):
            i = cell_list[ii]
            for jj in range(ii + 1, len(cell_list)):
                j = cell_list[jj]
                pair = (min(i, j), max(i, j))
                if pair in checked:
                    continue
                checked.add(pair)

                gi = geoms[i]
                gj = geoms[j]
                if gi.is_empty or gj.is_empty:
                    continue

                try:
                    if not gi.intersects(gj):
                        continue
                    intersection = gi.intersection(gj)
                    if intersection.is_empty or intersection.area < 1e-10:
                        continue
                except Exception:
                    continue

                # Trim the larger cell; keep the smaller cell intact.
                if areas[i] >= areas[j]:
                    gi = gi.difference(gj)
                    if not gi.is_valid:
                        gi = gi.buffer(0)
                    geoms[i] = gi
                    areas[i] = gi.area
                else:
                    gj = gj.difference(gi)
                    if not gj.is_valid:
                        gj = gj.buffer(0)
                    geoms[j] = gj
                    areas[j] = gj.area
                clipped += 1

    out = []
    for f, g in zip(features, geoms):
        if g.is_empty:
            continue
        f["geometry"] = mapping(g)
        # Also clip nucleusGeometry if it extends beyond the trimmed cell
        if "nucleusGeometry" in f:
            try:
                ng = shape(f["nucleusGeometry"])
                ng = ng.intersection(g)
                if not ng.is_empty:
                    f["nucleusGeometry"] = mapping(ng)
                else:
                    del f["nucleusGeometry"]
            except Exception:
                del f["nucleusGeometry"]
        out.append(f)

    print(f"Overlap constraint: checked {len(checked)} pairs, clipped {clipped}, removed {n - len(out)} empty cells")
    return out


def rasterize_features_to_mask(
    features: List[Dict[str, Any]], height: int, width: int
) -> np.ndarray:
    """Rasterize cell feature polygons back to an integer label mask TIFF.

    Each cell is rasterized with its ``properties.id`` as the label value.
    Produces the same format as smooth_masks output — a 2-D integer label
    image where 0 is background.
    """
    from skimage.draw import polygon as draw_polygon

    max_id = max(
        (f["properties"].get("id", 0) for f in features if f["properties"].get("objectType") == "cell"),
        default=0,
    )
    dtype = np.int32 if max_id < 2**31 else np.int64
    mask = np.zeros((height, width), dtype=dtype)

    for feat in features:
        if feat["properties"].get("objectType") != "cell":
            continue
        cell_id = feat["properties"].get("id", 0)
        if cell_id <= 0:
            continue

        geom = shape(feat["geometry"])
        if geom.is_empty:
            continue

        polys = []
        if geom.geom_type == "Polygon":
            polys = [geom]
        elif geom.geom_type == "MultiPolygon":
            polys = list(geom.geoms)

        for poly in polys:
            ext = np.array(poly.exterior.coords)
            rr, cc = draw_polygon(ext[:, 1], ext[:, 0], shape=(height, width))
            mask[rr, cc] = cell_id

    return mask


def regularize_mask(
    mask: np.ndarray,
    seed_centroids: Optional[Dict[int, Tuple[float, float]]] = None,
) -> np.ndarray:
    """Re-partition a label mask using watershed from seed points.

    When overlapping instance masks (e.g. from CellSAM) are flattened to a
    single label plane, the last-written label 'wins' overlap pixels, creating
    irregular boundaries.  Watershed from seed points redistributes those
    pixels so boundaries between adjacent cells are equidistant.

    Parameters
    ----------
    mask : np.ndarray
        Integer label image where 0 is background.
    seed_centroids : dict, optional
        Mapping of label -> (row, col) centroid to use as seeds. If None,
        centroids are computed from the mask itself. Providing external
        centroids (e.g. from the nuclear mask) avoids using corrupted
        centroids from a flattened overlapping mask.
    """
    props = label_props_dict(mask)
    if not props:
        return mask

    seeds = np.zeros_like(mask, dtype=mask.dtype)
    for lab in props:
        if seed_centroids and lab in seed_centroids:
            r, c = seed_centroids[lab]
        else:
            r, c = props[lab]["centroid"]
        r = max(0, min(int(round(r)), mask.shape[0] - 1))
        c = max(0, min(int(round(c)), mask.shape[1] - 1))
        seeds[r, c] = lab

    occupied = mask > 0
    dist = ndi.distance_transform_edt(seeds == 0)
    result = watershed(dist, markers=seeds, mask=occupied).astype(mask.dtype)
    return result


def main() -> int:
    args = parse_args()
    if args.tile_size <= 0:
        raise ValueError("tile-size must be > 0")
    if args.tile_overlap < 0:
        raise ValueError("tile-overlap must be >= 0")
    if args.pixel_size_microns <= 0:
        raise ValueError("pixel-size-microns must be > 0")
    if args.dist_threshold <= 0:
        raise ValueError("dist-threshold must be > 0")
    if args.downsample_factor <= 0:
        raise ValueError("downsample-factor must be > 0")
    step = int(round(args.downsample_factor))
    if args.downsample_factor > 1.0 and step < 2:
        print(
            f"Warning: downsample-factor {args.downsample_factor} rounds to step={step}; no downsampling will be applied"
        )
    if tile_flags_explicit(sys.argv[1:]):
        print("Warning: --tile-size/--tile-overlap are parsed for compatibility but not used in this Python implementation")

    whole = load_label_mask(args.whole_cell_mask)
    nuc = load_label_mask(args.nuclear_mask)
    img_cyx, ch_names = load_image(args.tiff_file)

    img_cyx, nuc, whole = maybe_downsample(img_cyx, nuc, whole, args.downsample_factor)

    if whole.shape != nuc.shape:
        raise ValueError(f"Mask shapes differ: whole={whole.shape}, nuclear={nuc.shape}")
    if img_cyx.shape[1:] != whole.shape:
        raise ValueError(f"Image shape {img_cyx.shape[1:]} does not match mask shape {whole.shape}")

    print(f"Loaded whole cell mask: {whole.shape}")
    print(f"Loaded nuclear mask: {nuc.shape}")

    cell_labels, nuc_labels, records, match_stats, bbox_map = match_cells(
        nuc, whole, args.dist_threshold, args.estimate_cell_boundary_dist
    )

    unique_cells = [int(x) for x in np.unique(cell_labels) if x > 0]
    bbox_ids = set(bbox_map.keys())
    if len(unique_cells) != len(bbox_ids) or set(unique_cells) != bbox_ids:
        raise RuntimeError(
            "Internal mismatch between labeled cells and bbox map keys; "
            f"labels={len(unique_cells)}, bboxes={len(bbox_ids)}"
        )
    print(f"Total path objects: {len(unique_cells)}")
    print(
        "Matching summary: "
        f"nuclei={match_stats['nucleus_count']}, whole_cells={match_stats['whole_cell_count']}, "
        f"matched={match_stats['matched_cells']}, unmatched_whole={match_stats['unmatched_whole_cells']}, "
        f"dropped_estimated={match_stats['dropped_synth_cells']}"
    )

    records_by_id = {r.cell_id: r for r in records}

    percentiles = parse_csv_numbers(args.percentiles, cast=float)
    erosion_steps_raw = parse_csv_numbers(args.erosion_steps, cast=int, positive_only=True)
    erosion_steps = sorted(set(erosion_steps_raw))
    if erosion_steps_raw and erosion_steps_raw != erosion_steps:
        print(f"Warning: normalized erosion steps from {erosion_steps_raw} to {erosion_steps}")

    expansion_steps_raw = parse_csv_numbers(args.expansion_steps, cast=int, positive_only=True)
    expansion_steps = sorted(set(expansion_steps_raw))
    if expansion_steps_raw and expansion_steps_raw != expansion_steps:
        print(f"Warning: normalized expansion steps from {expansion_steps_raw} to {expansion_steps}")

    if percentiles:
        print(f"Will add intensity percentiles: {percentiles}")
    if erosion_steps:
        print(f"Will add erosion measurements at steps: {erosion_steps}")
    if expansion_steps:
        print(f"Will add expansion measurements at steps: {expansion_steps}")

    features = []
    total = len(unique_cells)
    task_iter = iter_tasks(
        unique_cells,
        bbox_map,
        records_by_id,
        cell_labels,
        nuc_labels,
        img_cyx,
        args,
        ch_names,
        percentiles,
        erosion_steps,
        expansion_steps,
    )

    if args.threads > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            max_inflight = max(1, args.threads * 4)
            future_map = {}
            for _ in range(max_inflight):
                try:
                    task = next(task_iter)
                except StopIteration:
                    break
                future_map[ex.submit(feature_for_cell, *task)] = task[0]

            done = 0
            while future_map:
                completed, _ = wait(future_map.keys(), return_when=FIRST_COMPLETED)
                for fut in completed:
                    cid = future_map.pop(fut)
                    try:
                        feat = fut.result()
                    except Exception as exc:
                        print(f"Warning: cell {cid} failed: {exc}", file=sys.stderr)
                        feat = None

                    if feat is not None:
                        features.append(feat)
                    done += 1
                    if done % 1000 == 0 or done == total:
                        print(f"Progress: {done}/{total} cells ({int(done * 100.0 / max(total, 1))}%)")

                while len(future_map) < max_inflight:
                    try:
                        task = next(task_iter)
                    except StopIteration:
                        break
                    future_map[ex.submit(feature_for_cell, *task)] = task[0]
    else:
        for i, t in enumerate(task_iter, 1):
            feat = feature_for_cell(*t)
            if feat is not None:
                features.append(feat)
            if i % 1000 == 0 or i == total:
                print(f"Progress: {i}/{total} cells ({int(i * 100.0 / max(total, 1))}%)")

    # Keep output deterministic regardless of parallel execution completion order.
    features.sort(key=lambda f: f["properties"].get("id", -1))

    # Resolve overlapping cell geometries (equivalent to QuPath CellTools.constrainCellOverlaps)
    features = constrain_cell_overlaps(features)

    # Top-level annotation feature for whole image extent.
    h, w = whole.shape
    annotation = {
        "type": "Feature",
        "id": "annotation-whole-image",
        "geometry": mapping(Polygon([(0, 0), (w, 0), (w, h), (0, h), (0, 0)])),
        "properties": {
            "objectType": "annotation",
            "type": "annotation",
            "name": "whole_image",
        },
    }

    out = {
        "type": "FeatureCollection",
        "features": [annotation] + features,
    }

    out_path = Path(args.output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        if args.pretty_json:
            json.dump(out, f, indent=2)
        else:
            json.dump(out, f, separators=(",", ":"))

    print(f"Exported to GeoJSON: {out_path}")

    if args.output_mask:
        mask_out = rasterize_features_to_mask(features, h, w)
        mask_path = Path(args.output_mask)
        mask_path.parent.mkdir(parents=True, exist_ok=True)
        tifffile.imwrite(str(mask_path), mask_out)
        print(f"Exported label mask: {mask_path} ({mask_out.dtype}, {len(np.unique(mask_out)) - 1} labels)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
