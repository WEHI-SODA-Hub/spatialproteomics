#!/usr/bin/env python3
"""
KRONOS2 per-cell embedding extraction.

Extracts a 768-d KRONOS2 CLS embedding for every cell in a CELLMEASUREMENT
GeoJSON and writes them back into that GeoJSON as per-cell measurements.

Everything except the two model calls lives in ``kronos2_common``; this script
is the KRONOS2-specific part:

    normed = model.preprocess(patches, markers, preferred_dapi=...)
    cls    = model(tensor(normed), markers)              # (B, 768)

The marker-aware z-score lives *inside* the model, so this pipeline must not
reimplement it -- it only supplies patches already divided by the image dtype's
scaling factor, which is the domain ``preprocess`` expects.

Usage:
    kronos2_embeddings.py \
        --tiff image.ome.tiff \
        --geojson cells.geojson.gz \
        --model-path /path/to/KRONOS2 \
        --output-geojson sample.geojson.gz \
        [--patch-size 64] [--batch-size 16] \
        [--marker-mapping map.json] [--nuclear-marker DAPI-01] \
        [--no-isolate-cell]
"""

import argparse
import logging
import os
import sys
import time
import warnings

import numpy as np

from kronos2_common import (
    ChannelReader,
    CellPatchDataset,
    apply_marker_mapping,
    get_channel_names,
    load_marker_mapping,
    log,
    numpy_collate,
    read_cells,
    scaling_factor,
    write_geojson_with_embeddings,
    write_marker_report,
)

warnings.filterwarnings("ignore")


def load_model(model_path, device=None):
    """Load KRONOS2 via AutoModel(trust_remote_code=True).

    The vendored DINOv2 + xFormers stack emits a caught Triton-probe traceback
    and an "xFormers is available" warning at import. Both look like failures
    and are not, so they are silenced here -- output only, no effect on the
    model or its embeddings. The xformers logger is raised to ERROR so a
    genuine error still surfaces.
    """
    logging.getLogger("xformers").setLevel(logging.ERROR)

    import torch
    from transformers import AutoModel

    dev = device or ("cuda" if torch.cuda.is_available() else "cpu")
    log(f"Loading KRONOS2 from {model_path} onto {dev}")
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="xFormers is available")
        model = AutoModel.from_pretrained(
            str(model_path), trust_remote_code=True, device=dev
        )
    log(f"  model ready in {time.time() - t0:.1f}s")
    return model, dev


def novel_markers(model, markers, nuclear_marker=None):
    """Markers absent from KRONOS2's vocabulary, using the model's own matcher.

    Reproduces the resolution inside ``model.preprocess`` -- ``clean_marker_name``
    plus ``marker_match_key`` against ``_marker_index``, with the DAPI
    special-case -- rather than reimplementing it, so this can never disagree
    with what the model actually does. Such markers are not an error: they fall
    back to a default (mean, std), which is why they are reported loudly.
    """
    stats = getattr(model, "_marker_stats", None) or {}
    if not stats:
        return []
    module = sys.modules[type(model).__module__]
    clean = module.clean_marker_name
    match_key = module.marker_match_key
    index = getattr(model, "_marker_index", None) or module.build_marker_index(stats)

    dapi = clean(nuclear_marker) if nuclear_marker else None
    out = []
    for marker in markers:
        if dapi is not None and clean(marker) == dapi:
            continue
        if match_key(marker) not in index:
            out.append(marker)
    return out


def extract(model, dataset, markers, nuclear_marker, batch_size, num_workers, device):
    """Run KRONOS2 over every cell patch, returning ``(n_cells, 768)``."""
    import torch

    n = len(dataset)
    if n == 0:
        return np.zeros((0, 768), dtype=np.float32)

    try:
        from torch.utils.data import DataLoader

        loader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=False,
            num_workers=num_workers,
            collate_fn=numpy_collate,
        )
    except Exception:
        loader = (
            numpy_collate([dataset[i] for i in range(s, min(s + batch_size, n))])
            for s in range(0, n, batch_size)
        )

    out = np.zeros((n, model.config.embed_dim), dtype=np.float32)
    done = 0
    t0 = time.time()
    total_batches = (n + batch_size - 1) // batch_size

    for batch_idx, (patches, indices) in enumerate(loader):
        # The z-score is the model's; we only hand it dtype-scaled patches.
        normed = model.preprocess(patches, markers, preferred_dapi=nuclear_marker)
        with torch.inference_mode():
            tensor = torch.from_numpy(np.ascontiguousarray(normed)).to(device)
            cls = model(tensor, markers)
        out[np.asarray(indices)] = cls.float().cpu().numpy().astype(np.float32)

        done += len(indices)
        if batch_idx == 0 or (batch_idx + 1) % max(1, total_batches // 10) == 0:
            rate = done / max(time.time() - t0, 1e-6)
            log(f"  {done}/{n} cells ({rate:.0f} cells/s)")

    log(f"  inference complete: {n} cells in {time.time() - t0:.1f}s")
    return out


def main():
    parser = argparse.ArgumentParser(description="KRONOS2 per-cell embeddings")
    parser.add_argument("--tiff", required=True, help="Multiplexed TIFF")
    parser.add_argument(
        "--geojson", required=True, help="CELLMEASUREMENT GeoJSON defining the cells"
    )
    parser.add_argument(
        "--model-path", required=True, help="KRONOS2 snapshot directory or HF repo id"
    )
    parser.add_argument("--output-geojson", required=True, help="Output GeoJSON path")
    parser.add_argument("--marker-report", default=None, help="Marker report path")
    parser.add_argument("--patch-size", type=int, default=64)
    parser.add_argument(
        "--batch-size",
        type=int,
        default=16,
        help="Batch size (16 reproduces the published gold standard; fp32 "
        "cuBLAS kernels are batch-dependent, so this shifts results slightly)",
    )
    parser.add_argument("--num-workers", type=int, default=4)
    parser.add_argument(
        "--max-value",
        type=float,
        default=None,
        help="Override the intensity divisor; default derives it from the "
        "image dtype (uint8=255, uint16=65535, float=400)",
    )
    parser.add_argument("--marker-mapping", default=None)
    parser.add_argument(
        "--nuclear-marker",
        default=None,
        help="Channel to treat as the nuclear stain; the model's only alias "
        "mechanism, so a DRAQ5/Hoechst slide needs this to get DAPI stats",
    )
    parser.add_argument("--isolate-cell", dest="isolate", action="store_true")
    parser.add_argument("--no-isolate-cell", dest="isolate", action="store_false")
    parser.set_defaults(isolate=True)
    parser.add_argument(
        "--allow-novel-defaults",
        action="store_true",
        help="Proceed when markers are outside KRONOS2's vocabulary (they take "
        "default normalisation stats). Without this the run fails loudly.",
    )
    args = parser.parse_args()

    log("=" * 60)
    log("KRONOS2 Embedding Extraction")
    log("=" * 60)
    for key, value in vars(args).items():
        log(f"    {key}: {value}")

    # --- channels + marker names ------------------------------------------
    channel_names = get_channel_names(args.tiff)
    log(f"Found {len(channel_names)} channels")
    mapping = load_marker_mapping(args.marker_mapping)
    markers, applied = apply_marker_mapping(channel_names, mapping)
    if applied:
        log(f"  applied {len(applied)} marker mapping(s)")

    # --- cells -------------------------------------------------------------
    log(f"Reading cells from {args.geojson}")
    geojson, features, centroids, rings = read_cells(args.geojson)
    log(f"  {len(features)} cells")

    # --- image -------------------------------------------------------------
    reader = ChannelReader(args.tiff, list(range(len(channel_names))))
    divisor = args.max_value if args.max_value else scaling_factor(reader.dtype)
    log(
        f"  image {reader.n_channels}x{reader.height}x{reader.width} "
        f"dtype={reader.dtype} divisor={divisor} access={reader.mode}"
    )

    # --- model -------------------------------------------------------------
    model, device = load_model(args.model_path)

    unseen = novel_markers(model, markers, args.nuclear_marker)
    report_extra = []
    if unseen:
        msg = (
            f"{len(unseen)} of {len(markers)} markers are outside KRONOS2's "
            f"vocabulary and will use DEFAULT normalisation stats: {unseen}"
        )
        report_extra.append("Novel markers (default stats used):\n  " + "\n  ".join(unseen))
        if args.allow_novel_defaults:
            log(f"WARNING: {msg}")
        else:
            log(f"ERROR: {msg}")
            log(
                "  Most novel markers are naming variants, not new biology -- "
                "map them onto vocabulary names with --marker-mapping. Re-run "
                "with --allow-novel-defaults to accept default stats instead."
            )
            sys.exit(1)

    if args.marker_report:
        write_marker_report(
            args.marker_report, args.tiff, channel_names, markers, applied, report_extra
        )

    # --- extract -----------------------------------------------------------
    if not features:
        log("No cells found; writing GeoJSON unchanged")
        write_geojson_with_embeddings(
            geojson, [], np.zeros((0, 1), dtype=np.float32), args.output_geojson
        )
        return

    dataset = CellPatchDataset(
        reader, centroids, rings, args.patch_size, divisor, isolate=args.isolate
    )
    log(
        f"Extracting {len(dataset)} cell patches "
        f"(patch={args.patch_size}, isolate={args.isolate}, batch={args.batch_size})"
    )
    embeddings = extract(
        model,
        dataset,
        markers,
        args.nuclear_marker,
        args.batch_size,
        args.num_workers if device != "cpu" else 0,
        device,
    )
    log(f"  embeddings: {embeddings.shape}")

    write_geojson_with_embeddings(geojson, features, embeddings, args.output_geojson)
    reader.close()
    log("=" * 60)
    log("KRONOS2 Embedding Extraction - Done")
    log("=" * 60)


if __name__ == "__main__":
    main()
