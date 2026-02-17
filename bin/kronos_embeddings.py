#!/usr/bin/env python3
"""
KRONOS Embedding Extraction for Segmented Cells

This script extracts KRONOS foundation model embeddings for each segmented cell
in a multiplexed image. It takes the original TIFF image and whole-cell segmentation
mask, extracts cell-centered patches, matches marker channels to KRONOS marker
metadata (including fuzzy matching for unresolved markers), runs the KRONOS model,
and outputs per-cell embeddings as a CSV file.

Usage:
    python kronos_embeddings.py \
        --tiff image.tiff \
        --mask whole_cell_mask.tiff \
        --model-path /path/to/kronos_model.pt \
        --marker-metadata /path/to/marker_metadata.csv \
        --output embeddings.csv \
        [--patch-size 64] \
        [--batch-size 32] \
        [--max-value 65535] \
        [--marker-mapping '{"CD3e": "CD3E"}'] \
        [--geojson annotations.geojson] \
        [--merge-geojson]
"""

import argparse
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd
import torch
import torch.utils.data
import tifffile

from difflib import SequenceMatcher

warnings.filterwarnings("ignore")


def get_channel_names(tif_path):
    """
    Extract channel/marker names from a TIFF file's metadata.

    Tries OME-XML metadata first, then ImageJ metadata, then falls back
    to numbered channel names.

    Args:
        tif_path: Path to the TIFF file.

    Returns:
        List of channel name strings.
    """
    with tifffile.TiffFile(tif_path) as tif:
        img = tif.asarray()
        if img.ndim == 2:
            n_channels = 1
        elif img.ndim == 3:
            n_channels = img.shape[0]
        else:
            n_channels = img.shape[0]

        channel_names = []

        # Try OME-XML metadata
        if tif.is_ome:
            try:
                from xml.etree import ElementTree as ET
                ome_xml = tif.ome_metadata
                root = ET.fromstring(ome_xml)
                ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
                channels = root.findall(".//ome:Channel", ns)
                if not channels:
                    channels = root.findall(".//Channel")
                channel_names = [ch.get("Name") or ch.get("ID") for ch in channels]
            except Exception:
                pass

        # Try ImageJ metadata
        if not channel_names and hasattr(tif, "imagej_metadata") and tif.imagej_metadata:
            if "Labels" in tif.imagej_metadata:
                channel_names = tif.imagej_metadata["Labels"]

        # Fallback to numbered channels
        if not channel_names:
            channel_names = [f"Channel_{i}" for i in range(n_channels)]

    return channel_names


def match_markers(image_marker_names, kronos_metadata_df, user_mapping=None, top_suggestions=5):
    """
    Match image channel marker names to KRONOS marker metadata.

    Uses exact matching (case-insensitive), user-provided mappings, and
    fuzzy string matching (SequenceMatcher) for suggestions.

    Args:
        image_marker_names: List of marker names from the image channels.
        kronos_metadata_df: DataFrame of KRONOS marker metadata (marker_name, marker_id, marker_mean, marker_std).
        user_mapping: Optional dict mapping image marker names to KRONOS marker names.
        top_suggestions: Number of fuzzy-match suggestions to report for unmatched markers.

    Returns:
        Tuple of (matched_df, unmatched_report):
            - matched_df: DataFrame with columns [image_channel, marker_name, marker_id, marker_mean, marker_std]
              for all matched markers.
            - unmatched_report: String describing unmatched markers and suggestions.
    """
    # Build lookup by uppercase name
    kronos_metadata_df = kronos_metadata_df.copy()
    kronos_upper_lookup = {}
    for _, row in kronos_metadata_df.iterrows():
        kronos_upper_lookup[row["marker_name"].upper()] = row

    kronos_names_upper = list(kronos_upper_lookup.keys())

    matched_records = []
    unmatched_markers = []

    for ch_idx, marker_name in enumerate(image_marker_names):
        resolved_name = marker_name

        # Check user mapping first
        if user_mapping and marker_name in user_mapping:
            resolved_name = user_mapping[marker_name]

        upper_name = resolved_name.upper()

        if upper_name in kronos_upper_lookup:
            row = kronos_upper_lookup[upper_name]
            matched_records.append({
                "image_channel": ch_idx,
                "marker_name": marker_name,
                "kronos_marker_name": row["marker_name"],
                "marker_id": int(row["marker_id"]),
                "marker_mean": float(row["marker_mean"]),
                "marker_std": float(row["marker_std"]),
            })
        else:
            unmatched_markers.append((ch_idx, marker_name))

    # Build unmatched report with fuzzy suggestions
    report_lines = []
    if unmatched_markers:
        report_lines.append(f"WARNING: {len(unmatched_markers)} marker(s) could not be matched to KRONOS metadata:")
        report_lines.append(f"  (Matched {len(matched_records)} of {len(image_marker_names)} markers)")
        report_lines.append("")
        for ch_idx, marker_name in unmatched_markers:
            similarities = [
                (kn, SequenceMatcher(None, marker_name.upper(), kn).ratio())
                for kn in kronos_names_upper
            ]
            similarities.sort(key=lambda x: x[1], reverse=True)
            suggestions = [s[0] for s in similarities[:top_suggestions]]
            report_lines.append(f"  Channel {ch_idx}: '{marker_name}' -> suggestions: {suggestions}")
        report_lines.append("")
        report_lines.append("  To resolve, provide --marker-mapping as JSON, e.g.:")
        report_lines.append('    --marker-mapping \'{"MyMarker": "KRONOS_MARKER_NAME"}\'')

    matched_df = pd.DataFrame(matched_records)
    unmatched_report = "\n".join(report_lines)

    return matched_df, unmatched_report


def extract_cell_patches(image, mask, patch_size=64):
    """
    Extract cell-centered patches from the multiplexed image using the segmentation mask.

    For each unique cell ID in the mask, finds the centroid and extracts a
    (num_channels, patch_size, patch_size) patch centered on it. Patches near
    image boundaries are zero-padded.

    Args:
        image: numpy array of shape (C, H, W) - the multiplexed image.
        mask: numpy array of shape (H, W) - the whole-cell segmentation mask.
        patch_size: Size of the square patch to extract around each cell centroid.

    Returns:
        Tuple of (cell_ids, centroids, patches):
            - cell_ids: list of integer cell IDs
            - centroids: list of (y_center, x_center) tuples
            - patches: numpy array of shape (num_cells, C, patch_size, patch_size)
    """
    cell_ids_unique = np.unique(mask)
    cell_ids_unique = cell_ids_unique[cell_ids_unique > 0]  # skip background

    C, H, W = image.shape
    half = patch_size // 2

    cell_ids = []
    centroids = []
    patches = []

    for cell_id in cell_ids_unique:
        # Find centroid of cell
        ys, xs = np.where(mask == cell_id)
        y_center = int(np.mean(ys))
        x_center = int(np.mean(xs))

        # Calculate patch bounds
        y1 = y_center - half
        y2 = y_center + half
        x1 = x_center - half
        x2 = x_center + half

        # Create padded patch
        patch = np.zeros((C, patch_size, patch_size), dtype=image.dtype)

        # Calculate source and destination slices accounting for boundaries
        src_y1 = max(0, y1)
        src_y2 = min(H, y2)
        src_x1 = max(0, x1)
        src_x2 = min(W, x2)

        dst_y1 = src_y1 - y1
        dst_y2 = dst_y1 + (src_y2 - src_y1)
        dst_x1 = src_x1 - x1
        dst_x2 = dst_x1 + (src_x2 - src_x1)

        patch[:, dst_y1:dst_y2, dst_x1:dst_x2] = image[:, src_y1:src_y2, src_x1:src_x2]

        cell_ids.append(int(cell_id))
        centroids.append((y_center, x_center))
        patches.append(patch)

    if len(patches) > 0:
        patches = np.stack(patches, axis=0)
    else:
        patches = np.zeros((0, C, patch_size, patch_size), dtype=image.dtype)

    return cell_ids, centroids, patches


def load_kronos_model(model_path, config_path=None, device="cpu"):
    """
    Load the KRONOS model from a local checkpoint.

    Args:
        model_path: Path to the .pt model weights file.
        config_path: Path to config.json. If None, searches in the same directory as model_path.
        device: torch device string.

    Returns:
        Tuple of (model, precision, embedding_dim).
    """
    # Try to find config.json alongside the model
    if config_path is None:
        model_dir = os.path.dirname(model_path)
        candidate = os.path.join(model_dir, "config.json")
        if os.path.exists(candidate):
            config_path = candidate

    if config_path and os.path.exists(config_path):
        with open(config_path) as f:
            cfg = json.load(f)
    else:
        cfg = {"model_type": "vits16", "token_overlap": False}

    # Import kronos - it should be installed in the container/environment
    from kronos import create_model_from_pretrained

    model, precision, embedding_dim = create_model_from_pretrained(
        checkpoint_path=model_path,
        cfg_path=config_path,
        cfg=cfg,
    )

    model = model.to(device)
    model.eval()

    return model, precision, embedding_dim


def extract_embeddings(model, patches, marker_ids, marker_means, marker_stds,
                       max_value=65535.0, batch_size=32, num_workers=4,
                       device="cpu", precision=torch.float32):
    """
    Run KRONOS inference on cell patches to extract embeddings.

    Uses a DataLoader with configurable num_workers for parallel data
    loading and preprocessing, which significantly speeds up inference
    on HPC systems with many CPU cores.

    Args:
        model: Loaded KRONOS model.
        patches: numpy array of shape (N, C, H, W) with raw intensity values.
        marker_ids: numpy array of marker IDs for each channel (length C).
        marker_means: numpy array of marker means for normalization (length C).
        marker_stds: numpy array of marker stds for normalization (length C).
        max_value: Maximum intensity value for initial normalization.
        batch_size: Batch size for inference.
        num_workers: Number of DataLoader workers for parallel data loading.
        device: torch device string.
        precision: torch dtype for model input.

    Returns:
        Dict with keys:
            - 'patch_embeddings': numpy array of shape (N, embedding_dim)
            - 'marker_embeddings': numpy array of shape (N, num_markers, embedding_dim) or None
    """
    N = patches.shape[0]
    if N == 0:
        return {"patch_embeddings": np.array([]), "marker_embeddings": None}

    # Create a simple TensorDataset and DataLoader for parallel loading
    means = torch.tensor(marker_means, dtype=precision)
    stds = torch.tensor(marker_stds, dtype=precision)
    marker_ids_tensor = torch.tensor(marker_ids, dtype=torch.long)

    class _PatchDataset(torch.utils.data.Dataset):
        def __init__(self, patches_np, means, stds, max_value, precision):
            self.patches_np = patches_np
            self.means = means
            self.stds = stds
            self.max_value = max_value
            self.precision = precision

        def __len__(self):
            return len(self.patches_np)

        def __getitem__(self, idx):
            patch = torch.tensor(self.patches_np[idx], dtype=self.precision)
            # Normalize: scale by max value, then standardize
            patch = patch / self.max_value
            patch = (patch - self.means[:, None, None]) / self.stds[:, None, None]
            return patch

    dataset = _PatchDataset(patches, means, stds, max_value, precision)
    dataloader = torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=(device != "cpu"),
        prefetch_factor=2 if num_workers > 0 else None,
    )

    all_patch_embeddings = []
    all_marker_embeddings = []
    processed = 0
    total_batches = len(dataloader)

    for batch_idx, batch in enumerate(dataloader):
        batch = batch.to(device, non_blocking=True)
        batch_marker_ids = marker_ids_tensor.unsqueeze(0).expand(batch.shape[0], -1).to(device)

        with torch.no_grad():
            output = model(batch, marker_ids=batch_marker_ids)

        # Handle different output formats
        if isinstance(output, tuple):
            patch_emb = output[0].cpu().numpy()  # patch/cls embeddings
            marker_emb = output[1].cpu().numpy() if len(output) > 1 else None
        elif isinstance(output, dict):
            patch_emb = output.get("x_norm_clstoken", output.get("cls_token")).cpu().numpy()
            marker_emb = output.get("x_norm_patchtokens", output.get("patch_tokens"))
            if marker_emb is not None:
                marker_emb = marker_emb.cpu().numpy()
        else:
            patch_emb = output.cpu().numpy()
            marker_emb = None

        all_patch_embeddings.append(patch_emb)
        if marker_emb is not None:
            all_marker_embeddings.append(marker_emb)

        processed += batch.shape[0]
        print(f"  Processed {processed}/{N} cells (batch {batch_idx+1}/{total_batches})", end="\r")

    print()  # newline after progress

    result = {
        "patch_embeddings": np.concatenate(all_patch_embeddings, axis=0),
        "marker_embeddings": np.concatenate(all_marker_embeddings, axis=0) if all_marker_embeddings else None,
    }

    return result


def save_embeddings(cell_ids, centroids, embeddings, output_path, sample_id=None):
    """
    Save per-cell embeddings to a CSV file.

    Columns: cell_id, y_center, x_center, emb_0, emb_1, ..., emb_N

    Args:
        cell_ids: List of cell IDs.
        centroids: List of (y, x) centroid tuples.
        embeddings: numpy array of shape (num_cells, embedding_dim).
        output_path: Path for the output CSV file.
        sample_id: Optional sample identifier to include.
    """
    if len(cell_ids) == 0:
        # Write empty file with headers only
        with open(output_path, "w") as f:
            f.write("cell_id,y_center,x_center\n")
        print(f"No cells found. Empty file written to {output_path}")
        return

    embedding_dim = embeddings.shape[1]
    emb_columns = [f"emb_{i}" for i in range(embedding_dim)]

    data = {
        "cell_id": cell_ids,
        "y_center": [c[0] for c in centroids],
        "x_center": [c[1] for c in centroids],
    }

    if sample_id:
        data["sample_id"] = [sample_id] * len(cell_ids)

    for i, col in enumerate(emb_columns):
        data[col] = embeddings[:, i]

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    print(f"Saved {len(cell_ids)} cell embeddings ({embedding_dim}D) to {output_path}")


def merge_embeddings_into_geojson(geojson_path, cell_ids, centroids, embeddings, output_path, distance_threshold=5.0):
    """
    Merge KRONOS embeddings into a GeoJSON file by matching cells via centroid proximity.

    The cellmeasurement GeoJSON uses UUIDs (not mask labels) for cell IDs, so we
    match cells by finding the nearest GeoJSON feature centroid to each embedding
    centroid. A distance threshold prevents incorrect matches.

    Args:
        geojson_path: Path to the input GeoJSON from cellmeasurement.
        cell_ids: List of cell IDs (mask labels) from embedding extraction.
        centroids: List of (y, x) centroid tuples from embedding extraction.
        embeddings: numpy array of shape (num_cells, embedding_dim).
        output_path: Path for the output merged GeoJSON.
        distance_threshold: Maximum distance (pixels) for centroid matching.

    Returns:
        Tuple of (num_matched, num_unmatched, num_geojson_cells).
    """
    with open(geojson_path) as f:
        geojson = json.load(f)

    # Extract centroids from GeoJSON cell features using Shoelace formula
    cell_features = []
    cell_geojson_centroids = []
    for feature in geojson.get("features", []):
        props = feature.get("properties", {})
        if props.get("objectType") != "cell":
            continue
        cell_features.append(feature)
        # Compute proper area-weighted centroid using Shoelace formula
        # GeoJSON coordinates are [x, y] where x=column, y=row in image space
        coords = feature["geometry"]["coordinates"][0]  # outer ring
        n = len(coords) - 1  # last point == first point (closed ring)
        A = 0.0
        cx = 0.0
        cy = 0.0
        for i in range(n):
            x0, y0 = coords[i]
            x1, y1 = coords[i + 1]
            cross = x0 * y1 - x1 * y0
            A += cross
            cx += (x0 + x1) * cross
            cy += (y0 + y1) * cross
        A = A / 2.0
        if abs(A) > 1e-10:
            cx = cx / (6.0 * A)
            cy = cy / (6.0 * A)
        else:
            # Degenerate polygon, fall back to vertex average
            cx = sum(c[0] for c in coords) / len(coords)
            cy = sum(c[1] for c in coords) / len(coords)
        cell_geojson_centroids.append((cy, cx))  # (y, x) to match embedding convention

    if not cell_features or len(cell_ids) == 0:
        # Nothing to merge, write unchanged
        with open(output_path, "w") as f:
            json.dump(geojson, f)
        return 0, 0, len(cell_features)

    # Build arrays for vectorised distance computation
    emb_centroids = np.array(centroids)  # (N, 2) as (y, x)
    geo_centroids = np.array(cell_geojson_centroids)  # (M, 2) as (y, x)

    embedding_dim = embeddings.shape[1]
    emb_columns = [f"kronos_emb_{i}" for i in range(embedding_dim)]

    num_matched = 0
    num_unmatched = 0

    # For each GeoJSON cell, find nearest embedding centroid
    for feat_idx, feature in enumerate(cell_features):
        gc = geo_centroids[feat_idx]
        # Compute distances to all embedding centroids
        dists = np.sqrt(np.sum((emb_centroids - gc) ** 2, axis=1))
        best_idx = np.argmin(dists)
        best_dist = dists[best_idx]

        if best_dist <= distance_threshold:
            # Add embeddings as measurements
            measurements = feature["properties"].setdefault("measurements", {})
            measurements["kronos_cell_id"] = int(cell_ids[best_idx])
            for i, col_name in enumerate(emb_columns):
                measurements[col_name] = float(embeddings[best_idx, i])
            num_matched += 1
        else:
            num_unmatched += 1

    with open(output_path, "w") as f:
        json.dump(geojson, f)

    return num_matched, num_unmatched, len(cell_features)


def main():
    parser = argparse.ArgumentParser(
        description="Extract KRONOS embeddings for segmented cells"
    )
    parser.add_argument("--tiff", required=True, help="Path to multiplexed TIFF image")
    parser.add_argument("--mask", required=True, help="Path to whole-cell segmentation mask TIFF")
    parser.add_argument("--model-path", required=True, help="Path to KRONOS model weights (.pt)")
    parser.add_argument("--config-path", default=None, help="Path to KRONOS config.json (auto-detected if not specified)")
    parser.add_argument("--marker-metadata", required=True, help="Path to KRONOS marker_metadata.csv")
    parser.add_argument("--output", required=True, help="Output CSV path for cell embeddings")
    parser.add_argument("--patch-size", type=int, default=64, help="Patch size for cell-centered crops (default: 64)")
    parser.add_argument("--batch-size", type=int, default=32, help="Batch size for model inference (default: 32)")
    parser.add_argument("--num-workers", type=int, default=4, help="Number of DataLoader workers for parallel data loading (default: 4)")
    parser.add_argument("--max-value", type=float, default=65535.0, help="Maximum intensity value for normalization (default: 65535)")
    parser.add_argument("--marker-mapping", default=None, help="JSON string or file path mapping image markers to KRONOS markers")
    parser.add_argument("--sample-id", default=None, help="Sample identifier to include in output")
    parser.add_argument("--geojson", default=None, help="Path to GeoJSON from cellmeasurement (required for --merge-geojson)")
    parser.add_argument("--merge-geojson", action="store_true", default=False,
                        help="Merge embeddings into GeoJSON as per-cell properties")
    parser.add_argument("--distance-threshold", type=float, default=5.0,
                        help="Max pixel distance for centroid matching between embeddings and GeoJSON cells (default: 5.0)")
    args = parser.parse_args()

    if args.merge_geojson and not args.geojson:
        parser.error("--merge-geojson requires --geojson to be specified")

    # Determine device
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")
    
    # Disable multiprocessing workers on CPU or NFS to avoid cleanup issues
    if device == "cpu" and args.num_workers > 0:
        print(f"  Warning: Disabling num_workers (was {args.num_workers}) on CPU to avoid NFS/multiprocessing issues")
        args.num_workers = 0

    # Load marker metadata
    print(f"Loading KRONOS marker metadata from {args.marker_metadata}")
    kronos_metadata = pd.read_csv(args.marker_metadata)

    # Parse user marker mapping
    user_mapping = None
    if args.marker_mapping:
        if os.path.isfile(args.marker_mapping):
            with open(args.marker_mapping) as f:
                user_mapping = json.load(f)
        else:
            try:
                user_mapping = json.loads(args.marker_mapping)
            except json.JSONDecodeError:
                print(f"WARNING: Could not parse --marker-mapping as JSON: {args.marker_mapping}")

    # Read image
    print(f"Reading image: {args.tiff}")
    image = tifffile.imread(args.tiff)
    if image.ndim == 2:
        image = image[np.newaxis, ...]  # add channel dim
    print(f"  Image shape: {image.shape} (C={image.shape[0]}, H={image.shape[1]}, W={image.shape[2]})")

    # Get channel names from image metadata
    channel_names = get_channel_names(args.tiff)
    print(f"  Found {len(channel_names)} channel names: {channel_names[:10]}{'...' if len(channel_names) > 10 else ''}")

    # Match markers
    matched_df, unmatched_report = match_markers(channel_names, kronos_metadata, user_mapping)
    if unmatched_report:
        print(unmatched_report)
    print(f"  Matched {len(matched_df)} markers to KRONOS metadata")

    if len(matched_df) == 0:
        print("ERROR: No markers could be matched to KRONOS metadata. Exiting.")
        sys.exit(1)

    # Select only matched channels from the image
    matched_channels = matched_df["image_channel"].values
    image_matched = image[matched_channels]
    marker_ids = matched_df["marker_id"].values
    marker_means = matched_df["marker_mean"].values
    marker_stds = matched_df["marker_std"].values

    print(f"  Using {len(matched_channels)} channels for embedding extraction")

    # Read segmentation mask
    print(f"Reading segmentation mask: {args.mask}")
    mask = tifffile.imread(args.mask)
    if mask.ndim == 3:
        mask = mask[0]  # take first channel if multi-channel
    print(f"  Mask shape: {mask.shape}, unique cells: {len(np.unique(mask)) - 1}")

    # Extract cell-centered patches
    print(f"Extracting cell patches (patch_size={args.patch_size})...")
    cell_ids, centroids, patches = extract_cell_patches(image_matched, mask, args.patch_size)
    print(f"  Extracted {len(cell_ids)} cell patches")

    if len(cell_ids) == 0:
        print("WARNING: No cells found in segmentation mask.")
        save_embeddings([], [], np.array([]), args.output, args.sample_id)
        return

    # Load KRONOS model
    print(f"Loading KRONOS model from {args.model_path}")
    model, precision, embedding_dim = load_kronos_model(args.model_path, args.config_path, device)
    print(f"  Model loaded: precision={precision}, embedding_dim={embedding_dim}")

    # Extract embeddings
    print(f"Extracting embeddings (batch_size={args.batch_size}, num_workers={args.num_workers})...")
    results = extract_embeddings(
        model, patches, marker_ids, marker_means, marker_stds,
        max_value=args.max_value, batch_size=args.batch_size,
        num_workers=args.num_workers, device=device, precision=precision,
    )

    patch_embeddings = results["patch_embeddings"]
    print(f"  Patch embeddings shape: {patch_embeddings.shape}")

    # Save embeddings
    save_embeddings(cell_ids, centroids, patch_embeddings, args.output, args.sample_id)

    # Save marker matching report
    report_path = args.output.replace(".csv", "_marker_report.txt")
    with open(report_path, "w") as f:
        f.write("KRONOS Marker Matching Report\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Image: {args.tiff}\n")
        f.write(f"Total image channels: {len(channel_names)}\n")
        f.write(f"Matched markers: {len(matched_df)}\n\n")
        f.write("Matched markers:\n")
        for _, row in matched_df.iterrows():
            f.write(f"  Ch {row['image_channel']}: {row['marker_name']} -> {row['kronos_marker_name']} (ID={row['marker_id']})\n")
        if unmatched_report:
            f.write(f"\n{unmatched_report}\n")
    print(f"  Marker matching report saved to {report_path}")

    # Optionally merge embeddings into GeoJSON
    if args.merge_geojson:
        merged_geojson_path = args.output.replace("_kronos_embeddings.csv", "_kronos_merged.geojson")
        if not merged_geojson_path.endswith(".geojson"):
            merged_geojson_path = args.output.rsplit(".", 1)[0] + "_kronos_merged.geojson"
        print(f"Merging embeddings into GeoJSON: {args.geojson}")
        num_matched, num_unmatched, num_geo_cells = merge_embeddings_into_geojson(
            args.geojson, cell_ids, centroids, patch_embeddings,
            merged_geojson_path, distance_threshold=args.distance_threshold
        )
        print(f"  GeoJSON merge: {num_matched}/{num_geo_cells} cells matched, {num_unmatched} unmatched")
        print(f"  Merged GeoJSON saved to {merged_geojson_path}")

    print("Done!")


if __name__ == "__main__":
    main()
