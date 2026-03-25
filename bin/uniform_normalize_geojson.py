#!/usr/bin/env python3
"""
UniFORM-style measurement normalization for cellmeasurement GeoJSON outputs.

This script applies a cohort-level distribution alignment to numeric cellmeasurement
features across samples and writes per-sample `*_uniform.geojson` outputs.
It intentionally excludes KRONOS-derived measurements.
"""

import argparse
import json
import math
import os
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.signal import correlate


def parse_args():
    parser = argparse.ArgumentParser(description="Normalize cellmeasurement GeoJSON measurements across samples")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input GeoJSON files")
    parser.add_argument("--num-bins", type=int, default=1024, help="Histogram bin count")
    parser.add_argument("--min-value", type=float, default=1.0, help="Minimum value before log transform")
    parser.add_argument(
        "--exclude-pattern",
        default=r"^(kronos_|emb_)",
        help="Regex for measurement keys to exclude from normalization",
    )
    parser.add_argument(
        "--output-suffix",
        default="_uniform",
        help="Suffix to append before .geojson for normalized outputs",
    )
    return parser.parse_args()


def is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and np.isfinite(value)


def log_transform(values, min_value=1.0):
    arr = np.asarray(values, dtype=float)
    valid = arr >= min_value
    out = np.full(arr.shape, np.nan, dtype=float)
    out[valid] = np.log(arr[valid])
    return out


def choose_reference(hist_matrix):
    mean_hist = np.mean(hist_matrix, axis=0)
    distances = np.linalg.norm(hist_matrix - mean_hist, axis=1)
    return int(np.argmin(distances))


def compute_fft_shifts(reference_hist, hist_list):
    shifts = []
    for hist in hist_list:
        corr_fft = correlate(hist.flatten(), reference_hist.flatten(), mode="full", method="fft")
        shift_fft = np.argmax(corr_fft) - (len(reference_hist) - 1)
        shifts.append(int(shift_fft))
    return shifts


def load_geojsons(paths):
    records = []
    for path in paths:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        sample_id = Path(path).stem
        cell_features = [
            (idx, feat)
            for idx, feat in enumerate(data.get("features", []))
            if feat.get("properties", {}).get("objectType") == "cell"
        ]
        records.append({
            "path": path,
            "sample_id": sample_id,
            "data": data,
            "cell_features": cell_features,
        })
    return records


def collect_measurement_keys(records, exclude_pattern):
    regex = re.compile(exclude_pattern)
    key_counts = defaultdict(int)
    sample_count = len(records)

    for rec in records:
        keys_here = set()
        for _, feat in rec["cell_features"]:
            measurements = feat.get("properties", {}).get("measurements", {})
            for key, value in measurements.items():
                if regex.search(key):
                    continue
                if is_number(value):
                    keys_here.add(key)
        for key in keys_here:
            key_counts[key] += 1

    # Use keys present in all samples for robust cohort alignment
    return sorted([key for key, count in key_counts.items() if count == sample_count])


def normalize_records(records, keys, num_bins, min_value):
    n_samples = len(records)
    if n_samples < 2:
        print("Only one sample found; writing pass-through normalized files.")
        return {key: [1.0] * n_samples for key in keys}

    scales_by_key = {}

    for key in keys:
        sample_values = []
        global_min = float("inf")
        global_max = float("-inf")

        for rec in records:
            vals = []
            for _, feat in rec["cell_features"]:
                measurements = feat.get("properties", {}).get("measurements", {})
                value = measurements.get(key)
                if is_number(value):
                    vals.append(float(value))
            log_vals = log_transform(vals, min_value=min_value)
            log_vals = log_vals[np.isfinite(log_vals)]
            if log_vals.size == 0:
                log_vals = np.array([0.0], dtype=float)
            sample_values.append(log_vals)
            global_min = min(global_min, float(np.min(log_vals)))
            global_max = max(global_max, float(np.max(log_vals)))

        if not np.isfinite(global_min) or not np.isfinite(global_max) or global_max <= global_min:
            scales_by_key[key] = [1.0] * n_samples
            continue

        hist_list = []
        for log_vals in sample_values:
            counts, _ = np.histogram(log_vals, bins=num_bins, range=(global_min, global_max))
            hist_list.append(counts.astype(float))

        hist_matrix = np.stack(hist_list, axis=0)
        ref_idx = choose_reference(hist_matrix)
        shifts = compute_fft_shifts(hist_matrix[ref_idx], hist_list)

        increment = (global_max - global_min) / max(1, (num_bins - 1))
        scales = [math.exp(-shift * increment) for shift in shifts]
        scales_by_key[key] = scales

    # Apply scaling in-place
    for s_idx, rec in enumerate(records):
        for _, feat in rec["cell_features"]:
            measurements = feat.get("properties", {}).get("measurements", {})
            for key in keys:
                value = measurements.get(key)
                if is_number(value):
                    measurements[key] = float(value) * scales_by_key[key][s_idx]

    return scales_by_key


def write_outputs(records, output_suffix):
    output_paths = []
    for rec in records:
        in_path = Path(rec["path"])
        out_path = in_path.with_name(f"{in_path.stem}{output_suffix}.geojson")
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(rec["data"], f)
        output_paths.append(str(out_path))
    return output_paths


def main():
    args = parse_args()
    records = load_geojsons(args.inputs)

    if not records:
        raise SystemExit("No input GeoJSON files provided.")

    keys = collect_measurement_keys(records, args.exclude_pattern)
    print(f"Loaded {len(records)} GeoJSON files")
    print(f"Normalizing {len(keys)} shared numeric measurement keys")

    if keys:
        normalize_records(records, keys, num_bins=args.num_bins, min_value=args.min_value)

    outputs = write_outputs(records, args.output_suffix)
    for out in outputs:
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
