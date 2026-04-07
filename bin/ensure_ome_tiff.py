#!/usr/bin/env python3
"""
Pre-processing script that ensures a TIFF file has valid OME-XML metadata.

If the input TIFF already has OME-XML (proper OME-TIFF), it is left unchanged.
Otherwise (e.g. ImageJ TIFFs), channel names are extracted from available
metadata (ImageJ labels, page descriptions, etc.) and the file is rewritten
as a proper OME-TIFF so that downstream tools like mesmer-segment can parse it.

Usage:
    ensure_ome_tiff.py INPUT_TIFF OUTPUT_TIFF
"""
import sys
import json
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import tifffile


def _local_tag_name(tag):
    """Return XML tag name without namespace."""
    return tag.split("}", 1)[-1]


def has_valid_ome_xml(tif):
    """Check whether the TIFF already contains parseable OME-XML with channel info."""
    try:
        ome_xml = tif.ome_metadata
        if ome_xml:
            root = ET.fromstring(ome_xml)
            # Check there is at least one Channel element
            for elem in root.iter():
                if _local_tag_name(elem.tag) == "Channel":
                    return True
    except Exception:
        pass

    # Also try the first-page description (mesmer-segment reads this)
    if tif.pages:
        desc = tif.pages[0].description
        if desc:
            try:
                root = ET.fromstring(desc)
                for elem in root.iter():
                    if _local_tag_name(elem.tag) == "Channel":
                        return True
            except ET.ParseError:
                pass
    return False


def get_channel_names(tif, n_channels):
    """
    Extract channel names from TIFF metadata.
    Mirrors the logic in cellsam_segment.py.

    Tries:
      1. MIBI JSON metadata (per-page JSON with channel.target)
      2. OME-XML metadata
      3. ImageJ metadata labels
      4. Fallback to numbered channels
    """
    channel_names = []

    # 1. MIBI JSON
    try:
        first_desc = json.loads(tif.pages[0].description)
        if "channel.target" in first_desc:
            for page in tif.pages:
                desc = json.loads(page.description)
                channel_names.append(desc["channel.target"])
    except (json.JSONDecodeError, TypeError, KeyError):
        pass

    # 2. OME-XML
    if not channel_names:
        try:
            ome_xml = tif.ome_metadata
            if not ome_xml and tif.pages:
                desc_tag = tif.pages[0].tags.get("ImageDescription")
                ome_xml = desc_tag.value if desc_tag else tif.pages[0].description
            if ome_xml:
                root = ET.fromstring(ome_xml)
                pixels = None
                for elem in root.iter():
                    if _local_tag_name(elem.tag) == "Pixels":
                        pixels = elem
                        break
                if pixels is not None:
                    for elem in pixels:
                        if _local_tag_name(elem.tag) == "Channel":
                            channel_names.append(
                                elem.get("Name") or elem.get("ID")
                            )
        except Exception:
            pass

    # 3. ImageJ metadata labels
    if not channel_names:
        if hasattr(tif, "imagej_metadata") and tif.imagej_metadata:
            if "Labels" in tif.imagej_metadata:
                channel_names = list(tif.imagej_metadata["Labels"])

    # 4. Fallback
    if not channel_names:
        channel_names = [f"Channel_{i}" for i in range(n_channels)]

    return channel_names[:n_channels]


def get_resolution(tif):
    """
    Try to extract pixel size in microns from ImageJ or TIFF resolution tags.
    Returns (physical_size, unit) or (None, None).
    """
    # Try ImageJ metadata for pixel size
    if hasattr(tif, "imagej_metadata") and tif.imagej_metadata:
        ij = tif.imagej_metadata
        # ImageJ sometimes stores 'spacing' or uses resolution tags
        if "unit" in ij and ij["unit"] in ("um", "µm", "micron", "\\u00B5m"):
            # Resolution is stored in TIFF tags
            page = tif.pages[0]
            tags = page.tags
            res_unit = tags.get("ResolutionUnit")
            x_res = tags.get("XResolution")
            if x_res is not None:
                # XResolution is pixels-per-unit as a rational (num, denom)
                val = x_res.value
                if isinstance(val, tuple):
                    ppu = val[0] / val[1] if val[1] else 0
                else:
                    ppu = float(val)
                if ppu > 0:
                    return 1.0 / ppu, "µm"

    return None, None


def main():
    parser = argparse.ArgumentParser(
        description="Ensure a TIFF has valid OME-XML metadata."
    )
    parser.add_argument("input", type=Path, help="Input TIFF file")
    parser.add_argument("output", type=Path, help="Output OME-TIFF file")
    args = parser.parse_args()

    with tifffile.TiffFile(args.input) as tif:
        if has_valid_ome_xml(tif):
            # Already a valid OME-TIFF — just copy/symlink
            if args.input.resolve() != args.output.resolve():
                import shutil
                shutil.copy2(args.input, args.output)
            print(
                f"TIFF already has valid OME-XML, copied to {args.output}",
                file=sys.stderr,
            )
            return

        # Read image data
        data = tif.asarray()  # typically (C, Y, X) or (Y, X)

        n_pages = len(tif.pages)
        if data.ndim == 3 and n_pages == 1:
            # Single-page interleaved (Y, X, C) — transpose to (C, Y, X)
            if data.shape[2] < data.shape[0] and data.shape[2] < data.shape[1]:
                data = np.transpose(data, (2, 0, 1))

        n_channels = data.shape[0] if data.ndim == 3 else 1
        channel_names = get_channel_names(tif, n_channels)
        physical_size, _ = get_resolution(tif)

    print(
        f"Converting to OME-TIFF with channels: {channel_names}",
        file=sys.stderr,
    )

    # Build OME-TIFF metadata dict for tifffile
    metadata = {"Channel": {"Name": channel_names}}
    if physical_size is not None:
        metadata["PhysicalSizeX"] = physical_size
        metadata["PhysicalSizeY"] = physical_size
        metadata["PhysicalSizeXUnit"] = "µm"
        metadata["PhysicalSizeYUnit"] = "µm"

    tifffile.imwrite(
        args.output,
        data,
        ome=True,
        metadata=metadata,
        compression="deflate",
    )
    print(f"Wrote OME-TIFF to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
