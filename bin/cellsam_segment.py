#!/usr/bin/env python3
import numpy as np
import tifffile
from cellSAM.cellsam_pipeline import cellsam_pipeline
import cellSAM
import argparse
import sys
import os


def _local_tag_name(tag):
    '''Return XML tag name without namespace.'''
    return tag.split('}', 1)[-1]


def parse_ome_channel_names(ome_xml):
    '''Parse channel names from OME-XML (namespace-agnostic).'''
    from xml.etree import ElementTree as ET

    root = ET.fromstring(ome_xml)
    pixels = None

    for elem in root.iter():
        if _local_tag_name(elem.tag) == 'Pixels':
            pixels = elem
            break

    if pixels is None:
        return []

    channel_names = []
    for elem in pixels:
        if _local_tag_name(elem.tag) == 'Channel':
            channel_names.append(elem.get('Name') or elem.get('ID'))

    return channel_names


def download_model_weights():
    '''Download latest CellSAM model weights (v1.2) from users.deepcell.org.'''
    if 'DEEPCELL_ACCESS_TOKEN' in os.environ:
        print("Downloading/verifying latest CellSAM model weights (v1.2)...")
        cellSAM.get_model()
    else:
        print("Warning: DEEPCELL_ACCESS_TOKEN not set. Using default model weights.")


def get_channel_names(tif, n_channels):
    '''
    Extract channel names from TIFF metadata.

    Tries multiple strategies in order:
    1. MIBI JSON metadata (per-page JSON with channel.target)
    2. OME-XML metadata (including OPAL OME/QPTIFF exports)
    3. ImageJ metadata
    4. Fallback to numbered channels
    '''
    import json
    channel_names = []

    # Try MIBI JSON metadata (each page has a JSON description with channel.target)
    try:
        first_desc = json.loads(tif.pages[0].description)
        if 'channel.target' in first_desc:
            channel_names = []
            for page in tif.pages:
                desc = json.loads(page.description)
                channel_names.append(desc['channel.target'])
    except (json.JSONDecodeError, TypeError, KeyError):
        pass

    # Try OME-XML metadata
    if not channel_names:
        try:
            ome_xml = tif.ome_metadata
            if not ome_xml and tif.pages:
                desc_tag = tif.pages[0].tags.get("ImageDescription")
                ome_xml = desc_tag.value if desc_tag else tif.pages[0].description

            if ome_xml:
                channel_names = parse_ome_channel_names(ome_xml)
        except Exception:
            pass

    # Try ImageJ metadata
    if not channel_names and hasattr(tif, 'imagej_metadata') and tif.imagej_metadata:
        if 'Labels' in tif.imagej_metadata:
            channel_names = tif.imagej_metadata['Labels']

    # Fallback to numbered channels
    if not channel_names:
        channel_names = [f"Channel_{i}" for i in range(n_channels)]

    return channel_names


def find_channel_index(channel_name, channel_names):
    '''
    Find the index of a channel by name.

    First tries exact match, then case-insensitive match.
    Returns None if not found.
    '''
    # Exact match
    if channel_name in channel_names:
        return channel_names.index(channel_name)

    # Case-insensitive match
    for i, name in enumerate(channel_names):
        if name.lower() == channel_name.lower():
            return i

    return None


def create_cellsam_input(img, nuc_idx, mem_idx=None):
    '''
    Create a 3-channel array for CellSAM input.

    Format: [blank, nuclear, membrane]
    If membrane channel is not provided, it remains blank (zeros).
    '''
    if img.ndim == 2:
        h, w = img.shape
        cellsam_img = np.zeros((h, w, 3), dtype=img.dtype)
        cellsam_img[..., 1] = img
    else:
        h, w = img.shape[1], img.shape[2]
        cellsam_img = np.zeros((h, w, 3), dtype=img.dtype)
        cellsam_img[..., 1] = img[nuc_idx, ...]

        if mem_idx is not None:
            cellsam_img[..., 2] = img[mem_idx, ...]

    return cellsam_img


def remove_border_cells(mask):
    '''
    Remove cells that touch the image borders.

    Args:
        mask: 2D numpy array with cell labels

    Returns:
        Filtered mask with border-touching cells removed (set to 0)
    '''
    h, w = mask.shape
    border_cells = set()

    # Find cells touching top and bottom edges
    border_cells.update(mask[0, :])
    border_cells.update(mask[h-1, :])

    # Find cells touching left and right edges
    border_cells.update(mask[:, 0])
    border_cells.update(mask[:, w-1])

    # Remove background (0) from border cells set
    border_cells.discard(0)

    # Create filtered mask
    filtered_mask = mask.copy()
    for cell_id in border_cells:
        filtered_mask[mask == cell_id] = 0

    n_removed = len(border_cells)
    n_total = len(np.unique(mask)) - 1  # -1 for background
    print(f"Removed {n_removed} border-touching cells (kept {n_total - n_removed}/{n_total})")

    return filtered_mask


def remove_small_cells(mask, min_area):
    '''
    Remove cells smaller than a minimum pixel area.

    Args:
        mask: 2D numpy array with cell labels
        min_area: Minimum number of pixels for a cell to be retained

    Returns:
        Filtered mask with small cells removed (set to 0)
    '''
    if min_area <= 0:
        return mask

    labels, counts = np.unique(mask, return_counts=True)
    small = labels[(labels != 0) & (counts < min_area)]

    if small.size == 0:
        print(f"No cells below min-area {min_area}px")
        return mask

    filtered = mask.copy()
    filtered[np.isin(filtered, small)] = 0

    n_total = int((labels != 0).sum())
    print(f"Removed {len(small)} cells below {min_area}px (kept {n_total - len(small)}/{n_total})")
    return filtered


def extract_channels(tiff_path, nuclear_channel, membrane_channels, compartment):
    '''
    Extract and format channels from a multi-channel TIFF for CellSAM.

    Args:
        tiff_path: Path to the input TIFF file
        nuclear_channel: Name of the nuclear channel
        membrane_channels: List of membrane channel names (for whole-cell segmentation)
        compartment: 'nuclear' or 'whole-cell'

    Returns:
        3-channel numpy array formatted for CellSAM [blank, nuclear, membrane]
    '''
    # Load the TIFF image
    img = tifffile.imread(tiff_path)

    # Get channel names from metadata
    with tifffile.TiffFile(tiff_path) as tif:
        # Check if single-page interleaved TIFF (Y, X, C) format
        if img.ndim == 3 and tif.is_ome and len(tif.pages) == 1:
            # Try to get channel count from OME metadata
            try:
                from xml.etree import ElementTree as ET
                root = ET.fromstring(tif.ome_metadata)
                ns = {'ome': root.tag.split('}')[0].strip('{')}
                pixels = root.find('.//ome:Pixels', ns)
                n_channels = int(pixels.get('SizeC', img.shape[0]))
                # If last dimension matches channel count, it's interleaved (Y, X, C)
                if img.shape[2] == n_channels:
                    img = np.transpose(img, (2, 0, 1))  # Convert to (C, Y, X)
            except:
                pass

        n_channels = img.shape[0] if img.ndim > 2 else 1
        channel_names = get_channel_names(tif, n_channels)

    print(f"Available channels: {channel_names}")

    # Find nuclear channel
    nuc_idx = find_channel_index(nuclear_channel, channel_names)
    if nuc_idx is None:
        raise ValueError(f"Nuclear channel '{nuclear_channel}' not found in {channel_names}")

    # Find membrane channel (for whole-cell segmentation)
    mem_idx = None
    if compartment == "whole-cell" and membrane_channels:
        for mem_chan in membrane_channels:
            mem_idx = find_channel_index(mem_chan, channel_names)
            if mem_idx is not None:
                break

        if mem_idx is None:
            print(f"Warning: Membrane channel not found. Using mean of non-nuclear channels.")
            # Fallback: use mean of all non-nuclear channels
            if img.ndim > 2:
                all_channels = np.arange(img.shape[0])
                other_channels = all_channels[all_channels != nuc_idx]
                if len(other_channels) > 0:
                    # Create pseudo-membrane channel
                    mem_img = img[other_channels, ...].mean(axis=0)
                    # Temporarily add to image stack for processing
                    img = np.vstack([img, mem_img[np.newaxis, ...]])
                    mem_idx = img.shape[0] - 1

    return create_cellsam_input(img, nuc_idx, mem_idx)


def parse_arguments():
    '''Parse command-line arguments for CellSAM segmentation.'''
    parser = argparse.ArgumentParser(description='CellSAM segmentation')
    parser.add_argument('tiff', help='Input TIFF file path')
    parser.add_argument('--output', required=True, help='Output mask file path')
    parser.add_argument('--compartment', required=True, choices=['nuclear', 'whole-cell'],
                      help='Segmentation compartment')
    parser.add_argument('--nuclear-channel', required=True, help='Nuclear channel name')
    parser.add_argument('--membrane-channel', action='append', default=[],
                      help='Membrane channel name(s)')
    parser.add_argument('--bbox-threshold', type=float, default=0.4,
                      help='Bounding box threshold')
    parser.add_argument('--block-size', type=int, default=712,
                      help='Block size for tiling')
    parser.add_argument('--overlap', type=int, default=120,
                      help='Overlap between tiles')
    parser.add_argument('--iou-depth', type=int, default=120,
                      help='IOU depth parameter')
    parser.add_argument('--iou-threshold', type=float, default=0.5,
                      help='IOU threshold for merging')
    parser.add_argument('--use-wsi', action='store_true',
                      help='Use WSI mode')
    parser.add_argument('--gauge-cell-size', action='store_true',
                      help='Gauge cell size')
    parser.add_argument('--low-contrast-enhancement', action='store_true',
                      help='Apply low contrast enhancement')
    parser.add_argument('--model-path', type=str, default=None,
                      help='Custom model path (optional)')
    parser.add_argument('--remove-border-cells', action='store_true',
                      help='Remove cells touching image borders')
    parser.add_argument('--min-area', type=int, default=0,
                      help='Minimum cell area in pixels. Cells smaller than this are removed (default: 0 = disabled)')

    return parser.parse_args()


def run_segmentation(args):
    '''Run CellSAM segmentation pipeline.'''
    import torch
    print(f"PyTorch version: {torch.__version__}", file=sys.stderr)
    print(f"CUDA available: {torch.cuda.is_available()}", file=sys.stderr)
    print(f"CUDA version: {torch.version.cuda}", file=sys.stderr)
    if torch.cuda.is_available():
        print(f"GPU available: {torch.cuda.get_device_name(0)}", file=sys.stderr)
    else:
        print("WARNING: No GPU detected - CellSAM will run on CPU and may be very slow.", file=sys.stderr)

    # Extract and format channels
    img = extract_channels(
        args.tiff,
        args.nuclear_channel,
        args.membrane_channel,
        args.compartment
    )

    # Run CellSAM segmentation
    mask = cellsam_pipeline(
        img,
        bbox_threshold=args.bbox_threshold,
        block_size=args.block_size,
        overlap=args.overlap,
        iou_depth=args.iou_depth,
        iou_threshold=args.iou_threshold,
        use_wsi=args.use_wsi,
        gauge_cell_size=args.gauge_cell_size,
        low_contrast_enhancement=args.low_contrast_enhancement,
        model_path=args.model_path
    )

    return mask


def main():
    '''Main execution function.'''
    # Download/verify model weights
    download_model_weights()

    # Parse arguments
    args = parse_arguments()

    # Run segmentation
    mask = run_segmentation(args)

    # Always remove border-touching cells
    mask = remove_border_cells(mask)

    # Remove small cells if threshold set
    if args.min_area > 0:
        mask = remove_small_cells(mask, args.min_area)

    # Save segmentation mask
    tifffile.imwrite(args.output, mask.astype(np.uint32))
    print(f"Segmentation mask saved to {args.output}")


if __name__ == "__main__":
    main()
