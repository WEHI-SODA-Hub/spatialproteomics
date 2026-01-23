process CELLSAMSEGMENT {
    tag "$meta.id"
    label 'process_high'
    secret 'DEEPCELL_ACCESS_TOKEN'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(tiff), val(nuclear_channel), val(membrane_channels)
    val(compartment)

    output:
    tuple val(meta), path("*.tiff"), emit: segmentation_mask
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_channels = membrane_channels.first() != [] ? membrane_channels.first().split(":") : []
    def nuc_channel = nuclear_channel.first()
    
    """
    #!/usr/bin/env python3
    import numpy as np
    import tifffile
    from cellSAM.cellsam_pipeline import cellsam_pipeline
    import cellSAM
    import argparse
    import sys
    import os
    
    
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
        1. OME-XML metadata
        2. ImageJ metadata
        3. Fallback to numbered channels
        '''
        channel_names = []
        
        # Try OME-XML metadata
        if tif.is_ome:
            try:
                from xml.etree import ElementTree as ET
                ome_xml = tif.ome_metadata
                root = ET.fromstring(ome_xml)
                
                # Parse with OME namespace
                ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
                channels = root.findall('.//ome:Channel', ns)
                
                # Fallback: try without namespace
                if not channels:
                    channels = root.findall('.//Channel')
                
                channel_names = [ch.get('Name') or ch.get('ID') for ch in channels]
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
        parser.add_argument('--compartment', required=True, choices=['nuclear', 'whole-cell'],
                          help='Segmentation compartment')
        parser.add_argument('--nuclear-channel', required=True, help='Nuclear channel name')
        parser.add_argument('--membrane-channel', action='append', default=[],
                          help='Membrane channel name(s)')
        parser.add_argument('--bbox-threshold', type=float, default=${params.cellsam_bbox_threshold},
                          help='Bounding box threshold')
        parser.add_argument('--block-size', type=int, default=${params.cellsam_block_size},
                          help='Block size for tiling')
        parser.add_argument('--overlap', type=int, default=${params.cellsam_overlap},
                          help='Overlap between tiles')
        parser.add_argument('--iou-depth', type=int, default=${params.cellsam_iou_depth},
                          help='IOU depth parameter')
        parser.add_argument('--iou-threshold', type=float, default=${params.cellsam_iou_threshold},
                          help='IOU threshold for merging')
        parser.add_argument('--use-wsi', type=lambda x: x.lower() == 'true', 
                          default='${params.cellsam_use_wsi}', help='Use WSI mode')
        parser.add_argument('--gauge-cell-size', type=lambda x: x.lower() == 'true',
                          default='${params.cellsam_gauge_cell_size}', help='Gauge cell size')
        parser.add_argument('--low-contrast-enhancement', type=lambda x: x.lower() == 'true',
                          default='${params.cellsam_low_contrast_enhancement}',
                          help='Apply low contrast enhancement')
        parser.add_argument('--model-path', type=str, default=None,
                          help='Custom model path (optional)')
        
        return parser
    
    
    def run_segmentation(args):
        '''Run CellSAM segmentation pipeline.'''
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
        
        # Build argument list
        args_list = [
            "${tiff}",
            "--compartment", "${compartment}",
            "--nuclear-channel", "${nuc_channel}"
        ]
        
        # Add membrane channels
        membrane_channels_list = ${mem_channels.collect { "\"${it}\"" }}
        for mem_chan in membrane_channels_list:
            args_list.extend(["--membrane-channel", mem_chan])
        
        # Parse arguments
        parser = parse_arguments()
        args = parser.parse_args(args_list)
        
        # Run segmentation
        mask = run_segmentation(args)
        
        # Save segmentation mask
        output_path = "${prefix}_${compartment}.tiff"
        tifffile.imwrite(output_path, mask.astype(np.uint32))
        print(f"Segmentation mask saved to {output_path}")
        
        # Write versions
        with open("versions.yml", "w") as f:
            f.write('"${task.process}":\\n')
            f.write('    cellsam: 0.1.0\\n')
    
    
    if __name__ == "__main__":
        main()
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_${compartment}.tiff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsam: 0.1.0
    END_VERSIONS
    """
}
