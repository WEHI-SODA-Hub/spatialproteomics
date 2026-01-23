# CellSAM Integration

This document describes the integration of CellSAM as a segmentation option in the sp_segment pipeline.

## Overview

CellSAM is a foundation model for cell segmentation that can handle various imaging modalities including multiplexed imaging data from COMET and MIBI platforms. It has been integrated alongside Mesmer and Cellpose as a segmentation option.

## Usage

### Samplesheet Configuration

To use CellSAM segmentation, add the `run_cellsam` column to your samplesheet:

```csv
sample,run_backsub,run_cellpose,run_mesmer,run_cellsam,tiff,nuclear_channel,membrane_channels
sample1,false,false,false,true,/path/to/image.tif,DAPI,CD45:PanCK
sample2,true,false,false,true,/path/to/image2.tif,Hoechst,Ecad:Vim
```

### Parameters

CellSAM-specific parameters can be configured in your nextflow config or via command line:

- `--cellsam_bbox_threshold` (default: 0.4): Confidence threshold for cell detection. Lower values increase recall but may reduce precision.
- `--cellsam_block_size` (default: 400): Size of tiles for large image processing. Smaller values (256-512) work better for dense images.
- `--cellsam_overlap` (default: 56): Tile overlap region for label merging. Should encompass typical cell size.
- `--cellsam_iou_depth` (default: 56): Depth for IOU-based label merging across tiles.
- `--cellsam_iou_threshold` (default: 0.5): IOU threshold for merging overlapping labels.
- `--cellsam_use_wsi` (default: true): Enable tiling for large images (recommended for >3000 cells).
- `--cellsam_gauge_cell_size` (default: false): Perform preliminary segmentation to estimate cell sizes.
- `--cellsam_low_contrast_enhancement` (default: false): Apply contrast enhancement preprocessing.
- `--cellsam_model_path` (default: null): Path to custom model weights (uses latest generalist model if null).

### Example Command

```bash
nextflow run main.nf \\
    --input samplesheet.csv \\
    --outdir results \\
    --cellsam_bbox_threshold 0.4 \\
    --cellsam_block_size 400 \\
    -profile conda
```

## Channel Format

CellSAM expects channels in the format: `(blank, nuclear, membrane)`

The module automatically:
1. Extracts the specified nuclear channel from the input TIFF
2. For whole-cell segmentation, extracts the membrane channel(s)
3. Formats them into the 3-channel array expected by CellSAM

## Background Subtraction

CellSAM can be combined with background subtraction by setting both `run_backsub=true` and `run_cellsam=true` in the samplesheet. The workflow will:
1. Perform background subtraction
2. Run CellSAM segmentation on the cleaned image
3. Generate cell measurements and optional report

## Output

CellSAM produces:
- Nuclear segmentation mask: `{sample}_nuclear.tiff`
- Whole-cell segmentation mask: `{sample}_whole-cell.tiff`
- Cell measurements: `{sample}_annotations.parquet`
- Optional segmentation report: `{sample}_report.html` (if `--generate_report true`)

## Requirements

The CellSAM module requires:
- Python 3.10+
- PyTorch (with GPU support recommended)
- CellSAM package and dependencies
- Segment Anything Model (SAM)

These are automatically installed via the conda environment specified in `modules/local/cellsamsegment/environment.yml`.

## Notes

- GPU acceleration is automatically used if available
- For very large images or many samples, consider using a computing cluster with GPU resources
- The DEEPCELL_ACCESS_TOKEN secret is inherited from the base configuration but not required for CellSAM
- If DEEPCELL_ACCESS_TOKEN secret is present the pipeline will install the latest model weights(1.2), otherwise will default to base model weights
- CellSAM works well across different imaging modalities without fine-tuning due to its foundation model architecture

## Troubleshooting

- **Out of memory errors**: Reduce `cellsam_block_size` or ensure GPU has sufficient memory
- **Poor segmentation**: Adjust `cellsam_bbox_threshold` (lower for out-of-distribution images)
- **Missing cells**: Increase `cellsam_overlap` to ensure cells at tile boundaries are captured
- **False positives**: Increase `cellsam_bbox_threshold` for stricter filtering
