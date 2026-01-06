# Testing CellSAM Integration

This document describes the tests for the CellSAM segmentation integration.

## Test Structure

The CellSAM tests follow the same structure as the existing Mesmer and Cellpose tests:

### Module Tests

**Location**: `modules/local/cellsamsegment/tests/`

Tests the `CELLSAMSEGMENT` process in isolation:
- Nuclear segmentation
- Whole-cell segmentation  
- Stub test

**Run tests:**
```bash
nf-test test modules/local/cellsamsegment/tests/main.nf.test
```

### Subworkflow Tests

**Location**: `subworkflows/local/cellsam_segment/tests/`

Tests the complete `CELLSAM_SEGMENT` subworkflow including:
- CellSAM segmentation (nuclear + whole-cell)
- Cell measurements
- Report generation

**Run tests:**
```bash
nf-test test subworkflows/local/cellsam_segment/tests/main.nf.test
```

**Location**: `subworkflows/local/cellsam_segment_wbacksub/tests/`

Tests the `CELLSAM_SEGMENT_WBACKSUB` subworkflow with background subtraction.

**Run tests:**
```bash
nf-test test subworkflows/local/cellsam_segment_wbacksub/tests/main.nf.test
```

## Integration Tests

### Test Profile

**Profile**: `test_cellsam`

**Configuration**: `conf/test_cellsam.config`

**Samplesheet**: `tests/cellsam_samplesheet.csv`

**Test Data**: Uses existing test data in `tests/data/comet/synthetic_multichannel.tif`

**Run integration test:**
```bash
# With conda
nextflow run . -profile test_cellsam,conda --outdir results_cellsam_test

# With docker
nextflow run . -profile test_cellsam,docker --outdir results_cellsam_test

# With singularity
nextflow run . -profile test_cellsam,singularity --outdir results_cellsam_test
```

### Test Parameters

The test profile uses smaller values for faster testing:
- `cellsam_block_size = 256` (default: 400)
- `cellsam_overlap = 32` (default: 56)
- `cellsam_bbox_threshold = 0.4` (default: 0.4)
- `cellsam_use_wsi = true` (tests tiling functionality)

## Test Samplesheet

The test samplesheet (`tests/cellsam_samplesheet.csv`) includes:
```csv
sample,run_backsub,run_cellpose,run_mesmer,run_cellsam,tiff,nuclear_channel,membrane_channels
test,true,false,false,true,tests/data/comet/synthetic_multichannel.tif,DAPI,CD45
```

This tests:
- ✅ Background subtraction enabled
- ✅ CellSAM segmentation enabled
- ✅ Nuclear and membrane channel specification
- ✅ COMET platform data format

## Running All Tests

### Run all nf-tests
```bash
# Run all module and subworkflow tests
nf-test test
```

### Run specific test tags
```bash
# Run only CellSAM tests
nf-test test --tag cellsamsegment

# Run only subworkflow tests
nf-test test --tag subworkflows
```

### Run integration test
```bash
# Quick test with conda
nextflow run . -profile test_cellsam,conda --outdir test_output

# Test with specific parameters
nextflow run . -profile test_cellsam,docker \\
    --outdir test_output \\
    --cellsam_bbox_threshold 0.3 \\
    --cellsam_block_size 200
```

## Expected Outputs

After a successful test run, you should see:

### Segmentation Masks
- `test_nuclear.tiff` - Nuclear segmentation mask
- `test_whole-cell.tiff` - Whole-cell segmentation mask

### Measurements
- `test_annotations.parquet` - Cell measurements with spatial coordinates

### Reports (if `generate_report = true`)
- `test_report.html` - Segmentation QC report

### Background Subtraction
- `test_bgsub.tif` - Background-subtracted image (if `run_backsub = true`)

## Continuous Integration

Tests can be integrated into CI/CD pipelines:

### GitHub Actions Example
```yaml
- name: Run CellSAM tests
  run: |
    nf-test test --tag cellsamsegment
    nextflow run . -profile test_cellsam,docker --outdir test_results
```

## Test Data

The tests use existing synthetic test data located in:
- `tests/data/comet/synthetic_multichannel.tif` (COMET format)
- `tests/data/mesmer/test_data.tiff` (Mesmer format)

To add custom test data:
1. Place TIFF file in `tests/data/`
2. Update samplesheet with correct path and channels
3. Run test with new data

## Troubleshooting Tests

### GPU Tests
CellSAM will use GPU if available but falls back to CPU. For CI/CD environments without GPUs:
- Tests will run slower on CPU
- Consider using smaller test images or reduced `block_size`

### Memory Issues
If tests fail with OOM errors:
- Reduce `cellsam_block_size` in test config
- Increase memory allocation in test profile
- Use smaller test images

### Missing Dependencies
Ensure the CellSAM conda environment includes all dependencies:
- PyTorch
- CellSAM package
- Segment Anything Model (SAM)

Check `modules/local/cellsamsegment/environment.yml`

## Snapshot Testing

nf-test uses snapshots to verify outputs remain consistent. After making changes:

```bash
# Update snapshots
nf-test test --update-snapshot

# Review snapshot changes
git diff modules/local/cellsamsegment/tests/main.nf.test.snap
```

## Adding New Tests

To add a new test case:

1. **Module test**: Add to `modules/local/cellsamsegment/tests/main.nf.test`
2. **Subworkflow test**: Add to `subworkflows/local/cellsam_segment/tests/main.nf.test`
3. **Integration test**: Create new config in `conf/` and samplesheet in `tests/`

Example new integration test:
```nextflow-config
// conf/test_cellsam_no_backsub.config
params {
    input = 'tests/cellsam_no_backsub_samplesheet.csv'
    // ... other params
}
```

Then add to profiles in `nextflow.config`:
```nextflow
test_cellsam_no_backsub { includeConfig 'conf/test_cellsam_no_backsub.config' }
```
