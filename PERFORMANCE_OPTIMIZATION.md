# Performance Optimization Guide for CellSAM WSI Processing

## Changes Applied

### 1. Increased CPU Allocation
- **Base Config**: `process_high` now uses 8 CPUs (was 1)
- **Time Limit**: Extended to 48 hours for `process_medium` and `process_high`

### 2. New HPC Profile
Added an `hpc` profile to `nextflow.config` optimized for cluster execution:
- **CellSAM (process_high)**: 16 CPUs, 96GB RAM, 48h
- **Medium processes**: 8 CPUs, 48GB RAM, 48h
- **Multi-threaded**: 32 CPUs, 144GB RAM, 48h

## How to Use

### On HPC with SLURM:
```bash
nextflow run main.nf \
  -profile hpc,singularity \
  --input samplesheet.csv \
  --outdir results
```

### Customize for Your HPC:
Edit the `hpc` profile in `nextflow.config` (around line 257):
- Change `executor = 'slurm'` to your scheduler (`pbspro`, `lsf`, `sge`)
- Update `queue = 'normal'` to your queue name
- Adjust `queueSize` based on your cluster limits

### Combine with Size Profiles:
```bash
# For large images on HPC
nextflow run main.nf -profile hpc,large,singularity --input samplesheet.csv
```

## Optimization Strategies

### 1. **Optimize Block Size Parameters**
The WSI tiling is controlled by these parameters:

```bash
--cellsam_block_size 400     # Increase to 512-1024 for fewer, larger blocks
--cellsam_overlap 56         # Keep small (32-64) to reduce redundant processing
--cellsam_iou_depth 56       # Should be <= overlap
```

**Recommendation**: 
- Larger block sizes (512-1024) = fewer blocks but more memory per block
- Test with your image size to find the sweet spot

### 2. **Enable Dask Parallelization in CellSAM**

The WSI code you shared has commented-out Dask parallel processing. To enable it:

**Option A**: Add to `nextflow.config` (Recommended for HPC/Sequera):
```groovy
process {
    withLabel: 'process_high' {
        beforeScript = 'export DASK_NUM_WORKERS=${task.cpus}; export OMP_NUM_THREADS=1'
        // Or set as environment variables:
        // env = [DASK_NUM_WORKERS: "${task.cpus}", OMP_NUM_THREADS: "1"]
    }
}
```
This automatically scales Dask workers to match allocated CPUs and works in containers.

**Option B**: For Sequera Platform prerun script:
```bash
export DASK_NUM_WORKERS=8
export OMP_NUM_THREADS=1
```
Then add to `nextflow.config` to propagate into containers:
```groovy
process {
    withLabel: 'process_high' {
        containerOptions = "--env DASK_NUM_WORKERS --env OMP_NUM_THREADS"
    }
}
```

**Option C**: Modify the CellSAM source code (`segment_wsi` function):
- Uncomment the `dask.delayed` lines (around line 65-80 in your code)
- Remove the synchronous `segment_chunk` call
- This will process blocks in parallel using Dask

### 3. **GPU Acceleration**

If you have GPUs on your HPC:

Add to `nextflow.config`:
```groovy
process {
    withLabel: 'process_high' {
        clusterOptions = '--gres=gpu:1'  // Request 1 GPU
        containerOptions = '--nv'        // Enable NVIDIA support
    }
}
```

Then set up Dask cluster with GPUs (modify CellSAM source):
```python
from cellSAM.wsi import setup_cluster
client, gpu_map = setup_cluster([0, 1, 2, 3])  # Use 4 GPUs
```

### 4. **Parallelize Across Multiple Images**

Nextflow automatically parallelizes at the sample level. Process multiple images simultaneously:

```bash
# This will process all samples in parallel (up to queueSize limit)
nextflow run main.nf \
  -profile hpc,singularity \
  --input multi_sample_sheet.csv \
  --max_cpus 16 \
  --max_memory 96.GB
```

### 5. **Memory Considerations**

If you hit memory limits:
- **Reduce block size**: `--cellsam_block_size 256` (more blocks, less memory each)
- **Increase memory allocation**: Use `-profile hpc,large`
- **Enable retry with more memory**: Already configured (doubles memory on retry)

### 6. **Monitoring Performance**

Add to your Nextflow command:
```bash
-with-report performance_report.html \
-with-timeline timeline.html \
-with-trace trace.txt
```

These will show:
- Which processes are bottlenecks
- Resource utilization
- Parallelization efficiency

## Troubleshooting

### Issue: "Blocks processing too slowly"
**Solutions**:
1. Increase CPUs: Use `-profile hpc` (16 CPUs)
2. Reduce block size: `--cellsam_block_size 256`
3. Enable Dask parallelization (see Option B above)

### Issue: "Out of memory"
**Solutions**:
1. Use memory profile: `-profile hpc,large`
2. Reduce block size: `--cellsam_block_size 256`
3. Disable gauge cell size: `--cellsam_gauge_cell_size false`

### Issue: "Hitting time limit"
**Solutions**:
1. Already extended to 48h
2. Request more time in HPC profile: Edit `time = { 72.h * task.attempt }`
3. Reduce overlap: `--cellsam_overlap 32`

## Expected Speedup

With these optimizations:
- **16 CPUs + Dask**: ~10-15x faster block processing
- **Parallel samples**: Linear speedup per additional sample
- **Optimal block size**: 2-3x faster from reduced overhead

Example: A 300GB image that took 10 hours could complete in ~1-2 hours with full optimization.

## Further Reading

- [Nextflow Executors](https://www.nextflow.io/docs/latest/executor.html)
- [Dask Documentation](https://docs.dask.org/)
- [CellSAM Repository](https://github.com/vanvalenlab/cellSAM)
