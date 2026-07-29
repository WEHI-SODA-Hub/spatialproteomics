# WEHI-SODA-Hub/sp_segment: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Changed

- **Cellpose model weights are fetched once per run** by a new `CELLPOSEMODEL`
  process and staged as an input to every segmentation task. Each patch task
  previously downloaded the ~1.2 GB `cpsam` checkpoint from cellpose.org itself
  — about 18 identical requests for a single sample — which produced HTTP 429
  rate limiting in CI and truncated downloads locally. The staged copy is reused
  by `-resume`, and every patch task in a run now provably reads the same file.

- **A missing Cellpose model path is now a startup error.** Cellpose does not
  fail on a `--pretrained-model` path that does not exist; it falls back to its
  built-in weights, so a typo produced a complete, plausible, silently wrong
  run. `cellpose_pretrained_model` and `cellpose_models_dir` are both checked in
  `PIPELINE_INITIALISATION`, before any compute.

- **The `backgroundsubtract` snapshot no longer records the output TIFF's md5.**
  The upstream backsub container stamps each OME-TIFF with a fresh UUID, so that
  md5 changed on every run and the snapshot could never pass twice. The test now
  asserts the file is produced by name, matching the treatment `kronos_input`
  already gets.

- **KRONOS1 has been replaced by KRONOS2.** `bin/kronos_embeddings.py` and the
  `KRONOSEMBEDDINGS` module are removed; embeddings are now 768-d rather than 384-d.
  Outputs from the two versions are not comparable.

- **KRONOS now runs once at the top level**, not three times inside the Mesmer,
  Cellpose and CellSAM subworkflows. A cohort-wide operator inside a per-segmenter
  subworkflow only sees that segmenter's samples, which made any pooled statistic
  unsound.

- **Cell patches are isolated to the target cell.** Pixels outside the cell polygon
  are zeroed, so an embedding describes that cell rather than the ~18 um
  neighbourhood a 64 px box spans. Disable with `--kronos_isolate_cell false`.

- **The intensity divisor is derived from the image dtype** (uint8=255, uint16=65535,
  float=400) instead of a hardcoded 65535, matching KRONOS2's own scaling factor. A
  float32 image was previously divided by 65535 rather than 400.

- **The nuclear channel is passed to the model** as `preferred_dapi`. This is
  KRONOS2's only marker-alias mechanism, so slides stained with DRAQ5 or Hoechst
  previously received default normalisation statistics.

- `kronos_batch_size` now defaults to `16` (was `32`), which reproduces the published
  values. KRONOS2 runs fp32 and cuBLAS picks batch-dependent kernels.

- **The per-cell embeddings CSV is no longer written.** It duplicated the vectors
  already merged into the GeoJSON, which is now the only embedding artefact.

### Added

- `cellpose_models_dir` parameter, pointing at a pre-staged Cellpose model
  cache so sites with no outbound network on compute nodes skip the download.

- `bin/kronos2_common.py`, a shared scaffold for per-cell foundation-model
  extraction, so a future encoder does not become a second copy of the script.
- A Wave-built container pinned by digest
  (`community.wave.seqera.io/library/kronos2embeddings:70590359503eec08`),
  carrying torch 2.6.0+cu124. Verified to reproduce a standalone run's
  embeddings to cosine 1.000000.
- `kronos_nuclear_marker`, `kronos_isolate_cell` and `kronos_allow_novel_defaults`
  parameters; `test_mesmer_kronos2` profile.
- Unmatched markers now fail the run by default. KRONOS2 matches names exactly with
  no alias or fuzzy step, so a naming variant such as `Cytokeritin` would otherwise be
  normalised with default statistics silently. Use `--kronos_marker_mapping` to
  resolve them, or `--kronos_allow_novel_defaults` to accept the fallback.

- Per-cell footprints are filled with Pillow rather than rasterio, dropping the
  GDAL stack (which pulled a conda libtiff/libjpeg symbol clash into the
  container). The two rasterisers differ slightly at the cell boundary:
  embeddings shift by cosine 0.99974, roughly 15x less than the cell-isolation
  change above, and less than rasterio's own two `all_touched` modes differ
  from each other.

### Removed

- `kronos_marker_metadata` -- KRONOS2 carries its own 288-marker vocabulary.
- `kronos_config_path` -- KRONOS2 has no equivalent config file.
- `kronos_distance_threshold` -- the centroid-matching fallback it configured was
  never implemented; cells now join by GeoJSON feature order.

### Fixed

- Stub runs of the Cellpose path failed with `For input string: ""`. The
  `SOPA_PATCHIFYIMAGE` stub wrote an empty patch-count file that the caller parses
  with `.toInteger()`.
- `KRONOSEMBEDDINGS` passed `--nv` unconditionally, which is invalid under Docker.
  The replacement selects `--nv` or `--gpus all` by container engine.

## v0.4.0 - 2026-04-20

### Changed

- **cellmeasurement: Python replaces Groovy as the only implementation.** The original Groovy app
  (`https://github.com/WEHI-SODA-Hub/cellmeasurement`) is no longer used. All measurement logic
  is now in https://github.com/WEHI-SODA-Hub/cellmeasurement-py.

- **cellmeasurement: tiled measurement processing has been re-implemented.** User can specify tile
  size and overlap, and the app will process measurements for tiles in parallel.

- **cellmeasurement: erosion measurement column names have changed.** The Groovy app produced
  columns named `{Channel}: {Compartment}: Eroded_{N}px: Mean/Median` where N was a fixed pixel
  depth. The Python implementation uses 5 _equal-area_ bins and names them
  `{Compartment}: ErosionBin_{N}: Mean/Median`. If you have existing Groovy-generated GeoJSON
  that you compare with new Python output, expect different column names.

- **cellmeasurement: `--erosion-steps` no longer accepts pixel-depth values.** Passing numeric
  values such as `--erosion-steps=4,7,11,14,18` (Groovy API) is not supported. Use the boolean
  flag `--erosion-steps` to enable the 5 equal-area bins.

- pipeline can be run with whole-cell segmentation only (no nucleus segmentation)

### Added

- **KRONOS embedding output**: new `kronosembeddings` module extracts per-cell embeddings using
  the KRONOS foundation model. Produces `*_kronos_embeddings.csv` (cell IDs + 384 embedding
  dimensions), `*_marker_report.txt` (channel-to-marker match summary), and optionally
  `*_kronos_merged.geojson` (embeddings merged into the cellmeasurement GeoJSON).
  Controlled by `enable_kronos`, `kronos_model_path`, `kronos_marker_metadata`, and related
  `kronos_*` parameters. Disabled by default (`enable_kronos = false`).

- **CellSAM segmentation**: new `cellsam_segment` module adds support for the CellSAM foundation
  model as a third segmentation option alongside Mesmer and Cellpose. Supports tiled (WSI) and
  non-tiled inference via `cellsam_use_wsi`. Controlled by `cellsam_bbox_threshold`,
  `cellsam_block_size`, `cellsam_overlap`, `cellsam_iou_threshold`, and `cellsam_use_wsi`
  parameters.

- **Mask smoothing (`smooth_masks`)**: new optional `smoothmasks` module reduces polygon
  complexity. Two methods: `morphological` (disk close+open, default) and
  `shapely` (Douglas-Peucker polygon simplification). Enabled via `smooth_masks
= true`; tunable with `smooth_method` and `smooth_kernel_size`.

- `dist_threshold` pipeline parameter: controls the maximum centroid distance (pixels) for
  matching a nucleus to a whole-cell ROI in cellmeasurement (default: 10.0).
- `downsample_factor` pipeline parameter: integer downsample factor applied to image and masks
  before cellmeasurement to speed up processing on very large images (default: 1.0 = disabled).

### Removed

- `gradle_cache_dir` parameter and all references removed (leftover from the Groovy era).

## v0.3.0 - 2025-11-13

- Rename to sp_segment
- Fix issue with OME-XML processing
- Add background image preview in report
- Support for embedded report
- Add percentile calculation
- Update workflow diagram

## v0.2.0 - 2025-09-12

Add segmentation report.

## v0.1.0 - 2025-08-29

Initial release of WEHI-SODA-Hub/sp_segment, created with the [nf-core](https://nf-co.re/) template.
