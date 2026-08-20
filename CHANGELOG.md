# WEHI-SODA-Hub/sp_segment: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Changed

- **Cellpose preprocessing now matches Cellpose's own contract: sopa's gaussian
  filter and CLAHE are both disabled.** Cellpose normalises by rescaling to the 1st
  and 99th percentiles and does nothing else — it contains no CLAHE code, and
  Cellpose-SAM was trained on exactly that normalisation. sopa layered two extra
  stages on top of every patch: `gaussian_filter(sigma=1)` and
  `equalize_adapthist(clip_limit=0.2)`. New params `cellpose_clip_limit` (`0`),
  `cellpose_gaussian_sigma` (`0`), `cellpose_clahe_kernel_size` (`null`) and
  `cellpose_tile_norm_blocksize` (`0`, Cellpose's own default) control them; restore
  the old behaviour with `--cellpose_clip_limit 0.2 --cellpose_gaussian_sigma 1`.

  Both stages were computed per patch, and CLAHE's kernel is `patch.shape // 8`, so
  the contrast-enhancement scale tracked `patch_width_pixel` — at COMET's 0.28 µm/px,
  52 µm at 1500 px, 105 µm at 3000 and 210 µm at 6000 for identical tissue. Measured
  on a real COMET image, segmenting the same region as part of a 1500 px vs a 3000 px
  patch and matching cells by IoU, the fraction of cell boundaries surviving a
  patch-size change (IoU ≥ 0.90) went from 38.6% to 81.9% for `cpsam_v2`, 43.2% to
  83.9% for `cpsam` and 36.8% to 86.6% for `cpdino`.

  **This changes results and they should not be pooled with earlier runs.** Median
  cell area falls 23–43% for every model, and cell counts move by +8–35%
  (`cpsam_v2`), −3–10% (`cpsam`) and +28–54% (`cpdino`) — CLAHE was suppressing
  detections in the newer models and inflating them in `cpsam`. `cellpose_diameter`
  needs no adjustment: Cellpose-SAM trained across 0.25×–4× scale, so the smaller
  cells sit well inside its range.

  Verified end to end on the `test` profile (30 tasks, nine-patch path), and run on a
  COMET brain-panel ROI with `cpsam_v2`, `cpdino` and cellSAM then inspected in
  QuPath — all three sensible and closely comparable (1,742 / 1,751 / 1,648 cells),
  with cellSAM drawing larger boundaries around membrane signal.

  These numbers measure reproducibility, not accuracy; no annotations were scored,
  and the normalisation measurements are DAPI-only. See
  [docs/usage.md](docs/usage.md) for the full reasoning.

- **The Singularity/Apptainer pull timeout is raised from Nextflow's 20 minute
  default to 1 hour.** That default is sized for the small images most nf-core
  pipelines pull; ours are 3.5 GB (Cellpose) and 5.8 GB (CellSAM) compressed,
  with SIF conversion on top of the download. A CI run was lost to a pull being
  SIGTERMed part-way through (exit 143) and reported as a task failure. CI also
  now caches the converted SIFs between runs, so a green pipeline no longer
  depends on re-pulling ~6 GB from the registries every time.

- **:warning: `cellsam_min_area` now defaults to `200`, matching
  `cellpose_min_area` (was `0`).** CellSAM ran with no size floor at all while
  Cellpose dropped anything under 200 px², so the same slide was filtered
  differently depending only on which segmenter produced it. **This changes
  CellSAM results:** objects below 200 px² — 15.7 µm² at COMET's 0.28 µm/px,
  about a 4.5 µm disc — are now discarded as debris. Pass
  `--cellsam_min_area 0` to restore the old behaviour, and don't pool CellSAM
  results from before and after.

  The two are directly comparable but not identical: both are in pixels² and
  both keep cells with area >= the threshold, but CellSAM counts the label's
  pixels while sopa measures the vectorised, smoothed polygon. Measured across
  round and ragged cells the polygon is 0.93–1.00× the pixel count (median
  0.99), so 200 cuts at an equivalent ~202 px on the Cellpose side.
  `mesmer_min_nuclei_area` is deliberately left alone — it filters nuclei, not
  whole cells.

- **The Cellpose tuning parameters say what they do.** All four were previously
  described by restating their own name ("Flow threshold for cellpose"), which
  gave no indication of which way to turn them. The schema and `docs/usage.md`
  now give the direction of effect: `cellpose_flow_threshold` higher finds more
  cells and `0` disables the check outright, `cellpose_cellprob_threshold`
  lower finds more _and larger_ cells — and so moves cell areas and every
  intensity measured over them, unlike the flow threshold — and
  `cellpose_diameter` is only a resize factor (`rescale = 30 / diameter`), so
  the default of 30 resamples nothing and filters nothing.

- **Container images and model weights can be kept off the home directory.** A
  run on a cluster failed after an hour with
  `Failed to pull singularity image ... disk quota exceeded`: Apptainer unpacks
  the OCI layers it converts to SIF into `~/.apptainer/cache`, tens of GB across
  this pipeline's images, and a home quota does not hold them. A new
  `container_cache_dir` parameter sets `apptainer.cacheDir`/`singularity.cacheDir`
  so the converted images land on scratch.

  That covers only half of a pull. The layer cache is `$APPTAINER_CACHEDIR`,
  which Nextflow reads from the shell that launched it — image pulls run in
  Nextflow's own process, not in a task, so no config file can redirect it.
  Setting `container_cache_dir` without exporting it is now a startup error
  naming the exports to run, and using Apptainer with the layer cache still
  under `$HOME` warns. Both fire before any download starts.
  `docs/usage.md` has a "Caches on a cluster" section covering all four caches.

- **`deepcell_cache_dir` is now bind-mounted into the container.** It had never
  been, and `autoMounts` binds only `$HOME`, `/tmp` and the work directory — so
  pointing it at scratch redirected `$HOME` inside the container to a path that
  was not there. The 1.7 GB CellSAM download went to the container's ephemeral
  layer and every task repeated it, which is the cost the parameter exists to
  avoid. `MESMERSEGMENT` was affected the same way, with its model-download
  `flock` file container-local and therefore serialising nothing.

- **Cellpose is upgraded to 4.2.1.1 and now runs on the GPU.** The container
  moves to a Wave build carrying sopa 2.2.6, cellpose 4.2.1.1 and Meta's
  `dinov3`, replacing `sopa:2.1.11-cellpose` (cellpose 4.0.8). Segmentation
  passes `--gpu` and `SOPA_SEGMENTATIONCELLPOSE` carries `process_gpu`; sopa
  warns that cellpose >=4 "can be slow without a GPU", and the pipeline had
  never passed the flag.

  **The upgrade on its own does not move results.** Holding the model at
  `cpsam` across the old and new stacks gives foreground IoU 0.9987 with every
  cell matched above 0.90 and mean cell area identical to one decimal place. Any
  difference you see comes from the model you select.

- **`cellpose_pretrained_model` selects the model and now defaults to
  `cpsam_v2`** (was `null`). It accepts a built-in name — `cpsam_v2`, `cpsam`,
  `cpdino`, `cpdino-vitb` — or a path to a custom model. Against the previous
  pipeline, `cpsam_v2` gives IoU 0.9552 on the same 119 cells and `cpdino`
  gives 0.9251 while finding 123.

- **`cellpose_model_type` is removed and now raises an error.** Cellpose 4.0.1
  and later ignore `--model-type`, and sopa only reads it when
  `pretrained_model` is unset — which never happens on v4, because sopa defaults
  it to `cpsam`. The `cyto3` default was therefore discarded on every run, and
  results attributed to cyto3 were cpsam. Failing is better than repeating that
  silently.

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

- **The nf-test snapshot suite has been repaired.** Eight tests across five
  files had not passed in any current environment, and unmodified `dev` failed
  them identically. Two causes: `versions.yml` md5s recorded in an environment
  that no longer matches (one historical value is provably a `versions.yml`
  where the sopa version rendered empty), and the `sopa_segment_compartment`
  snapshots predating `patch_width_pixel = 250` — they recorded one patch where
  the test profile now produces nine, so every downstream parquet, shape and
  mask disagreed. Each snapshot was regenerated and then reproduced by a
  separate run without `--update-snapshot`, which is the only pass that proves
  reproducibility.

- **`sopa_patchifyimage - zarr` runs again.** Its setup block called
  `SOPA_CONVERT` with one input after that process gained a second in
  `b015b6f`, so it had been dying on an arity error. Nothing caught it because
  the snapshot suite has never run in CI.

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

- **`cellpose_normalize_global`, which normalises intensity globally across the
  whole image instead of per patch, removing the residual tiling grid.** Even
  with sopa's CLAHE and gaussian disabled, Cellpose still rescales every patch to
  that patch's own 1–99 percentile, so an identical cell lands on a different
  scale in a dense versus a sparse patch and a faint grid persists at patch seams.
  Measured on one identical synthetic cell across patches of differing content,
  its post-normalisation peak varied by 284% under per-patch percentile and 1.8%
  under a shared global scale. When enabled, the pipeline computes one
  `[low, high]` pair per channel from the whole image (at the new
  `cellpose_normalize_percentiles`, default `1,99`) and passes it to Cellpose as
  `--method-kwargs '{"normalize": {"lowhigh": ...}}'`, ordered to match the
  channels Cellpose receives. Off by default; assumes `cellpose_clip_limit = 0`
  and supersedes `cellpose_tile_norm_blocksize`. The Cellpose normalize dict is
  now single-sourced in the `SOPA_SEGMENTATIONCELLPOSE` module rather than
  `conf/modules.config`, since the global bounds are computed per sample at run
  time.

- **`kronos_exclude_markers`, which withholds channels from KRONOS2 so they
  contribute nothing to the embedding.** A COMET panel's `Autofluorescence`
  channel is not a marker and is not in the 288-marker vocabulary, so the run
  failed with only one escape -- `kronos_allow_novel_defaults` -- which still
  _feeds_ the channel to the model, normalised with default statistics, the
  opposite of what is wanted. This is scoped to the embedding step and deletes
  nothing: the reader opens only the kept channel indices, so a withheld channel
  is never read, while it stays in the image, remains available to every other
  process, and keeps the intensity measurements CELLMEASUREMENT wrote for it.
  (Contrast `remove_markers`, which does strip channels out of the
  background-subtracted image.) Names are comma separated and case sensitive,
  matched against the image's own channel names before any
  `kronos_marker_mapping`; a name matching nothing fails the run and lists what
  is present, and the nuclear channel cannot be withheld because it is the
  model's `preferred_dapi`. The marker report shows withheld channels in place.
  **KRONOS2 attends across the channel set, so this changes every channel's
  embedding** -- hold the list constant across a cohort.

- **The `cpdino` and `cpdino-vitb` models are usable.** Cellpose needs Meta's
  `dinov3` for these, imports it in a bare `try`/`except` that only warns, and
  then dies with `NameError: name 'dinov3_vitl16' is not defined` when the model
  is built — after fetching a 1.13 GB checkpoint. The container now carries
  `dinov3`, pinned to commit `6876159a`; it has no releases or tags, so tracking
  `main` would have been an unpinned dependency. Only its architecture code is
  used, because cellpose builds the backbone with `pretrained=False` and loads
  its own weights, so Meta's gated backbone weights are not required. See
  `CITATIONS.md` for the licence terms, which apply to redistributing the
  container.

- `cellpose_models_dir` parameter, pointing at a pre-staged Cellpose model
  cache so sites with no outbound network on compute nodes skip the download.

- A `CELLPOSE_MODELS_DIR` hook in `tests/nextflow.config`. Export it and the
  sopa tests reuse a pre-staged Cellpose cache instead of re-fetching 1.2 GB per
  run, which is what makes the suite tractable to run locally.

- `docs/nf-test-status.md`, recording the state of the test suite — what was
  repaired, and for the rest what is known versus assumed. Notably, whether
  CellSAM is deterministic on CPU is **unresolved**: `accelerator = 1` in
  `conf/base.config` makes Nextflow pass `--gpus all`, so runs that looked like
  CPU runs were not.

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

- **:warning: KRONOS2 embedded one channel's pixels under every marker's name on
  OME-TIFFs that store each channel as its own 2D series.** `ChannelReader` sized
  itself from `series[0]` alone, so an 8-series image read as 1 channel while the
  patch buffer was still sized from the 8-name marker list -- and numpy broadcast
  that single channel across all eight. The run completed, reported the right
  channel count, and produced 768-d vectors computed from eight copies of DAPI.
  The layout is now detected up front and each channel read from its own series,
  the reader refuses to return fewer channels than were asked of it, and the run
  fails if the channel-name count and the image's channel count disagree.
  **Embeddings produced from such images by earlier versions are invalid and must
  be regenerated.** The COMET `_fixed` images written by the inForm unmixing
  export are affected.
- A channel whose OME `Channel` element carried neither `Name` nor `ID` was
  dropped from the marker list but not from the image the reader opened, sliding
  every later marker onto the wrong channel. Names are now built positionally,
  with one entry per channel always.
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
