# WEHI-SODA-Hub/sp_segment: Usage

## Introduction

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated or YAML file with a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

Ensure that each row has a unique `sample` name. Even if you want to run two both Mesmer and Cellpose on the same sample, you will have to give them different sample names to run both methods without collisions.

If you only want too perform background subtraction, the minimal sample sheet is:

```csv
sample,run_backsub,tiff
sample1,true,/path/to/sample1.tiff
sample2,true,/path/to/sample2.tiff
```

### Full samplesheet

A full sample sheet is shown below:

```csv
sample,run_backsub,run_mesmer,run_cellpose,run_cellsam,tiff
sample1,true,true,false,false,/path/to/sample1.tiff
sample2,true,false,true,false,/path/to/sample2.tiff
sample3,false,false,false,true,/path/to/sample3.tiff
```

You may also prefer to use YAML for your samplesheet, either is supported:

`samplesheet.yml`:

```yaml
- sample: sample1
  run_backsub: true
  run_mesmer: true
  run_cellpose: false
  run_cellsam: false
  tiff: /path/to/sample1.tiff
- sample: sample2
  run_backsub: true
  run_mesmer: false
  run_cellpose: true
  run_cellsam: false
  tiff: /path/to/sample2.tiff
- sample: sample3
  run_backsub: false
  run_mesmer: false
  run_cellpose: false
  run_cellsam: true
  tiff: /path/to/sample3.tiff
```

| Column         | Description                                                                                                           |
| -------------- | --------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Custom sample name.                                                                                                   |
| `run_backsub`  | Run background subtraction on the image                                                                               |
| `run_mesmer`   | Run Mesmer segmentation on the image (only one of `run_mesmer`, `run_cellpose`, `run_cellsam` can be true per row).   |
| `run_cellpose` | Run Cellpose segmentation on the image (only one of `run_mesmer`, `run_cellpose`, `run_cellsam` can be true per row). |
| `run_cellsam`  | Run CellSAM segmentation on the image (only one of `run_mesmer`, `run_cellpose`, `run_cellsam` can be true per row).  |
| `tiff`         | OME-TIFF for COMET or multi-channel TIFF from MIBI                                                                    |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run WEHI-SODA-Hub/sp_segment \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run WEHI-SODA-Hub/sp_segment -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: "./samplesheet.csv"
outdir: "./results/"
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Background subtraction parameters

| Parameter Name | Description                                                     |
| -------------- | --------------------------------------------------------------- |
| remove_markers | Marker channels to remove from the background subtracted image. |

### Cell processing parameters

| Parameter Name      | Description                                                                        |
| ------------------- | ---------------------------------------------------------------------------------- |
| use_whole_cell_only | Use only the whole-cell segmentation to process cells (skip nuclear segmentation). |

### Combine channel parameters

Works for both mesmer and sopa segmentation.

| Parameter Name | Description                                                |
| -------------- | ---------------------------------------------------------- |
| combine_method | Method used to combine membrane channels (max or product). |

### Mesmer parameters

The following Mesmer parameters can be set:

| Parameter Name             | Description                                                                                                                                              |
| -------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| mesmer_segmentation_level  | Segmentation level (legacy parameter).                                                                                                                   |
| mesmer_maxima_threshold    | Controls segmentation level directly in mesmer, (lower values = more cells, higher values = fewer cells).                                                |
| mesmer_interior_threshold  | Controls how conservative model is in distinguishing cell from background (higher values = more conservative interior therefore smaller nuclei & cells). |
| mesmer_maxima_smooth       | Smooths signal peaks (higher values = less irregular shapes/nuclei).                                                                                     |
| mesmer_min_nuclei_area     | Minimum area of nuclei to keep in square pixels.                                                                                                         |
| mesmer_remove_border_cells | Remove cells that touch the image border.                                                                                                                |
| mesmer_pixel_expansion     | Manual pixel expansion after segmentation.                                                                                                               |
| mesmer_padding             | Number of pixels to crop the image by on each side before segmentation.                                                                                  |

### Cellpose parameters

| Parameter Name              | Default  | Description                                                                                  |
| --------------------------- | -------- | -------------------------------------------------------------------------------------------- |
| cellpose_diameter           | 30       | How wide your cells are, in px. Rescales the image so they reach the 30 px Cellpose targets. |
| cellpose_min_area           | 200      | Discard cells below this many px². **Higher = fewer cells**; 0 keeps every mask.             |
| cellpose_flow_threshold     | 0.4      | Max flow error for a mask to survive. **Higher = more cells** (some ill-shaped); 0 disables. |
| cellpose_cellprob_threshold | 0.0      | Cell-probability cutoff. **Lower = more and larger cells**; higher = fewer and smaller.      |
| cellpose_pretrained_model   | cpsam_v2 | Model to segment with: a built-in name or a path to a custom model.                          |
| cellpose_models_dir         | null     | Directory of pre-staged Cellpose weights; skips the download entirely.                       |

#### Tuning the thresholds

The two thresholds do different jobs, and only one of them changes your
measurements:

- **`cellpose_flow_threshold` rejects on shape.** Nothing stops the network
  predicting flows that correspond to no real shape, and where it is uncertain
  it sometimes does. So Cellpose recomputes the flow gradients each finished
  mask implies, takes the mean squared error against the flows the network
  predicted, and drops masks whose error exceeds the threshold. Raise it
  (0.6, 0.8) when cells you can plainly see are being missed; lower it
  (0.2, 0.3) when you are getting ragged blobs. Setting `0` switches the check
  off and keeps everything — Cellpose guards the step with
  `if flow_threshold > 0`. Not used for 3D data.

- **`cellpose_cellprob_threshold` moves the boundaries.** It is the cutoff on
  the third network output, cell probability, and the pixels above it are the
  ones fed into the dynamics step that forms ROIs — so it decides not only
  which cells exist but where each one ends. Lower it (−1, −2) to pick up faint
  cells and grow the ones you have; raise it (1, 2) to stop segmenting dim
  background and shrink them. Because cell areas move, so does every per-cell
  intensity measured over those masks — treat a change here as a change to your
  results, not just to your cell count.

  It is **not a 0–1 probability**: the values are inputs to a sigmoid centred at
  zero, so they run from about **−6 to +6** and 0 is the midpoint rather than
  the floor. Both directions are meaningful and the useful range is wide.

`cellpose_diameter` tells Cellpose the scale of your data. Cellpose resizes the
image by `rescale = 30 / diameter` before the network sees it — 30 being the
cell diameter the model targets — and resizes the masks back afterwards, so a
cell you declare as 60 px is presented to the network at 30. The default of 30
asserts your cells are already at that scale, which is why it happens to
resample nothing; that is the consequence, not the meaning.

Cellpose-SAM tolerates diameters from 7.5 to 120 px (mean 30), so unlike earlier
Cellpose versions this is optional for typical data and getting it somewhat
wrong is not fatal. It still matters if your cells sit outside that range, and
Cellpose's own documented use is downsampling very large cells: `90` shrinks the
image 3× and speeds the run up, at the cost of resolving small objects. It
filters nothing — use `cellpose_min_area` for that.

Note `cellpose_min_area` is in **pixels², not microns²**, so it does not track
`pixel_size_microns`. At COMET's 0.28 µm/px the default 200 px² is 15.7 µm²,
about a 4.5 µm disc; on a different pixel size the same number means a
different biological size. See [Minimum cell area](#minimum-cell-area) — the
same floor applies on the CellSAM path.

Intensity preprocessing has its own parameter group, **Cellpose preprocessing
options** — see [Preprocessing](#preprocessing-the-pipeline-matches-cellposes-own-contract)
below.

#### Most of these do not need adjusting on the SAM and DINO models

Every model this pipeline offers — `cpsam_v2`, `cpsam`, `cpdino`, `cpdino-vitb` —
is a Cellpose 4 transformer model, and the defaults below are already the values
its authors use at test time. **The expected workflow is to change nothing and
pick a model.** Tuning inherited from Cellpose 2/3 habits mostly does not apply.

| Parameter                      | Needs tuning?                         | Why                                                                                                                                                                                                        |
| ------------------------------ | ------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cellpose_diameter`            | **No** — leave at 30                  | Cellpose-SAM trained on scales 0.25×–4× around a 30 px mean diameter, so it is scale-robust over ~7–120 px. 30 means `rescale = 1.0`, no resize.                                                           |
| `cellpose_flow_threshold`      | **No** — 0.4 is the published default | The paper's test-time value. Raise only if you see fragmented or spurious masks.                                                                                                                           |
| `cellpose_cellprob_threshold`  | **No** — 0.0 is the published default | The paper's test-time value. Lower to find more and larger masks, raise for fewer.                                                                                                                         |
| `cellpose_gaussian_sigma`      | **No** — 0                            | sopa's addition; Cellpose applies no smoothing by default.                                                                                                                                                 |
| `cellpose_clip_limit`          | **No** — 0                            | sopa's CLAHE, off. See the preprocessing section.                                                                                                                                                          |
| `cellpose_clahe_kernel_size`   | **No** — inert while CLAHE is off     | —                                                                                                                                                                                                          |
| `cellpose_tile_norm_blocksize` | **No** — 0 is Cellpose's own default  | Only worth touching for a specific normalisation experiment.                                                                                                                                               |
| `cellpose_model_type`          | **Removed** — setting it is an error  | Cellpose 4 ignores `--model-type` entirely.                                                                                                                                                                |
| `cellpose_pretrained_model`    | **Yes** — this is the real choice     | Which model runs. See [Choosing a model](#choosing-a-model).                                                                                                                                               |
| `cellpose_min_area`            | **Maybe**                             | A pipeline-side size filter in px², not a model parameter. Set it for your cell type.                                                                                                                      |
| `patch_width_pixel`            | **Only for memory/throughput**        | Cellpose tiles internally at 256 px regardless, so this controls task size and parallelism. With sopa's preprocessing disabled it has much less effect on results — see the reproducibility numbers below. |

There is also **no automatic diameter estimation** in Cellpose 4: the paper states
plainly that "there was no diameter estimation and resizing performed like in
previous versions of Cellpose". If you know the QuPath Cellpose extension, note
that its `.diameter(0.0)` "automatic computation" belongs to Cellpose 2/3, which
shipped a `SizeModel`. That class no longer exists, and on Cellpose 4 a diameter
of `0` silently means "do not rescale" rather than "estimate it" — so carrying the
habit across produces a plausible but unintended run. Leave it at 30.

`cellpose_model_type` has been **removed**. Cellpose 4 ignores `--model-type`
(it logs `model_type argument is not used in v4.0.1+`), and sopa only reads it
on the cellpose<4 code path, so the old `cyto3` default was silently discarded
on every run — results attributed to `cyto3` were actually `cpsam`. Setting it
now raises an error rather than being ignored again. Use
`cellpose_pretrained_model`.

#### Choosing a model

| Model         | Notes                                                                              |
| ------------- | ---------------------------------------------------------------------------------- |
| `cpsam_v2`    | **Default.** The current Cellpose-SAM model.                                       |
| `cpsam`       | The previous Cellpose-SAM model. Use for continuity with pre-cellpose-4.2 results. |
| `cpdino`      | DINOv3 backbone (ViT-L, ~303M parameters), added in cellpose 4.2.                  |
| `cpdino-vitb` | Smaller ViT-B DINOv3 variant.                                                      |

Anything that is not one of those names is treated as a path to your own
trained model, and is passed to cellpose untouched.

How far apart are they? Measured on the synthetic COMET test image (119 cells),
against the pre-upgrade pipeline:

| Run                              | cells | foreground IoU | cells at IoU >= 0.90 |
| -------------------------------- | ----- | -------------- | -------------------- |
| `cpsam` on the upgraded pipeline | 119   | 0.9987         | 119/119 (100%)       |
| `cpsam_v2`                       | 119   | 0.9552         | 112/119 (94%)        |
| `cpdino`                         | 123   | 0.9251         | 106/119 (89%)        |

The first row is the useful one: upgrading cellpose, sopa and the device while
holding the model fixed is effectively a no-op, so any difference you see comes
from the model you choose, not from the upgrade. Note this is a small synthetic
image — enough to show the models differ, not enough to say which is better for
your tissue.

The DINO models need Meta's `dinov3` package. The pinned container carries it,
and only its architecture code is used, so Meta's gated backbone weights are not
required. See `CITATIONS.md` for the licence terms, which apply if you
redistribute the container or publish results using those models.

#### Preprocessing: the pipeline matches Cellpose's own contract

Cellpose normalises each image it is given by rescaling to the 1st and 99th
percentiles, and does nothing else — its `sharpen_radius` and `smooth_radius`
both default to `0`, and it contains no CLAHE code at all.

That is not an inference from the code alone — it is the documented training
recipe. The Cellpose-SAM paper states that "during training, all images were
normalized such that 0 was set to the first percentile of the image intensity and
1 was the 99th percentile", and that at test time "there was no diameter
estimation and resizing performed like in previous versions of Cellpose".

Sopa adds two preprocessing stages of its own on top of that, applied to every
patch before Cellpose sees it: a `scipy` gaussian filter (`gaussian_sigma=1`)
and CLAHE (`skimage.exposure.equalize_adapthist`, `clip_limit=0.2`). **This
pipeline disables both**, so Cellpose receives exactly what its released weights
were trained on:

| Stage                                    | Whose    | Default here   |
| ---------------------------------------- | -------- | -------------- |
| `gaussian_filter(sigma)`                 | sopa     | **off** (`0`)  |
| `equalize_adapthist(clip_limit, kernel)` | sopa     | **off** (`0`)  |
| percentile (1, 99) rescale               | Cellpose | on (unchanged) |

These four parameters form the **Cellpose preprocessing options** group:

| Parameter Name               | Default | Description                                                                            |
| ---------------------------- | ------- | -------------------------------------------------------------------------------------- |
| cellpose_gaussian_sigma      | `0`     | sopa's gaussian smoothing sigma. `0` disables it; `null` restores sopa's `1`.          |
| cellpose_clip_limit          | `0`     | sopa's CLAHE clip limit. `0` disables CLAHE; `null` restores sopa's `0.2`.             |
| cellpose_clahe_kernel_size   | `null`  | CLAHE kernel in pixels. `null` means `patch.shape // 8`. No effect while CLAHE is off. |
| cellpose_tile_norm_blocksize | `0`     | Cellpose tile-local normalisation block. `0` is Cellpose's whole-patch default.        |

To restore sopa's behaviour, set `--cellpose_clip_limit 0.2` and
`--cellpose_gaussian_sigma 1`.

##### Why CLAHE is off

Segmentation always runs on patches, and all of this normalisation is computed
**from the patch**. A cell's input therefore depends on what else landed in its
patch. Patch overlap does not help: it resolves cells straddling a boundary, not
cells normalised differently.

CLAHE was by far the largest contributor, for a specific reason: left at
skimage's default, the CLAHE kernel is `patch.shape // 8`, so the
contrast-enhancement length scale moves with `--patch_width_pixel`. At COMET's
0.28 µm/px that is 52 µm at `patch_width_pixel=1500`, 105 µm at 3000 and 210 µm
at 6000 — roughly 6 to 25 nuclear diameters for identical tissue, decided
entirely by a patching knob.

Cellpose-SAM's training augmentations are relevant to why this particular stage
matters and the rest do not. Training jittered each channel's brightness (pixel
mean, σ = 0.2) and rescaled its contrast (standard deviation), and degraded 50%
of each batch with Poisson noise, gaussian blur, downsampling and anisotropic
blur. So the model is explicitly robust to **global** intensity and contrast
shifts, and to blurring — which is why sopa's gaussian is harmless, and why
patch-to-patch differences in a single percentile rescale are tolerable. CLAHE is
different in kind: it is a **spatially varying** local remapping, which is not in
that augmentation set at any scale.

Measured on a real 34-channel COMET image (DAPI, 32934×18076) by holding a
1024² region fixed and varying only the patch layout, then recording how much
the values Cellpose receives for those same physical pixels changed. Mean
per-pixel spread, worst of two axes (patch size varied at fixed origin; patch
origin shifted at fixed size). A gaussian-only control measured exactly 0.0000,
which is what makes the rest meaningful:

| Nuclei/mm² | DAPI+ | sopa preprocessing | Cellpose contract |
| ---------- | ----- | ------------------ | ----------------- |
| 2,250      | 10.1% | 0.1211             | **0.0265**        |
| 2,932      | 14.6% | 0.0651             | **0.0052**        |
| 4,647      | 24.2% | 0.0817             | **0.0075**        |
| 12,006     | 66.2% | 0.0758             | **0.0553**        |

Sparse tissue was the worst case, which matters because a whole-slide run covers
every density.

##### What that costs you in cells

The above measures the network's input. The consequence is measurable in the
output. Holding a 1024² region fixed and segmenting it as part of a 1500 px
versus a 3000 px patch, then matching cells across the two layouts by IoU —
i.e. **how many of your cells survive changing `patch_width_pixel` and nothing
else**:

| Model                | sopa: IoU ≥ 0.90 | Cellpose: IoU ≥ 0.90 | sopa mean IoU | Cellpose mean IoU |
| -------------------- | ---------------- | -------------------- | ------------- | ----------------- |
| `cpsam_v2` (default) | 103/267 (38.6%)  | **267/326 (81.9%)**  | 0.869         | 0.928             |
| `cpsam`              | 155/359 (43.2%)  | **271/323 (83.9%)**  | 0.860         | 0.929             |
| `cpdino`             | 74/201 (36.8%)   | **317/366 (86.6%)**  | 0.855         | 0.936             |

Under sopa's preprocessing fewer than half of all cell boundaries survive a
patch-size change, and it is worst for `cpsam_v2`, the default. All three models
land at 82–87% once Cellpose's contract is used, with mean IoU rising from ~0.86
to ~0.93. Three architectures, same direction, similar magnitude.

##### What else changes

Disabling sopa's preprocessing is **not** result-neutral. On dense patches:

| Model      | Cell count       | Median cell area |
| ---------- | ---------------- | ---------------- |
| `cpsam_v2` | **+8% to +35%**  | −23% to −43%     |
| `cpsam`    | **−3% to −10%**  | −25% to −37%     |
| `cpdino`   | **+28% to +54%** | −25% to −42%     |

The direction of the count change is model-dependent — CLAHE suppresses
detections in the newer models and inflates them in `cpsam` — but **median cell
area falls 23–43% for every model**. Since per-cell marker intensities are
computed over these masks, downstream measurements shift accordingly, and results
produced before and after this change should not be pooled in one analysis.

Cells now come out smaller — measured equivalent diameters of ~24–27 px against
the configured 30 — but `cellpose_diameter` almost certainly does **not** need
changing for that. Cellpose-SAM was trained with images resized by a scale factor
log-distributed between 0.25× and 4× relative to a 30 px mean cell diameter, so a
1.2× scale discrepancy sits well inside the range it was trained to handle. The
default of 30 gives `rescale = 1.0`, i.e. no resizing, which is what the paper's
test-time protocol used.

A secondary effect, worth knowing but **not** a reason on its own: CLAHE
bypasses Cellpose's empty-patch guard. Cellpose zeroes a patch whose 1–99
percentile range falls below `1e-3`, so an off-tissue patch reaches the network
as zeros. `equalize_adapthist` expands an all-zero patch to the full range
(mean 0.928, std 0.258), so that check never fires. On the image above, 13 of 84
patches at `patch_width_pixel=3000` are entirely zero — unremarkable for exported
ROIs with off-tissue padding.

In practice this appears to be harmless: run on such a patch, `cpsam` returned
**0 masks either way**, so the amplified field does not manufacture cells. Treat
it as wasted work and a lost safety net rather than a correctness problem.

##### Visual check on real tissue

These defaults were run through the pipeline on a COMET brain-panel ROI
(1683×2349 at 0.28 µm/px, DAPI plus four markers combined by `max`) with
`cpsam_v2`, `cpdino` and cellSAM, and the resulting GeoJSONs were inspected in
QuPath. All three loaded cleanly and looked closely comparable:

| Method              | Cells | Median area        | Median equiv. diameter |
| ------------------- | ----- | ------------------ | ---------------------- |
| cellpose `cpsam_v2` | 1,742 | 864 px² (68 µm²)   | 33.2 px (9.3 µm)       |
| cellpose `cpdino`   | 1,751 | 916 px² (72 µm²)   | 34.2 px (9.6 µm)       |
| cellSAM             | 1,648 | 1,129 px² (89 µm²) | 37.9 px (10.6 µm)      |

Three architecturally unrelated models landing within 6% on cell count is decent
mutual corroboration. They differ on cell _size_: cellSAM draws noticeably larger,
more generous boundaries around membrane signal — visible by eye in QuPath and
matching the ~30% larger median area. `cpsam_v2` and `cpdino` are close to
indistinguishable from one another.

This is a plausibility check, not an accuracy measurement — there were no
annotations to score against, so it establishes that the defaults produce sensible
segmentation on real tissue, not that any one model is the most correct.

> **All of this measures reproducibility, not accuracy.** It shows that the same
> tissue yields the same cells regardless of how it happened to be patched. It
> does **not** show the cells are more correct — a config could be perfectly
> self-consistent and still segment badly, and nothing here was compared against
> annotations. Sopa presumably added CLAHE because dim IF channels benefit from
> it, and that trade is unevaluated. If your panel has weak markers, compare both
> settings before trusting either.
>
> **Judging accuracy needs visual inspection**, against annotations or by eye in
> QuPath, and no amount of normalisation measurement can substitute for it. The
> defaults have been eyeballed on real tissue across three models (see the visual
> check above) and produce sensible boundaries, but that is a plausibility check
> rather than a scored comparison. Treat them as the reproducible starting point
> and confirm the boundaries look right on your own panel. Both stages are one flag
> each to restore, so the comparison is cheap to run.
>
> Measured on nuclear (DAPI) segmentation on one image. The whole-cell compartment
> combines membrane markers, which are dimmer and sparser — precisely where CLAHE
> would be earning its keep — and was not tested.
>
> Worth knowing for context: cellSAM also CLAHEs every block and works well —
> but its released model was _trained_ with CLAHE applied, and it pins the kernel
> at a fixed 128 px instead of deriving it from the block size. Cellpose's
> weights were not trained that way. CLAHE is not inherently harmful; applying it
> to a model that never saw it is the problem.

##### Tile-local normalisation

`--cellpose_tile_norm_blocksize` makes Cellpose's percentile statistics local to
blocks within the patch rather than global to it. It defaults to `0`, Cellpose's
own default, meaning whole-patch percentiles.

If you enable it, keep the block much larger than a cell. At 0.28 µm/px with
`--cellpose_diameter 30`, 512 px is 143 µm or ~17 nuclear diameters; 128 px is
36 µm or ~4, small enough that a block in sparse tissue may contain no nuclei at
all. At 128 the minimum per-block 1–99 range fell to 5.0 against 1670 for the
whole patch, and since Cellpose only guards a flat block below `1e-3` it rescued
none of them — noise was stretched to full range and patch dependence got 3–24×
_worse_. Note also that enabling it while CLAHE is on regresses: the two stages
stack rather than cancel.

Neither option removes patch-locality entirely, because the block grid is still
laid out relative to the patch. Cellpose's `lowhigh` would — it applies fixed
absolute bounds with no patch statistics at all — but that is probably not an
improvement here. Cellpose-SAM was trained by normalising **each image** to its
own 1st/99th percentiles and then taking 256×256 crops, so per-patch percentile
normalisation is the closer analogue of the training setup; fixed cohort-wide
bounds would move further from it, not closer. Leave `lowhigh` alone unless you
have a specific reason.

#### GPU

Cellpose segmentation runs on the GPU. Sopa warns that cellpose >=4 "can be slow
without a GPU", and `SOPA_SEGMENTATIONCELLPOSE` carries the `process_gpu` label
so one is requested. On a machine without a GPU, override that label — the
models still run on CPU, just slowly.

#### Cellpose model weights

The pipeline fetches the Cellpose weights **once per run** and stages them as an
input to every segmentation task. Previously each patch task downloaded the
~1.2 GB `cpsam` checkpoint itself — roughly 18 identical requests for a single
sample — which caused HTTP 429 rate limiting and truncated downloads. The staged
copy is also reused by `-resume`.

On a cluster, or anywhere compute nodes have no outbound network, point the
pipeline at a shared copy instead.

> **Do not put the cache under `/opt`.** Nextflow bind-mounts the host
> directory containing a staged input, so `--cellpose_models_dir /opt/...`
> mounts the host's `/opt` over the container's, hiding `/opt/conda` where the
> Cellpose image is installed. Every task then fails with
> `sopa: command not found`. The pipeline rejects this at startup rather than
> letting you discover it that way.

```bash
nextflow run WEHI-SODA-Hub/sp_segment --cellpose_models_dir /shared/cellpose_models ...
```

That directory is the Cellpose model cache. Populate it once on a login node,
naming the model you intend to run:

```bash
CELLPOSE_LOCAL_MODELS_PATH=/shared/cellpose_models \
  python -c "from cellpose import models; models.CellposeModel(gpu=False, pretrained_model='cpsam_v2')"
```

Note that this pins _where_ the weights come from, not _which_ weights: the
container is pinned but cellpose.org serves whatever is current, so results are
only reproducible over time against a fixed `--cellpose_models_dir`.

`cellpose_pretrained_model` and `cellpose_models_dir` are both checked before
the pipeline starts any work — the former only as a path when it is not a
built-in model name. Cellpose silently falls back to
its built-in weights when given a `--pretrained-model` path that does not exist,
so a typo previously produced a complete and plausible but wrong run.

### CellSAM segmentation

CellSAM is an optional segmentation backend for whole-cell and nuclear masks.
Enable it per sample with `run_cellsam: true` in the samplesheet.

For gated model downloads, set your DeepCell token as a Nextflow secret:

```bash
nextflow secrets set DEEPCELL_ACCESS_TOKEN $YOUR_TOKEN
```

If no token is provided, CellSAM uses the bundled default model.

#### CellSAM parameters

| Parameter Name                     | Description                                                                                                |
| ---------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `cellsam_bbox_threshold`           | Confidence threshold for bounding-box detections (default: `0.4`).                                         |
| `cellsam_block_size`               | Tile size in pixels used when processing large images (default: `600`).                                    |
| `cellsam_overlap`                  | Overlap in pixels between adjacent tiles (default: `250`).                                                 |
| `cellsam_iou_depth`                | Search depth (pixels) for IoU-based duplicate removal at tile boundaries (default: `250`).                 |
| `cellsam_iou_threshold`            | IoU threshold for non-maximum suppression across tiles (default: `0.5`).                                   |
| `cellsam_use_wsi`                  | Enable whole-slide-image tiling mode (default: `true`).                                                    |
| `cellsam_gauge_cell_size`          | Automatically estimate cell size from the image before segmentation (default: `false`).                    |
| `cellsam_low_contrast_enhancement` | Apply contrast enhancement before segmentation for low-contrast images (default: `false`).                 |
| `cellsam_model_path`               | Path to a custom CellSAM model checkpoint. If `null` the built-in default model is used (default: `null`). |
| `cellsam_min_area`                 | Minimum cell area in square pixels; smaller objects are discarded (default: `200`, matching Cellpose).     |

#### Minimum cell area

`cellsam_min_area` and `cellpose_min_area` share a default of `200`, so a run is
not filtered differently just because of which segmenter produced it. Both are
in **pixels², not microns²** — at COMET's 0.28 µm/px that is 15.7 µm², roughly a
4.5 µm disc — and both keep cells whose area is greater than or equal to the
threshold. Set either to `0` to keep every mask.

They measure that area slightly differently, so they are close but not exactly
interchangeable:

|                     | Measures                                       | Applied                                      |
| ------------------- | ---------------------------------------------- | -------------------------------------------- |
| `cellsam_min_area`  | the label's pixel count                        | once, after the whole-slide tiles are merged |
| `cellpose_min_area` | the vectorised polygon's area, after smoothing | per patch, during segmentation               |

Comparing the two on round and ragged cells, the polygon comes out at 0.93–1.00×
the pixel count (median 0.99), so `200` cuts at an equivalent ~202 px on the
Cellpose side. The per-patch timing is safe because `patch_overlap_pixel` means
a cell truncated at one patch edge appears whole in its neighbour, and the
full-size copy is the one that survives to be resolved.

Mesmer's `mesmer_min_nuclei_area` is deliberately **not** aligned with these: it
filters nuclei rather than whole cells, and is applied inside the Mesmer
container rather than by this pipeline.

### KRONOS2 embeddings

KRONOS2 is an optional step that runs after cell measurement and writes a
768-dimensional embedding for every cell back into the cellmeasurement GeoJSON.

Cells are taken from the CELLMEASUREMENT GeoJSON polygons, so embedding row _i_
always corresponds to cell feature _i_ -- the join is exact by construction
rather than by spatial matching.

To enable KRONOS2:

```bash
nextflow run WEHI-SODA-Hub/sp_segment \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --enable_kronos true \
   --kronos_model_path /path/to/KRONOS2
```

#### Obtaining the model

The weights are gated on Hugging Face. Request access, then download the
**full** snapshot -- the loader adds the directory to `sys.path` and imports the
bundled `dinov2` package, so a hand-picked subset of files will not work:

```bash
hf auth login
hf download MahmoodLab/KRONOS2 --local-dir /path/to/KRONOS2 --exclude "demo_image/*"
```

**Do not stage the weights under a path the container also uses**, such as
`/opt`. Nextflow bind-mounts an input's parent directory into the container to
make it visible, so a host `/opt` shadows the container's own `/opt/conda` and
its Python disappears. The failure is reported as:

```
/usr/bin/env: 'python3': No such file or directory
```

from an image that demonstrably has `python3`, which points nowhere near the
cause. Somewhere like `/data/KRONOS2` or a project directory is safe. This
applies to any process input, not just KRONOS2 -- `/opt` is simply the common
case, because conda-based images live there.

No marker metadata file is needed: KRONOS2 carries its own 288-marker vocabulary
and applies the marker-aware z-score internally.

#### Marker names must match exactly

KRONOS2 matches marker names on a **separator-insensitive key with no alias or
fuzzy step** (`CD-8`, `CD_8` and `CD8` all match; `Cytokeritin` does not match
`CYTOKERATIN`). An unmatched marker is not dropped -- it is normalised with
_default_ statistics, which silently degrades that channel.

The run therefore **fails by default** when any marker is unmatched, listing the
offending names. Most are naming variants rather than new biology, so the fix is
usually a mapping file:

```json
{
  "Collagen 4": "COLLAGENiv",
  "Cytokeritin": "CYTOKERATIN",
  "VISA": "VISTA"
}
```

Pass it with `--kronos_marker_mapping /path/to/marker_mapping.json`. Mapping is
applied only to the names handed to the model; stored channel names and the
GeoJSON are untouched. See
[examples/kronos_marker_mapping.json](examples/kronos_marker_mapping.json).

The nuclear stain is a special case: it is passed as `preferred_dapi`, which is
KRONOS2's **only** alias mechanism. The samplesheet's `nuclear_channel` is used
automatically, so a slide stained with `DRAQ5` or `Hoechst2` still receives DAPI
statistics.

#### KRONOS2 embedding parameters

| Parameter Name                | Description                                                                                                                                                                             |
| ----------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `enable_kronos`               | Enable the KRONOS2 embedding step (default: `false`).                                                                                                                                   |
| `kronos_model_path`           | **Required** path to a full KRONOS2 snapshot directory (or a Hugging Face repo id).                                                                                                     |
| `kronos_patch_size`           | Side length in pixels of the square patch centred on each cell (default: `64`).                                                                                                         |
| `kronos_batch_size`           | Patches per inference batch (default: `16`). KRONOS2 runs fp32 and cuBLAS picks batch-dependent kernels, so changing this shifts embeddings slightly; `16` reproduces published values. |
| `kronos_num_workers`          | PyTorch DataLoader worker processes (default: `4`).                                                                                                                                     |
| `kronos_max_value`            | Override the intensity divisor. Unset (default) derives it from the image dtype -- `uint8`=255, `uint16`=65535, float=400 -- matching KRONOS2's own scaling factor.                     |
| `kronos_marker_mapping`       | JSON file or inline JSON mapping channel names onto KRONOS2 vocabulary names (default: `null`).                                                                                         |
| `kronos_nuclear_marker`       | Override the nuclear stain used as `preferred_dapi`. Defaults to the samplesheet's `nuclear_channel`.                                                                                   |
| `kronos_isolate_cell`         | Zero pixels outside the target cell so each embedding describes that cell rather than its surrounding neighbourhood (default: `true`).                                                  |
| `kronos_allow_novel_defaults` | Proceed when markers fall outside the vocabulary, accepting default normalisation stats (default: `false`).                                                                             |

### SOPA patching parameters

| Parameter Name      | Description                                                              |
| ------------------- | ------------------------------------------------------------------------ |
| technology          | Image type used for zarr conversion, only `ome_tif` is supported (COMET) |
| patch_width_pixel   | Width and height of image patch in pixels                                |
| patch_overlap_pixel | Number of pixels that image patches will overlap                         |

### Mask smoothing options

| Parameter Name     | Description                                                                                                                                                                                    |
| ------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| smooth_masks       | Enable mask smoothing before cell measurement to reduce polygon complexity (default: `false`). Prevents StackOverflowError in QuPath's GeoJSON export for images with complex cell boundaries. |
| smooth_method      | Smoothing method: `morphological` (close+open with disk kernel, conservative) or `gaussian` (blur+threshold, stronger smoothing). Default: `morphological`.                                    |
| smooth_kernel_size | Kernel size for smoothing. For morphological: disk radius (integer). For gaussian: sigma. Larger values = more smoothing. Default: `2`.                                                        |

### Cell measurement options

| Parameter Name              | Description                                                                                                                    |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| enable_measurements         | Calculate intensity and shape measurements for cell compartments (disabling will decrease execution time)                      |
| percentiles                 | Comma-separated list of percentiles to calculate per channel. Enable measurements must be set to `true` to use this parameter. |
| pixel_size_microns          | Pixel size in microns, use 0.28 for COMET and 0.390625 for MIBI                                                                |
| estimate_cell_boundary_dist | Where no matching membrane ROI exists, expand the nucleus by this many pixels                                                  |
| dist_threshold              | Maximum centroid distance in pixels for matching a nucleus to a whole-cell ROI (default: `10.0`).                              |
| downsample_factor           | Integer downsample factor applied to image and masks before measurement, `1` = disabled (default: `1.0`).                      |
| tile_size                   | Tile size in pixels for measurement image reads (default: `2048`). Aim for 2-4x more tiles than the thread count.              |
| tile_overlap                | Tile overlap in pixels for measurement image reads (default: `200`). Set at least as large as the largest cell diameter.       |
| neighbors                   | Number of nearest neighbours for neighbourhood feature aggregation, `0` = disabled (default: `0`).                             |
| erosion_steps               | Measure intensity in 5 equal-area erosion bins from the cell/nucleus boundary inward (default: `true`).                        |
| expansion_steps             | Measure intensity in 5 equal-area expansion bins within 20 µm outward from the cell boundary (default: `true`).                |
| environment_expansion       | Measure a pericellular 20 µm environment zone around each cell (default: `true`).                                              |
| gzip_geojson                | Gzip-compress the output GeoJSON (produces `.geojson.gz`). Recommended for large whole-slide images (default: `true`).         |
| geometry_checkpoint_dir     | Directory for resumable polygon-extraction checkpoints. Unset = disabled. See below.                                           |
| geometry_batch_size         | Cells per polygon-extraction batch, and the checkpoint granularity (default: `2000`).                                          |

#### Resuming polygon extraction after a killed job

On whole-slide images with millions of cells, polygon extraction is the longest
phase of `CELLMEASUREMENT`. Setting `geometry_checkpoint_dir` makes it write each
completed batch of cells to disk, so a task killed by wall-time or OOM resumes
from the last completed batch rather than starting over:

```bash
nextflow run . -profile medium --geometry_checkpoint_dir /vast/scratch/$USER/sp_segment_checkpoints
```

The directory **must be outside the Nextflow work directory**. A retry runs in a
freshly allocated work directory, so a checkpoint written there is discarded and
the retry restarts from zero. Checkpoints are namespaced per sample and per mask,
and are ignored if the run's settings (tolerance, batch size, mask dimensions,
cell count) differ from the ones that produced them. They are not deleted
automatically — remove the directory once a run has completed successfully.

### Report parameters

| Parameter Name  | Description                         |
| --------------- | ----------------------------------- |
| generate_report | Generate segmentation report for QC |

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull WEHI-SODA-Hub/sp_segment
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [WEHI-SODA-Hub/sp_segment releases page](https://github.com/WEHI-SODA-Hub/sp_segment/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility. Conda is currently not supported for the pipeline.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Not supported for this pipeline.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
