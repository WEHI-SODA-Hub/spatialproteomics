# KRONOS1 → KRONOS2 migration plan

**Branch:** `kronos-v2-update` (from `origin/dev` @ `ef47666`)
**Date:** 2026-07-28
**Background analysis:** `C:\Users\Mikem\CORAL\findings.md`

---

## 1. Decisions (locked)

| Decision | Choice |
|---|---|
| Retirement strategy | **Hard replace** — KRONOS1 script + module deleted; there is only one KRONOS |
| Patch semantics | **Isolate the target cell** — zero all pixels not belonging to the cell |
| Scaling divisor | **Derived from dtype** (uint8→255, uint16→65535, float→400) |
| Nuclear hint | **Plumbed through** to `model.preprocess(preferred_dapi=...)` |
| Novel markers | **Full prepare/pool graph**, built now |
| Modality | **Fluorescence only** (COMET / CODEX / OPAL) — MIBI is out of scope |

Two behaviour changes ship together and are **not** backward-comparable with existing v1 outputs: the model (384-d → 768-d) and the patch content (neighbourhood box → isolated cell). This is accepted — hard replace means there was no comparability to preserve.

---

## 2. What gets simpler — and what does not

*Revised 2026-07-28 after reading the actual model code (`marker_utils.py`, `modeling_kronos2.py`). An earlier draft of this section was wrong about D3; the correction is below.*

- **D1 (unmatched markers dropped)** — **gone.** v2 takes marker **names**; nothing is dropped. Unknown markers fall back to a default `(mean, std)` rather than being removed from the tensor. `--marker-metadata` disappears entirely.

- **D3 (literal marker matching)** — **NOT fixed by v2.** This was wrong in the first draft. KRONOS2 matching is **exact on the match key, by design**, with no alias, isotope-strip, or fuzzy step. From `marker_utils.py`:

  > "Matching is **exact on the match key** — there is no fuzzy / alias step — so distinct markers are never silently conflated (the source of the old `CD24 -> cd4` class of mis-normalizations)."

  `marker_match_key` only collapses separators (`- _ : ( ) space . /` and `α→a`). Verified against the real 288-marker vocabulary:

  | Channel name | match key | in vocab? |
  |---|---|---|
  | `CD8` | `cd8` | yes |
  | `HLA-DR` | `hladr` | yes |
  | `a-SMA` | `asma` | yes |
  | `PanCK` | `panck` | yes |
  | `CD8_141Pr` | `cd8141pr` | **no — novel** |
  | `dsDNA_89Y` | `dsdna89y` | **no — novel** |
  | `Hoechst1` | `hoechst1` | **no — novel** |
  | `Collagen 4` | `collagen4` | **no — novel** |

  **Consequence:** spelling variants, typos, and per-cycle stain names (`Hoechst1`, `DAPI-01`) land in the *novel* bucket. This does not lose data the way v1's dropping did, but it means name normalisation is **our** responsibility, not the model's.

  **Action:** carry the alias step forward rather than deleting it — `kronos_marker_mapping` (user JSON) is **retained**. Applied only to the names handed to KRONOS2, never to stored channel names, matching CORAL's scoped exception.

  **Scope note (2026-07-28):** KRONOS is used on **fluorescence data only** (COMET / CODEX / OPAL), not MIBI. Isotope-tagged names (`CD8_141Pr`) are therefore out of scope, and the proposed `kronos_strip_isotope` parameter is **dropped**. Should MIBI support ever be wanted, the isotope regex `_\d+[a-z]*$` is the one-line addition.

- **D5 (nuclear hint)** — **promoted from nice-to-have to essential.** `preferred_dapi` is the *only* alias mechanism in the entire model:

  ```python
  if clean_marker_name(m) == dapi: key = "dapi"
  else:                            key = self._marker_index.get(marker_match_key(m))
  ```

  Without it, a slide stained with `Hoechst1` or `DRAQ5-02` gets default stats instead of DAPI's. Already in scope; now known to be doing real work.

- **D2 (cell isolation)** and **D4 (dtype divisor)** — unchanged, both ours, both in scope.

**Bonus:** `conf/test_mesmer_kronos.config` references `tests/data/comet/kronos_marker_metadata.csv`, **which does not exist in the repo**. Together with the gated weights, the KRONOS1 path is currently unrunnable end-to-end in CI. v2 needs neither file, so this becomes testable for the first time.

---

## 3. The structural problem: a cohort barrier cannot live in a per-segmenter subworkflow

Today `KRONOSEMBEDDINGS` is invoked **three times**, once inside each of `mesmer_segment`, `sopa_segment`, and `cellsam_segment` (the `_wbacksub` variants only forward its emits).

Novel-marker stats must be pooled **across the whole cohort** — that is what makes `(mean, std)` meaningful. But a `collect()` inside `mesmer_segment` only ever sees Mesmer samples. If sample A runs Mesmer and sample B runs CellSAM, they land in different subworkflows and would pool separately, producing two different normalisations for the same marker.

**Therefore KRONOS must be hoisted out of the segmentation subworkflows into `workflows/sp_segment.nf`.**

Each segmentation subworkflow gains a uniform emit:

```groovy
emit:
kronos_input = ch_kronos_input   // [ meta, tiff, whole_cell_mask, geojson ]
```

and the top level mixes them:

```groovy
ch_kronos_input = MESMER_SEGMENT.out.kronos_input
    .mix(MESMER_SEGMENT_WBACKSUB.out.kronos_input)
    .mix(SOPA_SEGMENT.out.kronos_input)
    .mix(SOPA_SEGMENT_WBACKSUB.out.kronos_input)
    .mix(CELLSAM_SEGMENT.out.kronos_input)
    .mix(CELLSAM_SEGMENT_WBACKSUB.out.kronos_input)

if (params.enable_kronos) {
    KRONOS2(ch_kronos_input)      // new subworkflow, section 5
}
```

This is a better structure independently of the barrier — KRONOS is segmenter-agnostic, and the current arrangement duplicates the same ~30-line block three times.

**Side effect to verify, not assume:** each subworkflow currently does `ch_annotations = KRONOSEMBEDDINGS.out.merged_geojson`, so `SEGMENTATIONREPORT` consumes the embedding-merged GeoJSON. After hoisting, the report will consume the plain `CELLMEASUREMENT` output. That should be fine — the report is a segmentation-QC artefact and does not read embeddings — but confirm before merging. Publishing is unaffected: `KRONOS2EMBEDDINGS` publishes its merged GeoJSON into `${outdir}/cellmeasurement/` under the same filename, so the merged file still wins on disk.

---

## 4. Target process graph

```
                 ch_kronos_input  [meta, tiff, wc_mask, geojson]
                          │
                          ├──────────────────────────────┐
                          │ .collect()                   │ (per sample)
                          ▼                              │
              ┌────────────────────────┐                 │
              │ KRONOS2PREFLIGHT       │  once, CPU      │
              │  · marker_metadata.csv │  NO model       │
              │  · marker_utils.py     │  NO torch       │
              │  · channel names (meta │                 │
              │    only, no pixels)    │                 │
              │  · novel marker list   │                 │
              │  · validate CSV rows   │  ← FAIL FAST    │
              └───────────┬────────────┘                 │
                          │ novel_markers.json           │
              ┌───────────┴───────────┐                  │
        novel = []              novel = [...]            │
              │                       │                  │
              │                       ▼                  │
              │           ┌────────────────────────┐     │
              │           │ KRONOS2PREPARE         │◄────┘
              │           │  per sample, CPU       │
              │           │  no model — pure numpy │
              │           │  (n, μ, s²) per marker │
              │           └───────────┬────────────┘
              │                       │ .collect()
              │                       ▼
              │           ┌────────────────────────┐
              │           │ KRONOS2POOLSTATS       │  once, CPU
              │           │  pooled variance →     │  no model
              │           │  additional_markers.csv│
              │           └───────────┬────────────┘
              │                       │
              └───────────┬───────────┘
                          ▼
              ┌────────────────────────┐
              │ KRONOS2EMBEDDINGS      │  per sample, GPU
              │  register CSV → extract│  + model
              └────────────────────────┘
```

**Why a preflight process:** it moves the "you have a novel marker with no CSV row" error **before** any per-sample work runs, matching CORAL's pre-flight stage (`src/coral/features/prepare.py:152-177`).

**Confirmed: PREFLIGHT needs no model, no weights, no torch, no GPU.** `marker_utils.py` (4 kB) is self-contained pure-pandas — `clean_marker_name`, `marker_match_key`, `build_marker_index`, `load_marker_stats` — and the vocabulary is the plain `marker_metadata.csv` (288 rows). Both are downloadable without the 440 MB checkpoint. So novelty classification uses the model's *own* code without instantiating the model. Vendor these two files (or read them from `kronos_model_path`); never reimplement the matcher.

Staging every TIFF into `PREFLIGHT` is cheap on a shared filesystem (Nextflow symlinks rather than copies) and only channel-name metadata is read — no pixel decode. Mirrors the existing `bin/extract_markers.py` pattern.

**Early-exit path:** when `PREFLIGHT` reports zero novel markers, `PREPARE` and `POOLSTATS` are skipped entirely and `EMBEDDINGS` runs with no CSV. This is CORAL's early return (`prepare.py:157-158`) and will be the common case — the synthetic CI panel (`DAPI`, `CD45`) exercises exactly this path.

---

## 5. Files

### Delete

```
bin/kronos_embeddings.py
modules/local/kronosembeddings/{main.nf,environment.yml,meta.yml,tests/}
conf/test_mesmer_kronos.config
```

### Add

```
bin/kronos2_common.py          # shared: channel names, patch dataset, dtype scaling, IO
bin/kronos2_preflight.py
bin/kronos2_prepare.py
bin/kronos2_pool_stats.py
bin/kronos2_embeddings.py
bin/kronos2_list_markers.py    # vocabulary dump + "is my marker novel?" check
modules/local/kronos2preflight/{main.nf,environment.yml,meta.yml,tests/}
modules/local/kronos2prepare/{...}
modules/local/kronos2poolstats/{...}
modules/local/kronos2embeddings/{...}
subworkflows/local/kronos2/{main.nf,meta.yml,tests/}
conf/test_mesmer_kronos2.config
tests/data/comet/additional_markers.csv      # test fixture
docs/examples/kronos_marker_mapping.json     # worked example, section 6.1
docs/examples/additional_markers.csv         # worked example, section 6.1
# NOTE: no vendored marker vocabulary - generated locally only, see section 11
```

`bin/kronos2_common.py` is the CORAL four-slot scaffold in miniature (`src/coral/features/base.py:14`): everything except `transform`/`forward` is shared. Building it now is what stops a future KRONOS3 from becoming a third 1000-line script.

### Modify

```
workflows/sp_segment.nf                          # hoist KRONOS to top level
subworkflows/local/{mesmer,sopa,cellsam}_segment/main.nf          # remove KRONOS block, add kronos_input emit
subworkflows/local/{mesmer,sopa,cellsam}_segment_wbacksub/main.nf # forward kronos_input instead of embeddings/report
conf/modules.config                              # replace KRONOSEMBEDDINGS block with 4 new ones
nextflow.config                                  # params (section 6)
nextflow_schema.json                             # param schema
README.md, docs/usage.md, docs/output.md, CHANGELOG.md
```

---

## 6. Parameters

Hard replace means one KRONOS, so the `kronos_*` namespace is retained — existing `kronos_patch_size` / `kronos_batch_size` / `kronos_num_workers` settings keep working.

**Removed** (nf-schema `validate_params = true` will reject them with a clear error, which is the migration signal we want):

```
kronos_config_path            # v1 config.json — v2 has none
kronos_marker_metadata        # stats now live in the model
kronos_marker_mapping         # v2 resolves names internally
kronos_distance_threshold     # was dead code even in v1
```

**Retained:**

```groovy
enable_kronos                = false
kronos_model_path            = null      // local dir (preferred) or HF repo id
kronos_patch_size            = 64
kronos_num_workers           = 4
```

**Changed default:**

```groovy
kronos_batch_size            = 16        // was 32 — see section 9
kronos_max_value             = null      // was 65535 — null = derive from dtype
```

**New:**

```groovy
kronos_isolate_cell          = true      // zero non-target pixels
kronos_nuclear_marker        = null      // null = use samplesheet nuclear_channel
kronos_additional_markers    = null      // novel-marker CSV (text columns only)
kronos_stats_region          = 'cells'   // 'cells' | 'image'  — see section 8
kronos_allow_novel_defaults  = false     // true = warn instead of erroring on novel markers
```

No output-format parameter: the merged GeoJSON is the **only** embedding artefact (section 10). The v1 `--output` CSV argument and its `embeddings` output channel are removed, not made optional.

---

## 6.1 Worked example — what the user actually supplies

Grounded in the real 49-marker CODEX cHL panel and the real 288-marker KRONOS2 vocabulary. **45 of 49 matched; 4 did not.**

### The key insight: novel ≠ new biology

All four unmatched markers are **already in the vocabulary under a different spelling**:

| Panel name | match key | in vocab? | actually is | resolution |
|---|---|---|---|---|
| `Collagen 4` | `collagen4` | no | `COLLAGENiv` | **mapping** |
| `Cytokeritin` | `cytokeritin` | no | `CYTOKERATIN` (typo in panel) | **mapping** |
| `VISA` | `visa` | no | `VISTA` | **mapping** |
| `DAPI-01` | `dapi01` | no | `DAPI` | **`--kronos_nuclear_marker`** |

**Zero `additional_markers.csv` rows are needed for this panel.** Getting this wrong is silently harmful, not an error: writing a CSV row for `Cytokeritin` gives it *data-derived* stats and a *live-computed* BioLinkBERT embedding, instead of the *pretrained* stats and embedding it would inherit by mapping to `CYTOKERATIN`. Strictly worse, and invisible.

**Therefore PREFLIGHT must emit ranked suggestions for every novel marker**, and the docs must present mapping as the first resort. Suggestions can reuse CORAL's ranking (`src/coral/markers/normalize.py:157` — substring containment, then shared prefix, then `difflib` ratio), which on this panel returns exactly the right answer for all four. Advisory only, never auto-applied — an auto-applied fuzzy match is the `CD24 → CD4` bug the model authors deliberately designed out.

### Decision tree

```
marker not in vocabulary
   ├─ a near-match exists (preflight suggests it)   -> kronos_marker_mapping   [pretrained stats + embedding]
   ├─ it is the nuclear stain (DAPI-01, Hoechst2)  -> kronos_nuclear_marker   [maps to dapi stats]
   └─ genuinely new biology (FoxA1, NaKATPase)     -> additional_markers.csv  [data-derived stats + live embedding]
```

### 1. Samplesheet — unchanged

```csv
sample,run_backsub,run_mesmer,tiff,nuclear_channel,membrane_channels
cHL_01,true,true,/data/cHL_01.ome.tiff,DAPI-01,Cytokeritin:CD45
```

`nuclear_channel` now does double duty: the segmenter input, *and* the default `preferred_dapi` for KRONOS2.

### 2. Params

```groovy
params {
    enable_kronos             = true
    kronos_model_path         = '/data/models/KRONOS2'      // full snapshot minus demo_image/
    kronos_marker_mapping     = 'docs/examples/kronos_marker_mapping.json'
    kronos_additional_markers = null                        // this panel needs none
    kronos_patch_size         = 64
    kronos_batch_size         = 16
    kronos_isolate_cell       = true
}
```

### 3. `kronos_marker_mapping` — a JSON file or inline string

```json
{
  "Collagen 4":  "COLLAGENiv",
  "Cytokeritin": "CYTOKERATIN",
  "VISA":        "VISTA"
}
```

Applied **only** to the names handed to KRONOS2; the stored channel names and the GeoJSON are untouched. See `docs/examples/kronos_marker_mapping.json`.

### 4. `additional_markers.csv` — only for genuinely new biology

```csv
marker_name,compartment,family,compartment_desc,family_desc,marker_full_name,mean,std,pretraining
FoxA1,nucleus,transcription_factor,nuclear protein,Transcription factor,Forkhead box protein A1,,,no
NaKATPase,membrane,metabolic_transporter,membrane protein,Metabolic transporter,Sodium/potassium-transporting ATPase subunit alpha-1,,,no
```

- `mean` / `std` **left blank** — POOLSTATS fills them from your pixels
- `compartment` and `family` **must** come from the existing sets (section 7) or `register_additional_markers` raises
- must contain **only** not-yet-registered markers (partial re-registration raises)

Both example files are committed under `docs/examples/` and validated against the real vocabulary and category sets.

### 5. What PREFLIGHT prints

```
KRONOS2 preflight — 1 sample, 49 channels
  vocabulary: 288 markers (268 pretraining, 20 stats-only)
  matched:    45/49

  NOVEL (4) — resolve each before extraction:
    Collagen 4   -> suggestions: COLLAGEN, COLLAGENiv        [map?]
    Cytokeritin  -> suggestions: CYTOKERATIN, pancytokeratin [map?]
    VISA         -> suggestions: VISTA                       [map?]
    DAPI-01      -> suggestions: DAPI                        [nuclear marker?]

  ERROR: 4 novel marker(s) and no --kronos_additional_markers.
         Most are naming variants — try --kronos_marker_mapping first.
         Only add a CSV row for markers that are genuinely new biology.
```

---

## 7. Process specifications

### KRONOS2PREFLIGHT — once, `process_single`, **no model, no torch**

| | |
|---|---|
| Input | `.collect()` of all `[meta, tiff]`; optional `additional_markers.csv` |
| Reads | TIFF **metadata only** (OME-XML → ImageJ → numbered fallback) |
| Does | Reads `marker_metadata.csv` + `marker_utils.py` from `kronos_model_path`; applies `marker_match_key` to the union of channel names (after alias/isotope normalisation, section 2); classifies novel markers; validates CSV rows |
| Output | `novel_markers.json`, `preflight_report.txt` |
| Fails when | A novel marker has no CSV row, or a row leaves a required text column blank, or `kronos_additional_markers` is null and novel markers exist and `kronos_allow_novel_defaults` is false — each error naming the offending marker |

Required text columns (user-supplied; they feed BioLinkBERT and cannot be invented): `marker_name, compartment, family, family_desc, marker_full_name`. `mean`/`std` are left blank for us to fill.

**Two hard constraints, read off `modeling_kronos2.py:88-154` — both reject *before* any state change, so they are cheap to validate in PREFLIGHT:**

1. **Novel markers must reuse an existing `compartment` and `family`.** A new category has no row in the model's fixed embedding tables and raises `ValueError`. The allowed sets, extracted from the real `marker_metadata.csv`:

   - `compartment` (4): `cytoplasm`, `extracellular`, `membrane`, `nucleus`
   - `family` (28): `B_cell_marker`, `activation_marker`, `b_cell_marker`, `cell_cycle`, `chromatin_remodeler`, `cytoskeleton`, `dna_stain`, `ecm_matrix`, `endothelial_marker`, `epithelial_junction`, `epithelial_keratin`, `epithelial_marker`, `fibroblast_marker`, `histone_core`, `histone_mod`, `immune_checkpoint`, `metabolic_enzyme`, `metabolic_transporter`, `mhc_i`, `mhc_ii`, `myeloid_marker`, `neuron_marker`, `nk_cell_marker`, `pan_leukocyte`, `pericyte_marker`, `qc_mask`, `t_cell_marker`, `transcription_factor`

   (Note `B_cell_marker` and `b_cell_marker` both exist — a casing duplicate in their CSV. Match exactly; do not normalise.)

2. **Partial re-registration is unsupported.** If some markers in the CSV are already registered and others are not, it raises. `POOLSTATS` must therefore emit a CSV containing **only** the not-yet-registered novel markers, never the full panel. Registering an already-complete CSV is a safe no-op.

PREFLIGHT should validate both up front so the user learns at minute zero, not after the per-sample stats pass.

### KRONOS2PREPARE — per sample, `process_medium`, **no model**

| | |
|---|---|
| Input | `[meta, tiff, wc_mask]` + `novel_markers.json` |
| Does | For each novel marker present on this slide: take the region's pixels, scale by the dtype divisor to `[0,1]`, compute `(n, mean, var)` with `ddof=1` |
| Output | `<sample>_partials.json` — `{sample, region_key, scaling, stats: {marker: {n, mean, var}}}` |

Pure numpy, no GPU, no model. Mirrors `Kronos2Extractor.prepare_slide` (`src/coral/features/kronos2.py:192-243`).

### KRONOS2POOLSTATS — once, `process_single`, **no model**

Pools with the exact variance decomposition (`src/coral/markers/additional.py:315`), so the result equals computing over all pixels concatenated:

```
μ  = ( Σ nᵢ·μᵢ ) / ( Σ nᵢ )
s² = [ Σ (nᵢ−1)·sᵢ² + Σ nᵢ·(μᵢ−μ)² ] / ( Σ nᵢ − 1 )
```

Fills **only blank** `mean`/`std` cells — a user-supplied value is never overwritten, making the step idempotent. Writes via temp-file + `os.replace`.

### KRONOS2EMBEDDINGS — per sample, `process_gpu`, model on GPU

| | |
|---|---|
| Input | `[meta, tiff, wc_mask, geojson]` + filled `additional_markers.csv` (optional) |
| Does | `register_additional_markers(csv)` if present → centroids → isolated patches → `model.preprocess(...)` → `model(t, markers)` → 768-d CLS |
| Output | merged GeoJSON (`*.geojson{,.gz}`) + `_marker_report.txt` — **no CSV** |

Core loop, per `src/coral/features/kronos2.py:245-381`:

```python
scaled = np.asarray(patches, dtype=np.float32)      # already / scaling_factor(dtype)
normed = model.preprocess(scaled, markers, preferred_dapi=nuclear_marker)
with torch.inference_mode():
    cls = model(torch.from_numpy(np.ascontiguousarray(normed)).to(dev), markers)
emb = cls.float().cpu().numpy().astype(np.float32)  # (n, 768)
```

Use the `process_gpu` label rather than `process_multi` + hand-rolled directives, and **fix the container options**: the current block hardcodes `containerOptions = '--nv'`, which breaks under docker. `CELLSAMSEGMENT` already has the correct conditional (`--nv` for apptainer/singularity, `--gpus all` for docker) — copy it.

---

## 8. Stats region — a genuine divergence from CORAL, to be documented

CORAL computes novel-marker stats over the **tissue mask**. sp_segment has no tissue-detection stage, so we cannot reproduce that region.

`kronos_stats_region` default **`'cells'`** — the union of the whole-cell mask, i.e. exactly the pixels that get embedded. `'image'` computes over the full frame. Neither equals CORAL's tissue-masked region, so **stats computed here will not reproduce CORAL's numbers**, and the choice must be recorded in the CSV provenance and the marker report.

This is the strongest practical argument for eventually adding an Otsu tissue stage (`findings.md` §6.3.1) — it would make the two pipelines' novel-marker stats directly comparable. Out of scope here.

---

## 9. Model access, container, reproducibility

**Container.** A new `environment.yml` is required — do **not** try to share the v1 one, which pins cu128 and installs `kronos` from git:

These are the model repo's **own** `requirements.txt`, not inferred from CORAL:

```
--extra-index-url https://download.pytorch.org/whl/cu124
torch==2.6.0
xformers==0.0.29.post3 ; sys_platform == "linux"
transformers==4.56.0
timm==1.0.19
omegaconf>=2.3.0
huggingface_hub>=0.28.1
numpy
pandas
tifffile
safetensors
```

Note `xformers` is **linux-only** in their pin — relevant if anyone attempts a local Windows/macOS run. Add `setuptools>=70` (xformers' Triton probe wants it); `pandas` is already required by `marker_utils.py`.

**Weights.** Gated repo, access approved. Prefer a pre-downloaded local directory over a live fetch — `AutoModel.from_pretrained` accepts a path, and this is the only thing that works on egress-free compute nodes:

```bash
huggingface-cli download MahmoodLab/KRONOS2 --local-dir /path/to/KRONOS2
nextflow run . --enable_kronos --kronos_model_path /path/to/KRONOS2
```

Note `trust_remote_code=True` means the repo ships **custom Python** that must also be present offline; `huggingface-cli download` fetches it, but verify the module loads with `HF_HUB_OFFLINE=1` before trusting it on a cluster. `secret 'HF_TOKEN'` is the fallback for nodes with egress.

**Import noise.** KRONOS2's vendored DINOv2 + xFormers emit a caught Triton-probe traceback and "xFormers is available" warnings at import that read as failures but are harmless. Silence them exactly as CORAL does (`src/coral/features/kronos2.py:170-186`): raise the `xformers` logger to `ERROR` and filter that one warning. Skipping this guarantees bug reports.

**Batch size.** v2 runs fp32 with no autocast; on GPU cuBLAS picks batch-dependent kernels, so results shift ~5e-5 (batch 8 vs 16) to ~2e-4 (batch 4). The published gold standard is **batch 16**, hence the default change from 32. Document that changing it changes the numbers slightly.

---

## 10. Output — GeoJSON only

**Decision: the merged GeoJSON is the sole embedding artefact. The CSV is dropped.**

v1 emitted both a `*_kronos_embeddings.csv` and a GeoJSON carrying the same vectors as `properties.measurements` keys — the same data twice. One artefact, QuPath-native, joined by construction, is simpler and removes the CSV↔GeoJSON consistency question entirely.

**Concrete removals:**

- `--output` (CSV path) argument and the `save_embeddings()` writer
- the `embeddings` output channel and its `emit:` forwarding through all six subworkflows
- the `${params.outdir}/kronosembeddings` publish target for `*.csv`

`--output-geojson` becomes the primary output path; the marker report derives its name from it rather than from the CSV path (v1 did `args.output.replace(".csv", "_marker_report.txt")`, which must be rewritten).

**Size, for the record.** 768 floats per cell as JSON keys is roughly 29 KB/cell uncompressed — ~2.9 GB at 100k cells, ~29 GB at 1M, before `gzip_geojson` (typically 3-4× reduction). This is 2× the v1 384-d payload. If it becomes unworkable on a large slide, the fallback is to inline a reduced set (first *k* PCs) and write full vectors separately — deliberately **not** built now, since it adds back the second-artefact problem this decision removes.

**Keep from v1**, both genuinely good and easy to lose in a rewrite:

- clearing stale `kronos_*` keys from `measurements` before writing, so re-runs are idempotent rather than accumulating columns (`bin/kronos_embeddings.py:658-668`)
- the RSS-annotated `_log` progress output, which is more useful on long HPC jobs than a tqdm bar

**Also drop with the CSV:** v1's GeoJSON-derived full-resolution `uint32` mask and its `{prefix}_geojson_mask.tif` side-output. Under cell isolation we need per-cell footprints, but those should be rasterised **per patch** from the polygon, not once for the whole slide — a 30k² slide otherwise allocates and writes 3.6 GB for something the polygons already describe (`findings.md` §6.2.2).

---

## 11. Marker vocabulary reference — "is my marker novel?"

Users need to answer this **before** running the pipeline, not by reading a failure. Nothing available today provides it: CORAL vendors only its own `registry_v1.csv` and the KRONOS1 177-marker table (`_kronos_marker_meta.py`); the KRONOS2 vocabulary exists **only inside the model**, reachable via `model._marker_stats` / `model._marker_index`.

So we generate and commit it. Three parts:

### 11.1 Generator — `bin/kronos2_list_markers.py`

Loads the model and dumps its vocabulary. Also the `--check` mode, which is the question users actually have:

```bash
# dump the full vocabulary
kronos2_list_markers.py --model-path $KRONOS2_MODEL_DIR \
    --out-csv assets/kronos2_marker_vocabulary.csv \
    --out-md docs/kronos2_markers.md

# answer "are any of MY markers novel?" — from a panel list or straight off a TIFF
kronos2_list_markers.py --model-path $KRONOS2_MODEL_DIR --check-tiff sample.ome.tif
kronos2_list_markers.py --model-path $KRONOS2_MODEL_DIR --check CD8,FoxA1,HLA-DR
```

`--check` must classify using the model's **own** `clean_marker_name` / `marker_match_key` / `build_marker_index`, never a local reimplementation — same rule as `PREFLIGHT` (§7). A user-facing answer that disagrees with what the pipeline does later is worse than no answer.

### 11.2 Nothing is committed — generate locally (decided 2026-07-28)

`MahmoodLab/KRONOS2` is **CC-BY-NC-ND-4.0**. The *NoDerivatives* clause makes committing a vocabulary table derived from their `marker_metadata.csv` into this **public** repo a redistribution question we do not need to have. So:

- **No vendored artefacts.** `assets/kronos2_marker_vocabulary.csv` and `docs/kronos2_markers.md` are **not** created.
- The generator reads `marker_metadata.csv` from the user's own `kronos_model_path` and writes wherever they point it — a private location, an internal wiki, or the run's `--outdir`.
- `.gitignore` blocks `kronos2_marker_vocabulary.csv`, `kronos2_markers.md`, `marker_metadata.csv`, and `KRONOS2/` so a generated copy can never be committed by accident.
- Public docs **link** to the HF model page rather than mirroring its contents.

**This deletes the drift problem rather than managing it.** An earlier draft had a whole subsection on revision hashes and preflight drift warnings, needed only because a committed copy could go stale. Generating from the user's own model directory means the list is *by construction* the version they are running. That subsection is gone.

**What public docs may still state**, because it is an API constraint rather than a copy of their dataset — a user physically cannot author a valid `additional_markers.csv` without it:

- the headline counts (288 markers, 268 pretraining / 20 stats-only)
- the 4 legal `compartment` values and 28 legal `family` values (section 7)
- a link to https://huggingface.co/MahmoodLab/KRONOS2

This is the same reasoning as documenting an enum's permitted values. The 288-row stats table itself stays unmirrored.

**Still fine to commit:** `docs/examples/kronos_marker_mapping.json` and `docs/examples/additional_markers.csv`. These are *user-authored config*, not copies of KRONOS2 data — between them they name five markers, in the same way any config example names the API values it uses.

---

## 12. Testing

| Test | Covers |
|---|---|
| `nf-test` stub for each of the 4 modules | Wiring without weights — runs in CI |
| `conf/test_mesmer_kronos2.config` | Full path; needs `KRONOS2_MODEL_DIR`, so gated to a manual/nightly workflow |
| Synthetic panel `DAPI` + `CD45` | The **zero-novel-markers early exit** — the common path |
| Synthetic panel + one nonsense channel name | The novel path: preflight → prepare → pool → register |
| `tests/data/comet/additional_markers.csv` | Fixture with text columns filled, `mean`/`std` blank |
| Unit test on the pooling function | Pooled `(mean, std)` equals `np.concatenate` of the inputs — CORAL's own doctest (`additional.py:337-348`) transfers directly |

The last one is worth having as a plain pytest: it is pure arithmetic, needs no weights, and is the piece most likely to be silently wrong.

Also add the missing `tests/data/comet/kronos_marker_metadata.csv` situation to the CHANGELOG as a fixed bug — the v1 test profile could never have run.

---

## 13. Suggested PR sequence

Four reviewable PRs rather than one large one:

| PR | Contents | Reviewable without weights? |
|---|---|---|
| 1 | Hoist KRONOS out of the 3 subworkflows to top level; **no behaviour change**, still v1 | Yes — pure refactor, existing tests must pass unchanged |
| 2 | Delete v1; add `kronos2_common.py` + `KRONOS2EMBEDDINGS`; dtype divisor, nuclear hint, cell isolation; params + schema + docs | Stub tests yes; numbers need weights |
| 3 | `KRONOS2PREFLIGHT` + novel-marker detection and fail-fast | Yes (stub + unit) |
| 4 | `KRONOS2PREPARE` + `KRONOS2POOLSTATS` + `register_additional_markers` | Pooling unit test yes; end-to-end needs weights |

PR 1 being a pure refactor is the point — it makes the v2 diff in PR 2 readable, and it can be validated against the existing test suite before any model risk is introduced.

---

## 14. Open items

1. **`SEGMENTATIONREPORT` input** — confirm the report does not read embedding fields before accepting that it consumes the pre-KRONOS GeoJSON (section 3).
2. **Real panel novelty** — once the model is downloaded, run `KRONOS2PREFLIGHT` against a real COMET panel to find out how many markers are actually novel. If the answer is zero, PRs 3-4 are still correct but their urgency drops sharply.
3. **`kronos_stats_region` default** — `'cells'` assumes the whole-cell mask is a good stand-in for tissue. Worth sanity-checking on one real slide against `'image'`.
4. **Isolation and smoothing interact** — when `smooth_masks` is on, the isolated footprint comes from the smoothed polygon, not the raw segmentation. Probably desirable (it matches what `cellmeasurement` reports) but should be a conscious choice.
5. **`ro-crate-metadata.json`** references KRONOS; regenerate rather than hand-edit.
6. **Patch size 256 vs 64.** `config.json` declares `patch_size: 256` and `global_crops_size: 224`, and `modeling_kronos2.py` documents `forward` as taking `(B, n_markers, 256, 256)`. Our plan defaults `kronos_patch_size = 64` for cell patches. CORAL's tutorial 3 does run KRONOS2 on 64 px cell patches successfully (the ViT interpolates position embeddings, `interpolate_offset: 0.1`), so this works — but it is **off-nominal input** and should be validated against a 256 px run before being trusted for phenotyping. Consider making 256 the documented default for grid mode if that is ever added.
7. **`num_markers: 512`** in the backbone config is the hard ceiling on panel size. Not a concern for typical fluorescence panels, but worth a clear error rather than a deep crash if ever exceeded.
8. **`from_pretrained` needs the `dinov2/` package directory** from the snapshot (it does `sys.path.insert(0, d)` then imports `dinov2.models.inference`). An offline model dir must therefore be the **full snapshot minus `demo_image/`** — not a hand-picked file subset. Their own loader already uses `ignore_patterns=["demo_image/*"]`; mirror that.
9. ~~**Licensing.**~~ **Resolved 2026-07-28** — nothing derived from KRONOS2 is committed. The vocabulary is generated locally from the user's own model directory and kept private; `.gitignore` enforces it. See section 11.2.
