# nf-test suite status

Why this file exists: CI runs only `nextflow run`, so the snapshot suite has
never executed automatically. It had rotted accordingly. This records what was
repaired, what is still broken, and — for the open items — exactly what is known
and what is not, so nobody has to re-derive it.

## Repaired

The broken test and the `backgroundsubtract` snapshot are committed. The
regenerated snapshots land separately, once the phase-2 verification described
below has passed — a regenerated snapshot that has only been checked by the run
that wrote it is not evidence of anything.

| Test                                                  | Was                                                                        |
| ----------------------------------------------------- | -------------------------------------------------------------------------- |
| `sopa_patchifyimage - zarr`                           | Dead since `b015b6f`: setup called `SOPA_CONVERT` with 1 input, it takes 2 |
| `sopa_convert - tiff` / `- stub`                      | Stale `versions.yml` md5                                                   |
| `sopa_resolvecellpose - stub`                         | Stale `versions.yml` md5                                                   |
| `sopa_segmentationcellpose - stub`                    | Stale `versions.yml` md5                                                   |
| `sopa_patchifyimage - zarr - stub`                    | Stale since PR #35 fixed the empty patch-count file                        |
| `sopa_segment_compartment - whole-cell` / `- nuclear` | Snapshot predates `patch_width_pixel = 250`                                |
| `backgroundsubtract - tiff`                           | Snapshotted an OME-TIFF whose UUID changes every run                       |

The `versions.yml` md5s had flip-flopped across commits
(`e1cb9bf9` → `136a272d` → `451efc07` → `136a272d`). `451efc07` is the md5 of a
`versions.yml` where the sopa version rendered **empty**, so at least one past
regeneration ran somewhere `sopa --version` produced nothing, and was committed.
They are deterministic within an environment — three consecutive runs agree — so
this is staleness, not flakiness. Regenerate under `-profile test,docker`.

The compartment snapshots recorded 1 patch. `b015b6f` created them under the
3000 px default; `18c4b81` then set `patch_width_pixel = 250` in the test
profile, which yields 9 patches. Everything downstream — parquets, shapes, the
mask — differed as a consequence.

**Regenerating is a two-phase job.** `--update-snapshot` compares against what it
just wrote, so it always agrees with itself; only a second run _without_ it
proves reproducibility. That omission is what let the backsub-UUID defect
through originally.

## Unresolved

### 1. Is CellSAM deterministic on CPU?

**Status: unresolved.** Do not assume either answer.

The first attempt was **invalid**. `conf/base.config` sets `accelerator = 1` for
`withLabel:process_gpu`, and Nextflow's Docker executor expands that into
`--gpus all`. `CELLSAMSEGMENT` carries that label, so `-profile test,docker`
runs on the **GPU** despite no `--gpus` appearing anywhere in the configs.
Confirmed in the task's `.command.run`:

```
docker run -i --cpu-shares 4096 ... --gpus all --name $NXF_BOXID ...
```

What that (GPU) run does establish: `cellsam_segment - tiff` was byte-identical
across two independent runs and matched the committed snapshot —
`test_nuclear.tiff:md5,219b3fe8…`, `test_whole-cell.tiff:md5,7915618d…`,
`test.geojson.gz:md5,2125a4e0…`.

What it does **not** establish: anything about CPU, which is what a GPU-less CI
runner would use.

To resume, force CPU and assert it took effect:

```groovy
// tests/nextflow.config — for the check only, do not commit
env.CUDA_VISIBLE_DEVICES = ''
process { withLabel: process_gpu { accelerator = 0 } }
```

then run `cellsam_segment - tiff` twice and compare snapshot content directly,
rather than trusting a pass/fail against a committed value. Budget several
hours: see below for why.

### 2. CellSAM model weights are fetched per task

Not new, and not a regression — but it is why the check above is slow, and it
has almost no headroom.

Every `CELLSAMSEGMENT` task downloads the CellSAM v1.2 weights from DeepCell.
One task was measured pulling the full **1.70 GB in 41 m 38 s**. Another hit the
**1 h task limit at 78%** when throughput fell to ~50 kB/s, and failed.

`HOME=/root` inside the container and only the work directory is bind-mounted,
so `$HOME/.deepcell` is container-local and discarded per task. Nothing is
cached between tasks. `params.deepcell_cache_dir` exists to redirect it and
defaults to `null`; using it also requires bind-mounting that path into the
container.

This normally succeeds and so has gone unnoticed. It is a poor fit for CI.

### 3. `segmentationreport - basic functionality`

The output HTML md5 differs from the committed snapshot. **Not diagnosed.**

Hypothesis only: the report is produced by `quarto render --to html`, and Quarto
normally stamps output with a render date, which would make the md5
unreproducible by construction — the same category as the backsub UUID, and
equally not ours to fix. If that is confirmed, the fix is the same: assert the
file is produced, drop its md5 from the snapshot. **Verify before acting on
this** — run the test twice and compare.

### 4. `cellsamsegment - nuclear` / `- whole-cell` (module level)

Mask md5s differ from the committed snapshot. **Not diagnosed, and deliberately
not regenerated.**

These were left alone on purpose: regenerating them before (1) is settled risks
baking in a GPU-produced mask that a GPU-less CI runner could never reproduce.
Settle the device question first.

Note the module-level tests fail while the subworkflow-level `cellsam_segment`
tests pass. They use different fixtures, so stale snapshots are the more likely
explanation — but that is a guess, not a finding.

## Before enabling nf-test in CI

- `process_gpu` tasks request a GPU via `accelerator = 1`. GitHub's
  `ubuntu-latest` has none. Establish what Nextflow does there before running
  any test that touches `CELLSAMSEGMENT` or `MESMERSEGMENT` in CI.
- Cellpose weights are staged once per run by `CELLPOSEMODEL`, and
  `--cellpose_models_dir` skips the download entirely. Cache that directory in
  CI keyed on the sopa container tag.
- CellSAM has no equivalent yet; see (2).
