# KRONOS2 migration — handoff (2026-07-28)

Working notes for resuming. Delete before merging; the durable plan is
`docs/kronos2-migration-plan.md`.

## Where things are

Branch `refactor/kronos-hoist-toplevel` → PR #35 (open, against `dev`).

```
5040fe3 test: exclude kronos_input from snapshots (non-reproducible md5)
eb1cf55 fix: mark KRONOS2 bin scripts executable
713893e feat: replace KRONOS1 with KRONOS2 per-cell embeddings
79c9b5c refactor: invoke KRONOS once at top level instead of per segmenter
```

Both PR 1 (the top-level hoist) and PR 2 (the v1→v2 replacement) live on this
one branch, at your request.

## What is validated

- **Full pipeline, `-profile test,docker`** — exit 0, backsub → cellpose ×18 →
  resolve → parquet2tiff → cellmeasurement.
- **KRONOS2 end to end** — 119 cells → `(119, 768)` in the published GeoJSON.
  No NaN, no all-zero rows, no near-identical pairs; the 210 existing
  cellmeasurement columns preserved. Values match a standalone run of the
  script exactly, so the Nextflow path introduces no drift.
- **Cell isolation is real** — isolated vs raw-box embeddings differ
  (L2 mean 9.37, cosine mean 0.973) on identical cells and model.
- **Marker resolution** — the synthetic panel (DAPI, CD45) resolves with no
  mapping; a real 49-marker CODEX panel matched 45/49 against the true
  288-marker vocabulary.

## The one thing left to finish

**Snapshots.** They are regenerating in WSL at `/root/sp_segment` and are NOT
yet copied into the repo. To finish:

```bash
wsl -d Ubuntu-24.04 -u root
cd /root/sp_segment
grep -E "PASSED|FAILED|SUCCESS" nfregen.log      # confirm the regen passed
```

Then copy the six `.snap` files back and — this is the part that matters —
**verify with a second pass that does not use `--update-snapshot`**:

```bash
cp tests/nextflow.config /root/tests_nf.bak
echo "docker.runOptions = '-u \$(id -u):\$(id -g) -v /opt/cellpose_models:/tmp/cellpose_models'" >> tests/nextflow.config
nf-test test subworkflows/local/*/tests/main.nf.test --profile test,docker
cp /root/tests_nf.bak tests/nextflow.config       # do NOT commit that line
```

`--update-snapshot` compares against what it just wrote, so it always agrees
with itself. Only an independent run proves a snapshot is reproducible. That
omission is exactly what let the backsub-UUID defect through the first time.

## Blocked on you

1. **The container is a placeholder.** `modules/local/kronos2embeddings/main.nf`
   points at `community.wave.seqera.io/library/kronos2embeddings:PLACEHOLDER`.
   Build it via Seqera Wave from the module's `environment.yml` (as KRONOS1's
   was) and pin the digest. Two gotchas found the hard way:
   - The reference must be **fully qualified**, or `docker.registry = 'quay.io'`
     in `nextflow.config` silently rewrites it and you get a 401 from quay.
   - A **pip**-based build additionally needs `libexpat1` (rasterio) and
     `procps` (Nextflow's `ps` metrics). The conda `environment.yml` is fine
     as-is because conda-forge rasterio bundles its own libs.

2. **Rotate the DeepCell token** — it was pasted in plaintext during this
   session. It is registered as a Nextflow secret in WSL
   (`nextflow secrets list`).

3. **Two issues worth filing** (evidence below, neither belongs in this PR).

## Local environment built today

| Component | Detail |
|---|---|
| WSL2 | Ubuntu 24.04, systemd PID 1, `/dev/dxg` present |
| Java | OpenJDK 21 (Nextflow supports 17–21; the Windows JDK 25 is out of range) |
| Nextflow | 24.04.2, pinned to the CI matrix |
| nf-test | 0.9.5 |
| Docker | 29.1.3 + nvidia-container-toolkit, GPU verified in-container |
| KRONOS2 weights | `/opt/KRONOS2` (421 MB, full snapshot minus `demo_image/`) |
| KRONOS2 venv | `/opt/k2env` (torch 2.6.0+cu124) |
| Test image | `kronos2embeddings:local2` (9.45 GB, local only) |
| Cellpose model | `/opt/cellpose_models` (1.2 GB `cpsam`) |
| Repo clone | `/root/sp_segment` — run from here, **not** `/mnt/c` |

Two things that will bite if forgotten:

- Never rsync the repo from `/mnt/c`. Git's `autocrlf` gives Windows checkouts
  CRLF, and a `\r` in a shebang fails with exit 127. Use a git checkout.
- WSL `/tmp` is tmpfs and is cleared when the VM idles out. Write outputs to
  `/root/...`.

## Pre-existing bugs found (not caused by this work)

1. **Cellpose weights are unpinned.** The container is pinned to
   `sopa:2.1.11-cellpose`, but the 1.2 GB `cpsam` checkpoint is fetched from
   `cellpose.org` at runtime by each of ~18 tasks. This caused a CI failure on
   this PR (`HTTP 429`) and a local `PytorchStreamReader failed reading zip
   archive` from a truncated download. Caching it cut the `sopa_segment` test
   from 1548 s to 239 s. It also makes snapshots drift over time: unmodified
   `dev` fails its own committed snapshot in this environment.

2. **`cellpose_model_type` is a no-op.** The container ships cellpose 4.0.8,
   which logs `model_type argument is not used in v4.0.1+. Ignoring this
   argument` and always uses `cpsam`. The `cyto3` default has been silently
   discarded on every run, so past results attributed to cyto3 were cpsam.

3. **Stub runs were broken for the whole Cellpose path** — fixed here, since it
   blocked validating this work. `SOPA_PATCHIFYIMAGE`'s stub wrote an empty
   patch-count file that the caller parses with `.toInteger()`, so every stub
   run died with `For input string: ""`.

4. **`--nv` was passed unconditionally** to the KRONOS process, which is invalid
   under Docker — fixed here; the replacement selects by container engine.

## Next after this PR

Plan §7 stages 3 and 4: `KRONOS2PREFLIGHT` (novel-marker detection, needs no
model or torch — `marker_utils.py` is pure pandas and the vocabulary is a plain
CSV), then `KRONOS2PREPARE` / `KRONOS2POOLSTATS` for the cohort-pooled
novel-marker statistics.
