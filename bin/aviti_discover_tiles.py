#!/usr/bin/env python
'''
Module      : aviti_discover_tiles
Description : Discovers wells and tiles for an AVITI24 Teton/Teton Atlas
              cytoprofiling run and emits a flat manifest (one row per tile)
              that Nextflow can fan out on for parallel per-tile segmentation
              and later group back by well for stitching.

              Reads well/tile layout and stage coordinates from the run's
              RunParameters.json, and locates the raw per-channel projection
              images written by the instrument under:

                  <run_dir>/Projection/Well<WellLocation>/<batch>_<tile>_<Channel>.tif

              where <Channel> is one of Nucleus, Cell-Membrane, or (in
              3-channel cell paint mode) Actin. This mirrors the AVITI Cyto
              run output layout documented at
              https://docs.elembio.io/docs/elembio-cloud/runs/cyto-run-output/
              -- an independent implementation, not derived from Elembio's
              own (BSD-licensed) analysis notebook.
Copyright   : (c) WEHI SODA Hub, 2026
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import csv
import json
import sys
from enum import Enum
from pathlib import Path
from typing import Annotated, List, Optional

import typer

app = typer.Typer(add_completion=False)


class ChannelMode(str, Enum):
    AUTO = "auto"
    TWO_CHANNEL = "2ch"
    THREE_CHANNEL = "3ch"


NUCLEUS_SUFFIX = "Nucleus"
MEMBRANE_SUFFIX = "Cell-Membrane"
ACTIN_SUFFIX = "Actin"


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def load_run_parameters(run_dir: Path) -> dict:
    '''
    Load and return the parsed RunParameters.json for an AVITI run directory.
    '''
    run_parameters_path = run_dir / "RunParameters.json"
    if not run_parameters_path.is_file():
        raise FileNotFoundError(
            f"RunParameters.json not found under {run_dir}. Expected an AVITI "
            "run directory containing RunParameters.json and a Projection/ "
            "subdirectory."
        )
    with open(run_parameters_path) as fh:
        return json.load(fh)


def select_wells(run_parameters: dict, wells_filter: Optional[List[str]]) -> List[dict]:
    '''
    Return the list of well records to process, in RunParameters.json order.

    ``wells_filter`` restricts processing to specific WellLocation values (e.g.
    ["A1", "A2"]). Filtering fails loudly on a well that does not exist in the
    run, rather than silently processing fewer wells than requested -- a typo
    in --wells would otherwise look like an empty/short run.
    '''
    all_wells = run_parameters.get("Wells")
    if not all_wells:
        raise ValueError("RunParameters.json contains no 'Wells' entries.")

    if not wells_filter:
        return all_wells

    by_location = {well["WellLocation"]: well for well in all_wells}
    missing = [w for w in wells_filter if w not in by_location]
    if missing:
        raise ValueError(
            f"Requested well(s) not present in RunParameters.json: {missing}. "
            f"Available wells: {sorted(by_location)}"
        )
    return [by_location[w] for w in wells_filter]


def tile_channel_path(well_dir: Path, batch: str, tile_name: str, channel: str) -> Path:
    return well_dir / f"{batch}_{tile_name}_{channel}.tif"


def discover_manifest_rows(
    run_dir: Path,
    wells_filter: Optional[List[str]],
    cellpaint_batch: str,
    channel_mode: ChannelMode,
) -> List[dict]:
    '''
    Build one manifest row per (well, tile), validating that the raw
    per-channel projection TIFFs referenced by RunParameters.json actually
    exist on disk.

    channel_mode == AUTO is resolved once from the first tile encountered and
    then enforced for every subsequent tile: a run that silently mixed 2- and
    3-channel tiles would stitch and segment inconsistently, so a mismatch is
    an error rather than a per-tile fallback.
    '''
    run_parameters = load_run_parameters(run_dir)
    projection_dir = run_dir / "Projection"
    if not projection_dir.is_dir():
        raise FileNotFoundError(f"Projection directory not found: {projection_dir}")

    wells = select_wells(run_parameters, wells_filter)

    rows = []
    resolved_mode: Optional[ChannelMode] = (
        channel_mode if channel_mode != ChannelMode.AUTO else None
    )

    for well in wells:
        well_location = well["WellLocation"]
        well_dir = projection_dir / f"Well{well_location}"
        if not well_dir.is_dir():
            raise FileNotFoundError(
                f"Well directory not found for WellLocation {well_location}: {well_dir}"
            )

        tiles = well.get("Tiles") or []
        if not tiles:
            raise ValueError(f"Well {well_location} has no Tiles in RunParameters.json")

        for tile in tiles:
            tile_name = tile["Name"]
            nucleus = tile_channel_path(well_dir, cellpaint_batch, tile_name, NUCLEUS_SUFFIX)
            membrane = tile_channel_path(well_dir, cellpaint_batch, tile_name, MEMBRANE_SUFFIX)
            actin = tile_channel_path(well_dir, cellpaint_batch, tile_name, ACTIN_SUFFIX)

            missing = [p for p in (nucleus, membrane) if not p.is_file()]
            if missing:
                raise FileNotFoundError(
                    f"Missing required channel file(s) for well {well_location}, "
                    f"tile {tile_name}: {missing}"
                )

            actin_present = actin.is_file()
            this_tile_mode = ChannelMode.THREE_CHANNEL if actin_present else ChannelMode.TWO_CHANNEL

            if resolved_mode is None:
                resolved_mode = this_tile_mode
                log(
                    f"Auto-detected AVITI channel mode: {resolved_mode.value} "
                    f"(from well {well_location}, tile {tile_name})"
                )
            elif channel_mode == ChannelMode.AUTO and this_tile_mode != resolved_mode:
                raise ValueError(
                    f"Inconsistent channel mode across run: well {well_location}, "
                    f"tile {tile_name} looks like {this_tile_mode.value}, but an "
                    f"earlier tile resolved to {resolved_mode.value}. Pass "
                    "--channel-mode explicitly to override auto-detection."
                )
            elif channel_mode == ChannelMode.THREE_CHANNEL and not actin_present:
                raise FileNotFoundError(
                    f"--channel-mode 3ch requires an Actin file for well "
                    f"{well_location}, tile {tile_name}, but none was found: {actin}"
                )

            use_actin = resolved_mode == ChannelMode.THREE_CHANNEL and actin_present

            rows.append({
                "well": well_location,
                "tile": tile_name,
                "x_mm": tile["XMillimeters"],
                "y_mm": tile["YMillimeters"],
                "nucleus_tif": str(nucleus),
                "membrane_tif": str(membrane),
                "actin_tif": str(actin) if use_actin else "",
                "channel_mode": resolved_mode.value,
            })

    if not rows:
        raise ValueError("No tiles discovered -- check --wells and the run directory layout.")

    return rows


def write_manifest(rows: List[dict], output: Path) -> None:
    fieldnames = [
        "well", "tile", "x_mm", "y_mm",
        "nucleus_tif", "membrane_tif", "actin_tif", "channel_mode",
    ]
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


@app.command()
def main(
    run_dir: Annotated[Path, typer.Argument(
        exists=True, file_okay=False, dir_okay=True,
        help="Path to the AVITI run directory (containing RunParameters.json and Projection/)."
    )],
    output: Annotated[Path, typer.Option(help="Path to write the tile manifest CSV.")],
    wells: Annotated[Optional[str], typer.Option(
        help="Comma-separated list of WellLocation values to restrict processing to "
             "(e.g. 'A1,A2'). Unset processes every well in the run."
    )] = None,
    cellpaint_batch: Annotated[str, typer.Option(
        help="Imaging batch prefix holding the segmentation-relevant channels "
             "(nucleus/membrane/actin), e.g. 'CP01'."
    )] = "CP01",
    channel_mode: Annotated[ChannelMode, typer.Option(
        help="Force 2-channel (nucleus + membrane) or 3-channel "
             "(+ actin) mode, or auto-detect from the presence of the Actin file."
    )] = ChannelMode.AUTO,
):
    '''
    Discover AVITI wells/tiles and write a flat per-tile manifest CSV.
    '''
    if not run_dir.is_dir():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")

    wells_filter = None
    if wells:
        wells_filter = [w.strip() for w in wells.split(",") if w.strip()]

    rows = discover_manifest_rows(
        run_dir=run_dir,
        wells_filter=wells_filter,
        cellpaint_batch=cellpaint_batch,
        channel_mode=channel_mode,
    )
    write_manifest(rows, output)

    n_wells = len({row["well"] for row in rows})
    log(f"Discovered {len(rows)} tile(s) across {n_wells} well(s); manifest written to {output}")


if __name__ == "__main__":
    app()
