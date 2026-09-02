"""Unit tests for well/tile discovery in ``bin/aviti_discover_tiles.py``.

These exercise pure filesystem/JSON logic (RunParameters.json parsing, channel
mode auto-detection/enforcement, manifest row construction) with tiny
synthetic run directories -- no real AVITI data, no image content, no GPU.
"""

import json
from pathlib import Path

import pytest

from aviti_discover_tiles import (
    ChannelMode,
    discover_manifest_rows,
    select_wells,
)


def _make_run(
    tmp_path: Path,
    wells: dict,
    batch: str = "CP01",
    include_actin: bool = True,
) -> Path:
    """Build a minimal synthetic AVITI run directory.

    ``wells`` maps WellLocation -> list of (tile_name, x_mm, y_mm). Every tile
    gets empty Nucleus/Cell-Membrane files (and Actin, unless
    ``include_actin`` is False) -- discovery only checks for file existence,
    not image content.
    """
    run_dir = tmp_path / "run"
    projection_dir = run_dir / "Projection"
    projection_dir.mkdir(parents=True)

    run_parameters = {
        "WellLayout": "TwelveWellStandard",
        "Wells": [
            {
                "WellLocation": well,
                "Tiles": [
                    {"Name": name, "XMillimeters": x, "YMillimeters": y}
                    for name, x, y in tiles
                ],
            }
            for well, tiles in wells.items()
        ],
    }
    (run_dir / "RunParameters.json").write_text(json.dumps(run_parameters))

    for well, tiles in wells.items():
        well_dir = projection_dir / f"Well{well}"
        well_dir.mkdir(parents=True)
        for name, _x, _y in tiles:
            (well_dir / f"{batch}_{name}_Nucleus.tif").touch()
            (well_dir / f"{batch}_{name}_Cell-Membrane.tif").touch()
            if include_actin:
                (well_dir / f"{batch}_{name}_Actin.tif").touch()

    return run_dir


# ----------------------------------------------------------------------------
# well selection
# ----------------------------------------------------------------------------


def test_select_wells_returns_every_well_when_no_filter():
    run_parameters = {"Wells": [{"WellLocation": "A1"}, {"WellLocation": "A2"}]}
    wells = select_wells(run_parameters, None)
    assert [w["WellLocation"] for w in wells] == ["A1", "A2"]


def test_select_wells_filters_to_requested_wells_in_requested_order():
    run_parameters = {"Wells": [{"WellLocation": "A1"}, {"WellLocation": "A2"}, {"WellLocation": "B1"}]}
    wells = select_wells(run_parameters, ["B1", "A1"])
    assert [w["WellLocation"] for w in wells] == ["B1", "A1"]


def test_select_wells_fails_loudly_on_an_unknown_well():
    run_parameters = {"Wells": [{"WellLocation": "A1"}]}
    with pytest.raises(ValueError, match="not present"):
        select_wells(run_parameters, ["A1", "Z9"])


def test_select_wells_fails_on_empty_wells_list():
    with pytest.raises(ValueError, match="no 'Wells'"):
        select_wells({"Wells": []}, None)


# ----------------------------------------------------------------------------
# manifest discovery / channel mode
# ----------------------------------------------------------------------------


def test_discovers_one_row_per_tile_across_wells(tmp_path):
    run_dir = _make_run(
        tmp_path,
        {
            "A1": [("L1R01C01S1", 0.9, -0.08), ("L1R01C02S1", 2.663, -0.08)],
            "A2": [("L1R01C01S1", 0.9, -0.08)],
        },
    )
    rows = discover_manifest_rows(run_dir, None, "CP01", ChannelMode.AUTO)
    assert len(rows) == 3
    assert {r["well"] for r in rows} == {"A1", "A2"}
    assert all(r["channel_mode"] == "3ch" for r in rows)
    assert all(r["actin_tif"] for r in rows)


def test_auto_detects_2_channel_mode_when_actin_absent(tmp_path):
    run_dir = _make_run(
        tmp_path, {"A1": [("L1R01C01S1", 0.9, -0.08)]}, include_actin=False
    )
    rows = discover_manifest_rows(run_dir, None, "CP01", ChannelMode.AUTO)
    assert rows[0]["channel_mode"] == "2ch"
    assert rows[0]["actin_tif"] == ""


def test_forcing_3ch_without_actin_fails_loudly(tmp_path):
    run_dir = _make_run(
        tmp_path, {"A1": [("L1R01C01S1", 0.9, -0.08)]}, include_actin=False
    )
    with pytest.raises(FileNotFoundError, match="Actin"):
        discover_manifest_rows(run_dir, None, "CP01", ChannelMode.THREE_CHANNEL)


def test_inconsistent_channel_mode_across_run_fails_loudly(tmp_path):
    run_dir = _make_run(
        tmp_path,
        {"A1": [("L1R01C01S1", 0.9, -0.08)], "A2": [("L1R01C01S1", 0.9, -0.08)]},
    )
    # Remove Actin for the second well's tile only -- a mixed-mode run.
    (run_dir / "Projection" / "WellA2" / "CP01_L1R01C01S1_Actin.tif").unlink()

    with pytest.raises(ValueError, match="Inconsistent channel mode"):
        discover_manifest_rows(run_dir, None, "CP01", ChannelMode.AUTO)


def test_missing_required_channel_file_fails_loudly(tmp_path):
    run_dir = _make_run(tmp_path, {"A1": [("L1R01C01S1", 0.9, -0.08)]})
    (run_dir / "Projection" / "WellA1" / "CP01_L1R01C01S1_Cell-Membrane.tif").unlink()

    with pytest.raises(FileNotFoundError, match="Missing required channel"):
        discover_manifest_rows(run_dir, None, "CP01", ChannelMode.AUTO)


def test_missing_run_parameters_json_fails_loudly(tmp_path):
    run_dir = tmp_path / "empty_run"
    run_dir.mkdir()
    with pytest.raises(FileNotFoundError, match="RunParameters.json"):
        discover_manifest_rows(run_dir, None, "CP01", ChannelMode.AUTO)


def test_missing_projection_directory_fails_loudly(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "RunParameters.json").write_text(json.dumps({"Wells": [{"WellLocation": "A1", "Tiles": []}]}))
    with pytest.raises(FileNotFoundError, match="Projection"):
        discover_manifest_rows(run_dir, None, "CP01", ChannelMode.AUTO)


def test_stage_coordinates_are_carried_through_to_manifest_rows(tmp_path):
    run_dir = _make_run(tmp_path, {"A1": [("L1R01C01S1", 0.9, -0.08)]})
    rows = discover_manifest_rows(run_dir, None, "CP01", ChannelMode.AUTO)
    assert rows[0]["x_mm"] == 0.9
    assert rows[0]["y_mm"] == -0.08
