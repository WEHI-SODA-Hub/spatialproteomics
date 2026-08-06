"""Unit tests for the channel-handling half of ``bin/kronos2_common.py``.

These cover the code that decides *which pixels get which marker name* -- channel
naming, exclusion, the windowed reader, and the marker report. That is worth
pinning because its failure mode is silent: the run completes, the channel counts
look right, and the embeddings are wrong.

Pure python -- no torch, no GPU, no model weights.
"""

import numpy as np
import pytest
import tifffile

from kronos2_common import (
    ChannelReader,
    apply_marker_mapping,
    channel_layout,
    get_channel_names,
    select_channels,
    write_marker_report,
)

PANEL = ["DAPI", "CD8", "TCR", "GrB", "HLA-DR", "PD-1", "Pan-CK", "Autofluorescence"]


# ----------------------------------------------------------------------------
# fixtures: the two OME-TIFF layouts
# ----------------------------------------------------------------------------


def _plane(value, size=16):
    """A 2D plane whose every pixel is ``value``, so channels are separable."""
    return np.full((size, size), value, dtype=np.uint16)


@pytest.fixture
def stacked_tiff(tmp_path):
    """The usual layout: one CYX series holding every channel."""
    path = tmp_path / "stacked.ome.tif"
    data = np.stack([_plane(i + 1) for i in range(len(PANEL))])
    tifffile.imwrite(
        path, data, metadata={"axes": "CYX", "Channel": {"Name": PANEL}}
    )
    return path


@pytest.fixture
def split_tiff(tmp_path):
    """Each channel in its own 2D series -- the layout that broadcast DAPI.

    Written the way the inForm unmixing export writes it: a separate,
    non-contiguous series per channel.
    """
    path = tmp_path / "split.ome.tif"
    with tifffile.TiffWriter(path, ome=True) as writer:
        for i, name in enumerate(PANEL):
            writer.write(
                _plane(i + 1),
                contiguous=False,
                metadata={"axes": "YX", "Channel": {"Name": name}},
            )
    return path


# ----------------------------------------------------------------------------
# select_channels
# ----------------------------------------------------------------------------


def test_no_exclusion_keeps_every_channel_in_order():
    keep, kept, dropped = select_channels(PANEL)
    assert keep == list(range(len(PANEL)))
    assert kept == PANEL
    assert dropped == []


def test_excluding_the_last_channel():
    keep, kept, dropped = select_channels(PANEL, ["Autofluorescence"])
    assert keep == [0, 1, 2, 3, 4, 5, 6]
    assert kept == PANEL[:-1]
    assert dropped == ["Autofluorescence"]


def test_excluding_a_middle_channel_preserves_order_and_indices():
    keep, kept, dropped = select_channels(PANEL, ["TCR"])
    assert keep == [0, 1, 3, 4, 5, 6, 7]
    assert kept == ["DAPI", "CD8", "GrB", "HLA-DR", "PD-1", "Pan-CK", "Autofluorescence"]
    assert dropped == ["TCR"]


def test_repeated_and_blank_entries_are_tolerated():
    keep, kept, dropped = select_channels(PANEL, ["  TCR  ", "TCR", "", "   "])
    assert kept == ["DAPI", "CD8", "GrB", "HLA-DR", "PD-1", "Pan-CK", "Autofluorescence"]
    assert dropped == ["TCR"]
    assert 2 not in keep


def test_duplicate_channel_names_are_all_dropped():
    names = ["DAPI", "Blank", "CD8", "Blank"]
    keep, kept, dropped = select_channels(names, ["Blank"])
    assert keep == [0, 2]
    assert kept == ["DAPI", "CD8"]
    assert dropped == ["Blank", "Blank"]


def test_unknown_name_fails_and_lists_the_available_channels():
    with pytest.raises(ValueError) as excinfo:
        select_channels(PANEL, ["Autofluorecence"])
    message = str(excinfo.value)
    assert "Autofluorecence" in message
    for name in PANEL:
        assert name in message


def test_matching_is_case_sensitive():
    with pytest.raises(ValueError):
        select_channels(PANEL, ["autofluorescence"])


def test_the_nuclear_marker_cannot_be_excluded():
    with pytest.raises(ValueError, match="preferred_dapi"):
        select_channels(PANEL, ["DAPI"], protect="DAPI")


def test_protect_does_not_block_other_channels():
    _, kept, dropped = select_channels(PANEL, ["Autofluorescence"], protect="DAPI")
    assert dropped == ["Autofluorescence"]
    assert "DAPI" in kept


def test_excluding_everything_fails():
    with pytest.raises(ValueError, match="nothing to embed"):
        select_channels(PANEL, list(PANEL))


def test_kept_names_stay_aligned_with_the_mapped_marker_list():
    """The invariant the model depends on: one marker name per kept channel."""
    keep, kept, _ = select_channels(PANEL, ["Autofluorescence"])
    markers, _ = apply_marker_mapping(kept, {"TCR": "TCR_B", "Pan-CK": "CYTOKERATIN"})
    assert len(markers) == len(keep)


# ----------------------------------------------------------------------------
# channel naming and layout
# ----------------------------------------------------------------------------


def test_names_and_layout_of_a_stacked_image(stacked_tiff):
    assert get_channel_names(stacked_tiff) == PANEL
    with tifffile.TiffFile(stacked_tiff) as tif:
        assert channel_layout(tif) == ("stacked", len(PANEL))


def test_names_and_layout_of_a_split_image(split_tiff):
    """The regression: 8 one-channel series must count as 8 channels, not 1."""
    with tifffile.TiffFile(split_tiff) as tif:
        assert channel_layout(tif) == ("split", len(PANEL))
    assert len(get_channel_names(split_tiff)) == len(PANEL)


def test_an_unnamed_channel_gets_a_synthesised_name(tmp_path):
    """Names must stay positional; dropping one would shift every later marker."""
    path = tmp_path / "unnamed.ome.tif"
    data = np.stack([_plane(i + 1) for i in range(3)])
    tifffile.imwrite(
        path, data, metadata={"axes": "CYX", "Channel": {"Name": ["DAPI", "", "CD8"]}}
    )
    names = get_channel_names(path)
    assert len(names) == 3
    assert names[0] == "DAPI"
    assert names[2] == "CD8"
    assert names[1]


# ----------------------------------------------------------------------------
# ChannelReader
# ----------------------------------------------------------------------------


def _window(reader):
    return reader.read_window(0, 4, 0, 4)


def test_stacked_reader_returns_the_requested_channels(stacked_tiff):
    reader = ChannelReader(stacked_tiff, [0, 2, 5])
    try:
        window = _window(reader)
        assert window.shape == (3, 4, 4)
        assert [int(window[i].flat[0]) for i in range(3)] == [1, 3, 6]
    finally:
        reader.close()


def test_split_reader_returns_distinct_channels(split_tiff):
    """Before the fix this returned one channel, broadcast over all of them."""
    reader = ChannelReader(split_tiff, list(range(len(PANEL))))
    try:
        assert reader.layout == "split"
        assert reader.n_channels == len(PANEL)
        window = _window(reader)
        assert window.shape == (len(PANEL), 4, 4)
        assert [int(window[i].flat[0]) for i in range(len(PANEL))] == list(
            range(1, len(PANEL) + 1)
        )
    finally:
        reader.close()


def test_split_reader_honours_a_channel_subset(split_tiff):
    keep, _, _ = select_channels(PANEL, ["Autofluorescence"])
    reader = ChannelReader(split_tiff, keep)
    try:
        window = _window(reader)
        assert window.shape == (7, 4, 4)
        # Channel 8 (value 8) is the excluded one and must not appear.
        assert 8 not in {int(window[i].flat[0]) for i in range(7)}
    finally:
        reader.close()


def test_reader_refuses_to_return_too_few_channels(split_tiff, monkeypatch):
    """The guard against silently broadcasting one channel across many markers."""
    reader = ChannelReader(split_tiff, [0, 1, 2])
    try:
        monkeypatch.setattr(
            reader, "_read_window", lambda *_: np.zeros((1, 4, 4), dtype=np.uint16)
        )
        with pytest.raises(ValueError, match="refusing to broadcast"):
            _window(reader)
    finally:
        reader.close()


# ----------------------------------------------------------------------------
# marker report
# ----------------------------------------------------------------------------


def test_report_without_exclusion_is_unchanged(tmp_path):
    """keep_indices=None must reproduce the pre-exclusion report exactly."""
    path = tmp_path / "report.txt"
    markers, applied = apply_marker_mapping(PANEL, {"TCR": "TCR_B"})
    write_marker_report(path, "img.ome.tif", PANEL, markers, applied)
    text = path.read_text()

    assert f"Channels: {len(PANEL)}\n" in text
    assert "excluded" not in text
    assert "  TCR -> TCR_B  (mapped)\n" in text
    assert "  DAPI -> DAPI\n" in text
    assert "\nMappings applied via --marker-mapping:\n  TCR -> TCR_B\n" in text


def test_report_records_excluded_channels(tmp_path):
    path = tmp_path / "report.txt"
    keep, kept, _ = select_channels(PANEL, ["Autofluorescence"])
    markers, applied = apply_marker_mapping(kept, {"TCR": "TCR_B"})
    write_marker_report(
        path, "img.ome.tif", PANEL, markers, applied, keep_indices=keep
    )
    text = path.read_text()

    assert "Channels: 8 (7 embedded, 1 not shown)" in text
    assert "  Autofluorescence -> (not shown to KRONOS2)\n" in text
    assert "still present in the image" in text
    assert "\n  Autofluorescence\n" in text
    assert "  TCR -> TCR_B  (mapped)\n" in text


def test_report_attributes_rows_by_index_not_name(tmp_path):
    """Duplicate channel names must not mis-attribute a row."""
    path = tmp_path / "report.txt"
    names = ["DAPI", "Blank", "CD8", "Blank"]
    keep, kept, _ = select_channels(names, ["Blank"])
    markers, applied = apply_marker_mapping(kept, {})
    write_marker_report(path, "img.ome.tif", names, markers, applied, keep_indices=keep)
    text = path.read_text()

    assert text.count("Blank -> (not shown to KRONOS2)") == 2
    assert "  DAPI -> DAPI\n" in text
    assert "  CD8 -> CD8\n" in text


def test_report_carries_extra_lines(tmp_path):
    path = tmp_path / "report.txt"
    write_marker_report(
        path, "img.ome.tif", PANEL, PANEL, [], ["Novel markers (default stats used):"]
    )
    assert "Novel markers (default stats used):" in path.read_text()
