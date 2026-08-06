#!/usr/bin/env python3
"""
Shared scaffold for KRONOS2 per-cell embedding extraction.

Everything that is not the model itself lives here: channel-name reading,
marker-name resolution, cell geometry, lazy patch extraction, and GeoJSON
output. A KRONOS2 extractor script only has to supply the two model-specific
steps -- ``model.preprocess`` (the marker-aware z-score) and ``model(...)``
(the forward pass) -- which keeps a future encoder from becoming a second
thousand-line copy of this file.

Design notes
------------
Scaling: patches are divided by a divisor derived from the *image dtype*
(uint8 -> 255, uint16 -> 65535, float -> 400). This mirrors KRONOS2's own
``sp_image._scaling_factor``; a hardcoded 65535 silently mis-scales any
non-uint16 input (e.g. a float32 background-subtracted image).

Marker names: KRONOS2 matches marker names *exactly* on a
separator-insensitive key, with no alias or fuzzy step (see the model's
``marker_utils.marker_match_key``). Names such as ``Cytokeritin`` or
``Collagen 4`` therefore do NOT resolve to their vocabulary entries and would
silently fall back to default normalisation stats. Name resolution is
consequently the pipeline's responsibility, which is what ``--marker-mapping``
is for.

Cell geometry: cells come from the CELLMEASUREMENT GeoJSON, so the mapping
from embedding row to GeoJSON feature is exact by construction. Centroids are
computed analytically from the polygons (shoelace), and per-cell footprints
are rasterised lazily inside each patch window -- never as a whole-slide
label image, which for a 30k x 30k slide would cost several GB for data the
polygons already carry.
"""

import gzip
import json
import os
import sys
import time

import numpy as np
import tifffile

# ----------------------------------------------------------------------------
# logging
# ----------------------------------------------------------------------------


def _mem_gb():
    """Current process RSS in GB, or -1.0 where unavailable (e.g. Windows)."""
    try:
        import resource

        rss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # macOS reports bytes, Linux reports KB
        if sys.platform == "darwin":
            return rss_kb / (1024**3)
        return rss_kb / (1024**2)
    except Exception:
        return -1.0


def log(msg):
    """Timestamped, RSS-annotated progress line (useful on long HPC jobs)."""
    mem = _mem_gb()
    ts = time.strftime("%H:%M:%S")
    if mem >= 0:
        print(f"[{ts}] [RSS {mem:.1f} GB] {msg}", flush=True)
    else:
        print(f"[{ts}] {msg}", flush=True)


# ----------------------------------------------------------------------------
# image scaling
# ----------------------------------------------------------------------------


def scaling_factor(dtype):
    """Intensity divisor for an image dtype.

    Matches KRONOS2's bundled ``sp_image._scaling_factor`` exactly, so patches
    reach ``model.preprocess`` on the domain it expects.
    """
    dtype = np.dtype(dtype)
    if dtype == np.uint8:
        return 255.0
    if np.issubdtype(dtype, np.unsignedinteger):
        return 65535.0
    if np.issubdtype(dtype, np.floating):
        return 400.0
    raise ValueError(
        f"no scaling factor defined for image dtype {dtype!r}; "
        "pass --max-value to override"
    )


# ----------------------------------------------------------------------------
# channel names + marker mapping
# ----------------------------------------------------------------------------


def channel_layout(tif):
    """How a TIFF stores its channels: ``("stacked" | "split", n_channels)``.

    Two layouts occur in practice. The usual one keeps every channel in a single
    ``CYX`` series. The other gives each channel its own 2D series, which some
    OME writers produce -- notably the inForm unmixing export behind the COMET
    ``_fixed`` images. A reader that only consults ``series[0]`` sees the second
    as a one-channel image, and since the patch buffer is sized from the *marker*
    list, numpy then broadcasts that single channel's pixels under every marker's
    name: an embedding of one channel wearing eight labels, with no error. Both
    the name list and the pixel reader derive the channel count from here so they
    cannot disagree about what a channel is.
    """
    series = list(getattr(tif, "series", None) or [])
    if (
        len(series) > 1
        and all(len(s.shape) == 2 for s in series)
        and len({tuple(s.shape) for s in series}) == 1
    ):
        return "split", len(series)
    if series:
        shape = tuple(series[0].shape)
        return "stacked", (shape[0] if len(shape) >= 3 else 1)
    return "stacked", (len(tif.pages) if len(tif.pages) > 1 else 1)


def get_channel_names(tiff_path):
    """Channel names from OME-XML, then ImageJ labels, then numbered fallback.

    Reads metadata only -- no pixel decode. The result always holds exactly one
    entry per image channel: a channel the metadata leaves unnamed gets a
    synthesised name rather than being dropped. Dropping shortened the list while
    leaving the image unchanged, sliding every later name onto the wrong channel
    index -- and those indices are what ``ChannelReader`` reads pixels by.
    """
    with tifffile.TiffFile(tiff_path) as tif:
        _, n_channels = channel_layout(tif)

        names = []
        if tif.is_ome:
            try:
                from xml.etree import ElementTree as ET

                root = ET.fromstring(tif.ome_metadata)
                ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
                channels = root.findall(".//ome:Channel", ns) or root.findall(
                    ".//Channel"
                )
                names = [c.get("Name") or c.get("ID") for c in channels]
            except Exception:
                names = []

        if not names and getattr(tif, "imagej_metadata", None):
            if "Labels" in tif.imagej_metadata:
                names = list(tif.imagej_metadata["Labels"])

    names = [name or f"Channel_{i}" for i, name in enumerate(names[:n_channels])]
    names += [f"Channel_{i}" for i in range(len(names), n_channels)]
    return names


def load_marker_mapping(arg):
    """Parse ``--marker-mapping``: a JSON file path or an inline JSON string.

    Keys beginning with ``_`` are treated as comments and dropped, so the
    shipped example file can document itself.
    """
    if not arg:
        return {}
    if os.path.isfile(arg):
        with open(arg) as fh:
            raw = json.load(fh)
    else:
        try:
            raw = json.loads(arg)
        except json.JSONDecodeError:
            log(f"WARNING: --marker-mapping is neither a file nor valid JSON: {arg}")
            return {}
    return {k: v for k, v in raw.items() if not k.startswith("_")}


def apply_marker_mapping(names, mapping):
    """Rename channels for the model only; the stored image is untouched.

    Returns ``(resolved_names, applied)`` where ``applied`` records the
    substitutions actually made, for the marker report.
    """
    resolved, applied = [], []
    for name in names:
        if name in mapping:
            resolved.append(mapping[name])
            applied.append((name, mapping[name]))
        else:
            resolved.append(name)
    return resolved, applied


def select_channels(channel_names, exclude=(), protect=None):
    """Channel indices to show the model, skipping ``exclude``.

    Returns ``(keep_indices, kept_names, excluded_names)``. Matching is exact and
    case sensitive against the image's own channel names, *before* any
    ``--marker-mapping`` rename -- you cannot map a channel you intend to skip,
    and keeping both keyed on the same strings means the channel list is read
    once.

    Scoped to this step alone. A skipped channel is simply never read, so it
    contributes nothing to the embedding; it stays in the image, stays available
    to every other process, and keeps whatever intensity measurements
    CELLMEASUREMENT already wrote for it. Nothing is deleted from anything.

    This is for channels that are not markers: autofluorescence, blanks, empty
    cycles. Such a channel has no vocabulary entry and no mapping can give it
    one, and ``--allow-novel-defaults`` is the wrong instrument because it still
    *feeds* the channel to the model, normalised with default statistics.

    Raises ValueError rather than skipping anything silently: an unmatched name
    means a typo produced a full-panel embedding that looks perfectly fine.
    """
    wanted, seen = [], set()
    for name in exclude or ():
        name = name.strip()
        if name and name not in seen:
            seen.add(name)
            wanted.append(name)

    missing = [name for name in wanted if name not in channel_names]
    if missing:
        raise ValueError(
            f"cannot exclude {missing}: no such channel. This image has: "
            f"{list(channel_names)} (matching is case sensitive, and uses the "
            "image's own channel names, not mapped ones)"
        )

    if protect and protect in seen:
        raise ValueError(
            f"cannot exclude {protect!r}: it is the nuclear marker, which "
            "KRONOS2 takes as preferred_dapi. That is the model's only alias "
            "mechanism and the DAPI statistics anchor the normalisation of "
            "every other channel"
        )

    keep_indices, kept_names, excluded_names = [], [], []
    for index, name in enumerate(channel_names):
        if name in seen:
            excluded_names.append(name)
        else:
            keep_indices.append(index)
            kept_names.append(name)

    if not keep_indices:
        raise ValueError(
            f"excluding {len(excluded_names)} channel(s) would leave nothing to "
            "embed"
        )

    return keep_indices, kept_names, excluded_names


# ----------------------------------------------------------------------------
# cell geometry from GeoJSON
# ----------------------------------------------------------------------------


def _open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def _flatten_ring(coords):
    """Descend nested coordinate arrays until a list of [x, y] pairs remains."""
    while (
        len(coords) > 0
        and isinstance(coords[0], list)
        and len(coords[0]) > 0
        and isinstance(coords[0][0], list)
    ):
        coords = coords[0]
    return coords


def _polygon_centroid(ring):
    """Area-weighted polygon centroid (shoelace), with a degenerate fallback.

    Analytic, so no raster is needed to locate a cell.
    """
    n = len(ring) - 1  # closed ring: last point repeats the first
    area = cx = cy = 0.0
    for i in range(n):
        x0, y0 = ring[i][0], ring[i][1]
        x1, y1 = ring[i + 1][0], ring[i + 1][1]
        cross = x0 * y1 - x1 * y0
        area += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    area *= 0.5
    if abs(area) > 1e-10:
        return cx / (6.0 * area), cy / (6.0 * area)
    xs = [p[0] for p in ring]
    ys = [p[1] for p in ring]
    return sum(xs) / len(xs), sum(ys) / len(ys)


def read_cells(geojson_path):
    """Load cell features from a CELLMEASUREMENT GeoJSON.

    Returns ``(geojson, features, centroids, rings)`` where ``centroids`` are
    ``(x, y)`` in level-0 pixels and ``rings`` are the exterior rings used for
    per-patch footprint rasterisation. Row order is GeoJSON feature order, so
    embedding row i always corresponds to ``features[i]``.
    """
    with _open_maybe_gzip(geojson_path) as fh:
        geojson = json.load(fh)

    features, centroids, rings = [], [], []
    for feat in geojson.get("features", []):
        if feat.get("properties", {}).get("objectType") != "cell":
            continue
        ring = _flatten_ring(feat["geometry"]["coordinates"])
        if len(ring) < 4:  # not a closed polygon
            continue
        features.append(feat)
        centroids.append(_polygon_centroid(ring))
        rings.append(ring)

    return geojson, features, centroids, rings


# ----------------------------------------------------------------------------
# lazy windowed image access
# ----------------------------------------------------------------------------


class ChannelReader:
    """Windowed reader over a multi-channel TIFF.

    Prefers a lazy zarr view (``tifffile.aszarr``) so memory scales with the
    patch batch rather than the slide: a 40-marker 30k x 30k uint16 image is
    ~72 GB if materialised, which does not fit the process_multi allocation.
    Falls back to an eager read for layouts zarr cannot open, and says which
    path it took.

    Handles both layouts ``channel_layout`` distinguishes. In the split layout
    each channel is its own 2D series, so there is one zarr view per channel
    instead of one indexed by a channel axis; ``read_window`` returns the same
    ``(len(channels), h, w)`` block either way and callers need not care.

    ``channels`` is a list of channel *indices* -- a subset is normal (see
    ``select_channels``) and only those channels are ever read.
    """

    def __init__(self, tiff_path, channels):
        self.path = str(tiff_path)
        self.channels = list(channels)
        self._store = None
        self._stores = []
        self._zarr = None
        self._zarrs = None
        self._eager = None

        with tifffile.TiffFile(self.path) as tif:
            self.layout, self.n_channels = channel_layout(tif)
            series = tif.series[0]
            self.shape = tuple(series.shape)
            self.dtype = np.dtype(series.dtype)
            self.height, self.width = self.shape[-2], self.shape[-1]

        try:
            import zarr

            if self.layout == "split":
                # One store per channel, each owning its own file handle exactly
                # as the stacked path's imread() does -- borrowing them from a
                # TiffFile that closes here would leave their lifetime unclear.
                # Only the kept channels are opened.
                self._zarrs = []
                for index in self.channels:
                    store = tifffile.imread(self.path, aszarr=True, series=index)
                    self._stores.append(store)
                    self._zarrs.append(self._level0(zarr.open(store, mode="r")))
            else:
                self._store = tifffile.imread(self.path, aszarr=True)
                self._zarr = self._level0(zarr.open(self._store, mode="r"))
            self.mode = "lazy"
        except Exception as exc:  # noqa: BLE001 - fall back, never fail here
            log(f"  lazy zarr view unavailable ({exc}); falling back to full read")
            self.close()
            self._zarr, self._zarrs = None, None
            if self.layout == "split":
                # imread() would hand back series[0] alone, so name the series.
                with tifffile.TiffFile(self.path) as tif:
                    self._eager = np.stack(
                        [tif.series[index].asarray() for index in self.channels]
                    )
            else:
                arr = tifffile.imread(self.path)
                if arr.ndim == 2:
                    arr = arr[np.newaxis]
                arr = arr.reshape(-1, arr.shape[-2], arr.shape[-1])
                self._eager = arr[self.channels]
            self.mode = "eager"

    @staticmethod
    def _level0(view):
        """A pyramidal OME-TIFF opens as a group; level 0 is full resolution."""
        if hasattr(view, "keys") and not hasattr(view, "shape"):
            return view[sorted(view.keys())[0]]
        return view

    def read_window(self, y0, y1, x0, x1):
        """``(len(channels), y1-y0, x1-x0)`` window, clipped to the image.

        The channel count is checked rather than trusted: the patch buffer is
        sized from the marker list, so a window with too few channels would be
        broadcast into it by numpy and every marker would silently receive the
        same pixels. Failing here costs one shape comparison per patch.
        """
        window = self._read_window(y0, y1, x0, x1)
        if window.shape[0] != len(self.channels):
            raise ValueError(
                f"read {window.shape[0]} channel(s) for {len(self.channels)} "
                f"requested from a {self.layout}-layout image; refusing to "
                "broadcast one channel across several markers"
            )
        return window

    def _read_window(self, y0, y1, x0, x1):
        if self._eager is not None:
            return self._eager[:, y0:y1, x0:x1]
        if self._zarrs is not None:
            return np.stack(
                [np.asarray(arr[y0:y1, x0:x1]) for arr in self._zarrs], axis=0
            )
        arr = self._zarr
        if arr.ndim == 2:
            return np.asarray(arr[y0:y1, x0:x1])[np.newaxis]
        return np.stack(
            [np.asarray(arr[c, y0:y1, x0:x1]) for c in self.channels], axis=0
        )

    def close(self):
        for store in [self._store, *self._stores]:
            if store is None:
                continue
            try:
                store.close()
            except Exception:
                pass
        self._store, self._stores = None, []


# ----------------------------------------------------------------------------
# patch extraction
# ----------------------------------------------------------------------------


class CellPatchDataset:
    """Cell-centred patches, extracted lazily and optionally cell-isolated.

    ``isolate=True`` zeroes every pixel not inside the target cell's polygon,
    so the embedding describes that cell rather than its neighbourhood. At a
    64 px patch and ~0.28 um/px that box spans ~18 um and typically contains
    several neighbouring cells, so the distinction is material -- and it is
    what makes the output a *cell* embedding.

    Deliberately not a ``torch.utils.data.Dataset`` subclass so this module
    imports without torch; a DataLoader only needs ``__len__``/``__getitem__``.
    """

    def __init__(
        self,
        reader,
        centroids,
        rings,
        patch_size,
        divisor,
        isolate=True,
    ):
        self.reader = reader
        self.centroids = centroids
        self.rings = rings
        self.patch_size = int(patch_size)
        self.divisor = float(divisor)
        self.isolate = bool(isolate)
        self.half = self.patch_size // 2

    def __len__(self):
        return len(self.centroids)

    def _footprint(self, index, x0, y0):
        """Binary mask of cell ``index`` within its own patch window.

        Uses Pillow's C polygon filler rather than rasterio/GDAL. Only one small
        polygon is filled per patch, so the geospatial stack bought nothing --
        and pulling it in caused a libtiff/libjpeg symbol clash inside the conda
        container. Pillow is the same approach CORAL takes, for the same reason.
        """
        from PIL import Image, ImageDraw

        size = self.patch_size
        # Shift into patch-local coordinates so only this window is filled.
        local = [(float(px) - x0, float(py) - y0) for px, py in self.rings[index]]
        try:
            img = Image.new("L", (size, size), 0)
            ImageDraw.Draw(img).polygon(local, fill=1, outline=1)
            return np.asarray(img, dtype=bool)
        except Exception:
            # A malformed polygon must not drop the cell; keep the raw box.
            return np.ones((size, size), dtype=bool)

    def __getitem__(self, index):
        size = self.patch_size
        cx, cy = self.centroids[index]
        x0 = int(round(cx)) - self.half
        y0 = int(round(cy)) - self.half
        x1, y1 = x0 + size, y0 + size

        patch = np.zeros((len(self.reader.channels), size, size), dtype=self.reader.dtype)

        # Clip to the image; border cells are kept and zero-padded rather than
        # dropped, so every cell receives an embedding.
        sx0, sy0 = max(0, x0), max(0, y0)
        sx1, sy1 = min(self.reader.width, x1), min(self.reader.height, y1)
        if sx1 > sx0 and sy1 > sy0:
            win = self.reader.read_window(sy0, sy1, sx0, sx1)
            patch[:, sy0 - y0 : sy0 - y0 + (sy1 - sy0), sx0 - x0 : sx0 - x0 + (sx1 - sx0)] = win

        if self.isolate:
            patch = patch * self._footprint(index, x0, y0)[np.newaxis]

        return (patch.astype(np.float32) / self.divisor), index


def numpy_collate(items):
    """Stack ``(patch, index)`` pairs, keeping patches as numpy for the model."""
    patches = np.stack([p for p, _ in items])
    idx = np.asarray([i for _, i in items])
    return patches, idx


# ----------------------------------------------------------------------------
# output
# ----------------------------------------------------------------------------


def write_geojson_with_embeddings(
    geojson, features, embeddings, output_path, prefix="kronos_emb_"
):
    """Attach per-cell embeddings to their GeoJSON features and write.

    Row i maps to ``features[i]`` by construction (both derive from GeoJSON
    feature order), so no spatial matching is involved. Any embedding keys from
    a previous run are cleared first, making re-runs idempotent instead of
    accumulating stale columns.
    """
    dim = embeddings.shape[1]
    columns = [f"{prefix}{i}" for i in range(dim)]

    cleared = 0
    for feat in features:
        measurements = feat.get("properties", {}).get("measurements", {})
        stale = [k for k in measurements if k.startswith(prefix)]
        for k in stale:
            del measurements[k]
        if stale:
            cleared += 1
    if cleared:
        log(f"  cleared previous embeddings from {cleared} cells")

    for row, feat in enumerate(features):
        measurements = feat["properties"].setdefault("measurements", {})
        values = embeddings[row]
        for col, value in zip(columns, values):
            measurements[col] = float(value)

    with _open_maybe_gzip(output_path, "wt") as fh:
        json.dump(geojson, fh, separators=(",", ":"))

    log(f"  wrote {len(features)} cells x {dim}-d embeddings to {output_path}")
    return len(features)


def write_marker_report(
    path, tiff, channel_names, resolved, applied, extra_lines=(), keep_indices=None
):
    """Human-readable record of how each channel was presented to the model.

    ``channel_names`` is the image's full channel list, so the header count stays
    honest about the file; ``resolved`` covers only the channels actually
    embedded, and ``keep_indices`` says which ones those were. ``None`` means
    everything was kept, which reproduces the report as it was before exclusion
    existed.

    Withheld channels appear twice on purpose: in place, so "what happened to
    channel 8" needs no counting, and as a trailing block that stays greppable on
    a forty-marker panel. Both say "not shown to KRONOS2" rather than "removed",
    because that is the truth -- the channel is still in the image.
    """
    if keep_indices is None:
        keep_indices = list(range(len(channel_names)))
    by_index = dict(zip(keep_indices, resolved))
    excluded = [name for i, name in enumerate(channel_names) if i not in by_index]

    with open(path, "w") as fh:
        fh.write("KRONOS2 Marker Report\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Image: {tiff}\n")
        if excluded:
            fh.write(
                f"Channels: {len(channel_names)} "
                f"({len(by_index)} embedded, {len(excluded)} not shown)\n\n"
            )
        else:
            fh.write(f"Channels: {len(channel_names)}\n\n")
        fh.write("Channel -> name given to KRONOS2\n")
        for index, raw in enumerate(channel_names):
            if index not in by_index:
                fh.write(f"  {raw} -> (not shown to KRONOS2)\n")
                continue
            res = by_index[index]
            note = "  (mapped)" if raw != res else ""
            fh.write(f"  {raw} -> {res}{note}\n")
        if applied:
            fh.write("\nMappings applied via --marker-mapping:\n")
            for src, dst in applied:
                fh.write(f"  {src} -> {dst}\n")
        if excluded:
            fh.write(
                "\nWithheld from KRONOS2 via --exclude-marker "
                "(still present in the image):\n"
            )
            for name in excluded:
                fh.write(f"  {name}\n")
        for line in extra_lines:
            fh.write(f"\n{line}\n")
