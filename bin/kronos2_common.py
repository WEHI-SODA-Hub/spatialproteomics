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


def get_channel_names(tiff_path):
    """Channel names from OME-XML, then ImageJ labels, then numbered fallback.

    Reads metadata only -- no pixel decode.
    """
    with tifffile.TiffFile(tiff_path) as tif:
        if tif.series:
            shape = tif.series[0].shape
            n_channels = shape[0] if len(shape) >= 3 else 1
        else:
            n_channels = len(tif.pages) if len(tif.pages) > 1 else 1

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

        if not names:
            names = [f"Channel_{i}" for i in range(n_channels)]

    return [n for n in names if n is not None]


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
    """

    def __init__(self, tiff_path, channels):
        self.path = str(tiff_path)
        self.channels = list(channels)
        self._store = None
        self._zarr = None
        self._eager = None

        with tifffile.TiffFile(self.path) as tif:
            series = tif.series[0]
            self.shape = tuple(series.shape)
            self.dtype = np.dtype(series.dtype)

        if len(self.shape) == 2:
            self.n_channels, self.height, self.width = 1, self.shape[0], self.shape[1]
        else:
            self.n_channels, self.height, self.width = (
                self.shape[0],
                self.shape[-2],
                self.shape[-1],
            )

        try:
            import zarr

            self._store = tifffile.imread(self.path, aszarr=True)
            self._zarr = zarr.open(self._store, mode="r")
            # A pyramidal OME-TIFF opens as a group; level 0 is full resolution.
            if hasattr(self._zarr, "keys") and not hasattr(self._zarr, "shape"):
                self._zarr = self._zarr[sorted(self._zarr.keys())[0]]
            self.mode = "lazy"
        except Exception as exc:  # noqa: BLE001 - fall back, never fail here
            log(f"  lazy zarr view unavailable ({exc}); falling back to full read")
            arr = tifffile.imread(self.path)
            if arr.ndim == 2:
                arr = arr[np.newaxis]
            arr = arr.reshape(-1, arr.shape[-2], arr.shape[-1])
            self._eager = arr[self.channels]
            self.mode = "eager"

    def read_window(self, y0, y1, x0, x1):
        """``(len(channels), y1-y0, x1-x0)`` window, clipped to the image."""
        if self._eager is not None:
            return self._eager[:, y0:y1, x0:x1]
        arr = self._zarr
        if arr.ndim == 2:
            return np.asarray(arr[y0:y1, x0:x1])[np.newaxis]
        return np.stack(
            [np.asarray(arr[c, y0:y1, x0:x1]) for c in self.channels], axis=0
        )

    def close(self):
        if self._store is not None:
            try:
                self._store.close()
            except Exception:
                pass


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
        self._rasterize = None
        if self.isolate:
            from rasterio import features as _feat

            self._rasterize = _feat.rasterize

    def __len__(self):
        return len(self.centroids)

    def _footprint(self, index, x0, y0):
        """Binary mask of cell ``index`` within its own patch window."""
        from shapely.geometry import Polygon

        size = self.patch_size
        ring = self.rings[index]
        # Shift into patch-local coordinates so only this window is rasterised.
        local = [(px - x0, py - y0) for px, py in ring]
        try:
            poly = Polygon(local)
            if not poly.is_valid:
                poly = poly.buffer(0)
            if poly.is_empty:
                return np.ones((size, size), dtype=bool)
            return self._rasterize(
                [(poly, 1)],
                out_shape=(size, size),
                fill=0,
                dtype=np.uint8,
                all_touched=False,
            ).astype(bool)
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


def write_marker_report(path, tiff, channel_names, resolved, applied, extra_lines=()):
    """Human-readable record of how each channel was presented to the model."""
    with open(path, "w") as fh:
        fh.write("KRONOS2 Marker Report\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Image: {tiff}\n")
        fh.write(f"Channels: {len(channel_names)}\n\n")
        fh.write("Channel -> name given to KRONOS2\n")
        for raw, res in zip(channel_names, resolved):
            note = "  (mapped)" if raw != res else ""
            fh.write(f"  {raw} -> {res}{note}\n")
        if applied:
            fh.write("\nMappings applied via --marker-mapping:\n")
            for src, dst in applied:
                fh.write(f"  {src} -> {dst}\n")
        for line in extra_lines:
            fh.write(f"\n{line}\n")
