# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

## [0.2.0] - 2026-09-02

### Added

- `viewshed` (2026-08-24): line-of-sight visibility from an observer cell. A ray
  runs from the eye to every cell inside the radius, sampling the terrain
  bilinearly once per row or column it crosses, and the cell is hidden when a
  sample stands at a steeper elevation angle than the target. Cells are
  `VIEWSHED_VISIBLE` or `VIEWSHED_HIDDEN`, and nodata past the radius, on cells
  the DEM has no elevation for, and everywhere when the observer cell has none.
  A sample exactly on the line of sight grazes it, so a flat surface stays
  visible from an eye at ground level, and a nodata sample blocks nothing.

- COG validation in CI (2026-08-14): a `cog-validate` job writes a single-band
  f32 and a multi-band u8 COG through the `validation_cog` example, runs GDAL's
  `validate_cloud_optimized_geotiff.py --full-check` on both, and asserts the
  overviews via gdalinfo, since the validator only warns when a tiled file has
  none. The README's validation claim was one manual run until now.

- `writeCogBands` in `terrano-wasm` (2026-08-13): multi-band COG encoding from
  the browser, over one flat f64 buffer holding the bands end to end plus a
  band count, so band `b` starts at `b * width * height`. It fronts
  terrano-core's `write_cog_bands`, writing pixel-interleaved tiles with
  SamplesPerPixel set and overviews block-averaged per band. `writeCog` is now
  this function at one band, so both take the same georeferencing, sample
  format and nodata.

- COG sample formats (2026-08-13): `CogParams::format` writes u8, i8, u16,
  i16, u32, i32, f32 or f64 tiles, matching what the reader already decoded.
  It defaults to f64, so existing callers are unchanged, and an 8-bit image
  written as u8 is an eighth the size it was. `SampleFormat` grew from two
  variants to those eight and now also drives `write_geotiff_bands` and
  `read_geotiff_bands`, which read and write every one of them. `writeCog` in
  `terrano-wasm` takes the format by name, e.g. `"u16"`.

  Samples round to the nearest whole number and clamp to the format's range
  rather than wrapping, and overviews are averaged in f64 and rounded only
  when a tile is encoded. Nodata on an integer format is an ordinary sample
  value set aside to mean absent: it has to be whole and inside the range,
  and NaN samples are written as it. With no nodata declared, an integer file
  rejects NaN samples rather than storing a silent zero. An f32 nodata has to
  survive the narrowing, so a value like 0.1 is rejected. Output for every
  format passes `validate_cloud_optimized_geotiff.py --full-check=yes`, and
  gdalinfo reports the declared type, the nodata value and all three overview
  levels.

- `writeCog` in `terrano-wasm` (2026-08-13): COG encoding from the browser,
  taking a flat f64 buffer plus georeferencing and returning the file bytes.
  terrano-core builds for wasm32-unknown-unknown with no C dependencies,
  flate2 resolving to pure-Rust miniz_oxide, so deflate works there too.

### Fixed

- `write_geotiff` (2026-08-13): built and discarded its output buffer twice
  before the pass that counted, so every call did the work three times.

- `focal_stats`, `zonal_stats` and `rasterize` (2026-08-04): the neighbourhood
  and grouped summaries, plus the polygon burn that feeds them.
  `focal_stats(raster, radius, shape, stat)` reports min/max/mean/sum/std/
  median/majority/range over a square or circular window, clipping at the
  raster edge and skipping nodata neighbours, leaving nodata cells alone
  rather than growing the data area. `zonal_stats(values, zones)` returns one
  `ZoneStats` row per distinct zone label, counting a cell only where both
  grids carry data. `rasterize(polygons, ...)` burns `RegionPolygon`s onto a
  grid by cell centre with holes cut out, the inverse of `polygonize`, which
  is what turns a set of boundaries into a zone raster. Note the y axis
  differs: `polygonize` counts rows downward, `rasterize` reads north-up.
  All three are in `terrano-wasm`; `focalStats` takes no cell size because
  its window is measured in cells, and `zonalStats` returns flat rows of
  [zone, count, min, max, mean, sum, std, median].

- `polygonize` (2026-08-04): raster to vector, the classification counterpart
  to `contours`. Connected runs of exactly equal cells (4-connected, nodata
  bounding them) become `RegionPolygon`s whose rings follow cell corners,
  exterior first then its holes, with collinear runs collapsed so a rectangle
  is five vertices rather than one per cell edge. A region that pinches at a
  corner comes back as one polygon per lobe rather than a self-crossing ring.
  Exposed in `terrano-wasm` flat-encoded as [value, ring_count, (vertex_count,
  x, y, ...) per ring].

- `terrano-wasm` (2026-08-04): wasm-bindgen surface over terrano-core for the
  browser, free functions over flat f64 buffers (hillshade, slope, aspect,
  fill_sinks, reclassify, unary/binary algebra, normalized difference,
  contours flat-encoded as [level, count, x, y, ...]). Built with wasm-pack
  --target web, vendored into viewtopia (src/raster/wasm). `Raster::into_data`
  added so results move out without a copy.

- `RangeRead::read_ranges` reads several byte ranges in one call, results in request order, with a sequential default so existing transports keep working. `CogReader::read_window_bands` now collects every tile a window touches into one `read_ranges` call, so a multiplexing transport fetches them concurrently.
- Multi-band COGs: `write_cog_bands` writes a `BandedRaster` pixel-interleaved (SamplesPerPixel per band, PlanarConfiguration chunky, overviews block-averaged per band), `CogReader::read_window_bands` reads a window back as one `Raster` per band. `CogReader::open` accepts any sample count as long as every band shares one bit depth and sample format. `write_cog` and `read_window` are unchanged, single-band output is byte for byte the same, and `read_window` on a multi-band file errors and names `read_window_bands`.
- `CogReader` decodes real-world single-band COGs: deflate tiles, horizontal and floating-point predictors, uint/int/float sample types widened to f64, GDAL nodata mapped to NaN. `write_cog` optionally compresses tiles with deflate (`CogParams.deflate`).
- `CogReader`: windowed COG reading over a `RangeRead` byte-range seam with overview selection (`select_level`), fetching only the tiles a window touches. `RangeRead` ships for `std::fs::File` and `&[u8]`, remote streaming plugs in via HTTP `Range` requests.

### Fixed

- `write_cog` output is now a real cloud optimized GeoTIFF (2026-08-13), not
  just a tiled one. Overview IFDs carry NewSubfileType, without which GDAL
  read the pyramid as four unrelated pages and reported no overviews at all,
  and tile data is now laid out smallest overview first with full resolution
  last, the order the COG spec requires so a zoomed-out reader can stop early.
  Files pass `validate_cloud_optimized_geotiff.py --full-check` against GDAL
  3.11, single and multi-band, raw and deflate.
- `write_cog` writes GDAL_NODATA, set through the new `CogParams.nodata`, and
  substitutes it for NaN samples so the value the file declares absent is the
  one actually stored. Defaults to `nan`, which is what the writer already
  padded partial edge tiles with. The reader has always mapped the declared
  nodata back to NaN.
- `write_cog` GeoKeyDirectory gained GTModelTypeGeoKey and GTRasterTypeGeoKey
  alongside the CRS key, and multi-band output declares ExtraSamples, which
  silences a libtiff warning on every read of a multi-band terrano COG.
- `write_cog` and `write_cog_bands` no longer require `Seek`, so a browser can
  write straight into a `Vec<u8>`. Oversized images now error instead of
  silently truncating their 32-bit offsets.
- `aspect` (2026-08-06) now delivers its documented compass convention:
  degrees clockwise from north, where it previously returned the raw
  counterclockwise-from-east atan2 angle. Flat cells (zero gradient) keep
  nodata instead of reading as south-facing, the same treatment the
  nodata-gradient ring already got.
- `write_cog` produced files whose tile offsets pointed into IFD bytes (the per-IFD size bookkeeping undercounted the geo arrays), and `pack_geokeys` silently dropped the EPSG code. COGs now carry a valid GeoKeyDirectory and correct tile offsets, verified by read-back tests.

## [0.1.0] - 2026-05-30

### Added

- Initial release.
