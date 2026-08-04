# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

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

- `write_cog` produced files whose tile offsets pointed into IFD bytes (the per-IFD size bookkeeping undercounted the geo arrays), and `pack_geokeys` silently dropped the EPSG code. COGs now carry a valid GeoKeyDirectory and correct tile offsets, verified by read-back tests.

## [0.1.0] - 2026-05-30

### Added

- Initial release.
