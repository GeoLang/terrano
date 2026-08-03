# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

- `CogReader` decodes real-world single-band COGs: deflate tiles, horizontal and floating-point predictors, uint/int/float sample types widened to f64, GDAL nodata mapped to NaN. `write_cog` optionally compresses tiles with deflate (`CogParams.deflate`).
- `CogReader`: windowed COG reading over a `RangeRead` byte-range seam with overview selection (`select_level`), fetching only the tiles a window touches. `RangeRead` ships for `std::fs::File` and `&[u8]`, remote streaming plugs in via HTTP `Range` requests.

### Fixed

- `write_cog` produced files whose tile offsets pointed into IFD bytes (the per-IFD size bookkeeping undercounted the geo arrays), and `pack_geokeys` silently dropped the EPSG code. COGs now carry a valid GeoKeyDirectory and correct tile offsets, verified by read-back tests.

## [0.1.0] - 2026-05-30

### Added

- Initial release.
