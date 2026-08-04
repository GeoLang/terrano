# Terrano

[![CI](https://github.com/GeoLang/terrano/actions/workflows/ci.yml/badge.svg)](https://github.com/GeoLang/terrano/actions)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)

Raster algebra and terrain analysis engine for the GeoLang GIS stack.

[Documentation](https://geolang.github.io/terrano/) · [GitHub](https://github.com/GeoLang/terrano)

## Features

- **Terrain analysis** — Hillshade, slope (degrees, Horn's method), aspect (0–360°)
- **Contour generation** — Extract contour lines at configurable intervals with segment connectivity
- **Watershed delineation** — D8 watershed boundaries from pour points
- **Flow direction** — D8 single-direction flow routing from DEM
- **Flow accumulation** — Upstream area/cell count per pixel
- **Stream ordering** — Strahler stream order from flow accumulation
- **Sink filling** — Remove depressions for hydrologically-correct DEMs
- **Map algebra** — Unary (add, multiply, sqrt, abs, log) and binary (add, subtract, multiply, divide, min, max) operations
- **Reclassification** — Value-range-based class assignment
- **GeoTIFF I/O** — Read and write GeoTIFF rasters with CRS metadata
- **Multi-band rasters** — `BandedRaster` holds RGB/RGBA or any band set on one grid, written and read as a multi-band GeoTIFF with 8-bit or 64-bit float samples
- **Cloud Optimized GeoTIFF (COG)** — tiled writing with overview pyramids (raw or deflate), and windowed reads over a byte-range seam (`CogReader` fetches only the tiles a window touches, wire it to `Range` requests for remote streaming). Reads real-world single-band COGs: deflate, horizontal and floating-point predictors, integer and float sample types, GDAL nodata mapped to NaN. Multi-band COGs are pixel-interleaved through `write_cog_bands` and `CogReader::read_window_bands`
- **GRIB2 and NetCDF** — message scanning and variable reads for gridded weather and climate data
- **EO time-series** — `RasterStack` for multi-temporal analysis: composites (mean/median/max), linear trend fitting, change detection, anomaly z-scores, phenology metrics, normalized difference indices (NDVI, NDWI, etc.)

## Usage

```rust
use terrano_core::{
    Raster, slope, hillshade, aspect, contours, flow_direction,
    flow_accumulation, watershed, read_geotiff, write_geotiff,
};

// Read a DEM from bytes
let bytes = std::fs::read("elevation.tif").unwrap();
let (dem, meta) = read_geotiff(&bytes).unwrap();

// Terrain derivatives
let slopes = slope(&dem);
let hs = hillshade(&dem, 315.0, 45.0);
let asp = aspect(&dem);

// Contour lines every 10m, starting at 0
let lines = contours(&dem, 10.0, 0.0);

// Hydrology
let flow_dir = flow_direction(&dem);
let accumulation = flow_accumulation(&flow_dir);
let basins = watershed(&flow_dir);

// Write output
let mut out = std::fs::File::create("slopes.tif").unwrap();
write_geotiff(&slopes, &meta, &mut out).unwrap();
```

## CLI

The CLI is a demo harness over a synthetic DEM, it does not read or write raster files yet.
Use `terrano-core` directly for real work.

```sh
terrano stats --width 10 --height 10
terrano hillshade --azimuth 315 --altitude 45
```

## License

AGPL-3.0-or-later, see [LICENSE](LICENSE).

Copyright (C) 2026 Grok Image Compression Inc.
