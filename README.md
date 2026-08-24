# Terrano

[![CI](https://github.com/GeoLang/terrano/actions/workflows/ci.yml/badge.svg)](https://github.com/GeoLang/terrano/actions)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)

Raster algebra and terrain analysis engine for the GeoLang GIS stack.

[Documentation](https://geolang.github.io/terrano/) · [GitHub](https://github.com/GeoLang/terrano)

## Features

- **Terrain analysis** — Hillshade, slope (degrees, Horn's method), aspect (0–360°)
- **Contour generation** — Extract contour lines at configurable intervals with segment connectivity
- **Viewshed** — line-of-sight visibility from an observer cell, ray cast to every cell inside a radius
- **Watershed delineation** — D8 watershed boundaries from pour points
- **Flow direction** — D8 single-direction flow routing from DEM
- **Flow accumulation** — Upstream area/cell count per pixel
- **Stream ordering** — Strahler stream order from flow accumulation
- **Sink filling** — Remove depressions for hydrologically-correct DEMs
- **Map algebra** — Unary (add, multiply, sqrt, abs, log) and binary (add, subtract, multiply, divide, min, max) operations
- **Reclassification** — Value-range-based class assignment
- **Polygonize** — Connected runs of equal cells traced as polygon rings with holes, for a classified raster
- **Rasterize** — Polygons burnt onto a grid by cell centre, holes cut out, the inverse of polygonize
- **Focal statistics** — Moving-window min/max/mean/sum/std/median/majority/range over a square or circular neighbourhood
- **Zonal statistics** — Per-zone summary of one raster grouped by the labels of another
- **GeoTIFF I/O** — Read and write GeoTIFF rasters with CRS metadata
- **Multi-band rasters** — `BandedRaster` holds RGB/RGBA or any band set on one grid, written and read as a multi-band GeoTIFF in any `SampleFormat`
- **Cloud Optimized GeoTIFF (COG)** — tiled writing with overview pyramids (raw or deflate), validated in CI against GDAL's `validate_cloud_optimized_geotiff.py` (full check, overviews asserted), and windowed reads over a byte-range seam (`CogReader` fetches only the tiles a window touches, wire it to `Range` requests for remote streaming). Writes any `SampleFormat` (u8, i8, u16, i16, u32, i32, f32, f64) with the geo tags, GDAL_NODATA, and the IFD and tile ordering the COG spec calls for. Reads real-world single-band COGs: deflate, horizontal and floating-point predictors, integer and float sample types, GDAL nodata mapped to NaN. Multi-band COGs are pixel-interleaved through `write_cog_bands` and `CogReader::read_window_bands`. The writer runs in the browser too, via `writeCog` and `writeCogBands` in terrano-wasm
- **GRIB2 and NetCDF** — message scanning and variable reads for gridded weather and climate data
- **EO time-series** — `RasterStack` for multi-temporal analysis: composites (mean/median/max), linear trend fitting, change detection, anomaly z-scores, phenology metrics, normalized difference indices (NDVI, NDWI, etc.)

## Usage

```rust
use terrano_core::{
    Raster, slope, hillshade, aspect, contours, polygonize, reclassify, flow_direction,
    flow_accumulation, watershed, viewshed, read_geotiff, write_geotiff,
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

// Classified raster to polygons
let regions = polygonize(&reclassify(&dem, &[(0.0, 500.0, 1.0), (500.0, 2000.0, 2.0)]));

// Hydrology
let flow_dir = flow_direction(&dem);
let accumulation = flow_accumulation(&flow_dir);
let basins = watershed(&flow_dir);

// What an observer 2 m above cell (100, 100) can see within 5 km
let seen = viewshed(&dem, 100, 100, 2.0, 5000.0);

// Write output
let mut out = std::fs::File::create("slopes.tif").unwrap();
write_geotiff(&slopes, &meta, &mut out).unwrap();
```

## COG sample formats

`CogParams::format` picks the sample type, `SampleFormat::F64` by default. An
8-bit image written as `U8` is an eighth the size of the same image as `F64`,
which is what the browser path through `writeCog` cares about.

```rust
use terrano_core::{CogParams, SampleFormat, write_cog};

let params = CogParams {
    format: SampleFormat::U8,
    nodata: Some(255.0),
    ..CogParams::default()
};
```

Terrano holds every raster as f64 in memory, so writing a narrower format
converts on the way out. Values round to the nearest whole number and clamp to
the format's range, so an out-of-range sample is pinned to the nearest end
rather than wrapping. Overviews are averaged in f64 and only rounded when a
tile is encoded, which is why a u32 pyramid does not overflow.

Nodata means different things per format:

- On `F32` and `F64`, `nodata` is substituted for NaN samples and declared in
  GDAL_NODATA. `None` writes no tag and leaves NaN in the file. An `F32` nodata
  has to survive the narrowing, so `0.1` is rejected and `-9999.0` is fine.
- On the integer formats there is no NaN, so nodata is an ordinary sample value
  set aside to mean absent. It has to be whole and inside the format's range,
  and every NaN in the source is written as it. `None` declares no absent value
  at all: a NaN sample is then an error rather than a silent zero, and only the
  padding past the image edge is filled with zero.

One caveat inherent to integer rasters: a block average can round onto the
nodata value, which turns a cell of real data into an absent one in that
overview level. Pick a nodata at the edge of the range, not in the middle of
the data.

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
