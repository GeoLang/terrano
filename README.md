# Terrano

[![CI](https://github.com/GeoLang/terrano/actions/workflows/ci.yml/badge.svg)](https://github.com/GeoLang/terrano/actions)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)

Raster algebra and terrain analysis engine for the GeoLang GIS stack.

[Documentation](https://geolang.github.io/terrano/) · [GitHub](https://github.com/GeoLang/terrano)

## Features

- **Terrain analysis**: hillshade (0 to 255), slope in degrees and aspect in degrees clockwise from north, all from Horn's 3x3 gradient
- **Contours**: lines at a given interval and base, segments joined into connected lines
- **Viewshed**: line-of-sight visibility from an observer cell, ray cast to every cell inside a radius
- **Watershed delineation**: every cell of a D8 flow direction raster labelled by the pit, flat or edge cell it drains to. Pour points are found by tracing, not supplied
- **Flow direction**: D8 flow direction from a DEM
- **Flow accumulation**: count of upstream cells draining into each cell
- **Stream ordering**: Strahler order for cells above a flow accumulation threshold
- **Sink filling**: `fill_sinks` raises depressions so every cell drains to an edge
- **Map algebra**: unary (add, multiply, sqrt, abs, log) and binary (add, subtract, multiply, divide, min, max) operations
- **Reclassification**: value ranges mapped to class values
- **Polygonize**: connected runs of equal cells traced as polygon rings with holes, for a classified raster
- **Rasterize**: polygons burnt onto a grid by cell centre, holes cut out, the inverse of polygonize
- **Focal statistics**: moving-window min/max/mean/sum/std/median/majority/range over a square or circular neighbourhood
- **Zonal statistics**: per-zone summary of one raster grouped by the labels of another
- **GeoTIFF I/O**: read and write GeoTIFF rasters with origin, pixel size and EPSG code
- **Multi-band rasters**: `BandedRaster` holds RGB/RGBA or any band set on one grid, written and read as a multi-band GeoTIFF in any `SampleFormat`
- **Cloud Optimized GeoTIFF (COG)**: tiled writing with overview pyramids (raw or deflate), validated in CI against GDAL's `validate_cloud_optimized_geotiff.py` (full check, overviews asserted), and windowed reads through the `RangeRead` trait (`CogReader` fetches only the tiles a window touches, so an implementation backed by HTTP `Range` requests reads a remote file). Writes any `SampleFormat` (u8, i8, u16, i16, u32, i32, f32, f64) with the geo tags, GDAL_NODATA, and the IFD and tile ordering the COG spec calls for. Reads real-world single-band COGs: deflate, horizontal and floating-point predictors, integer and float sample types, GDAL nodata mapped to NaN. Multi-band COGs are pixel-interleaved through `write_cog_bands` and `CogReader::read_window_bands`. The writer runs in the browser too, via `writeCog` and `writeCogBands` in terrano-wasm
- **GRIB2 and NetCDF**: `grib::scan_grib` lists the messages in a GRIB2 file and `grib::decode_grib_message` decodes simple-packed ones (other packings return an error). `netcdf::read_netcdf_metadata` and `netcdf::read_netcdf_variable` read NetCDF classic and 64-bit offset files, not NetCDF-4
- **EO time-series**: `RasterStack` for multi-temporal analysis: composites (mean/median/min/max/standard deviation), linear trend fitting, change detection, anomaly z-scores, phenology metrics, normalized difference indices (NDVI, NDWI, etc.)
- **Browser build**: `terrano-wasm` is a wasm-bindgen surface over terrano-core taking flat f64 buffers: `hillshade`, `slope`, `aspect`, `fillSinks`, `reclassify`, `applyUnary`, `applyBinary`, `normalizedDifference`, `contours`, `polygonize`, `rasterize`, `focalStats`, `zonalStats`, `writeCog` and `writeCogBands`. terrano-core depends only on thiserror and flate2, so it builds for `wasm32-unknown-unknown` with no C toolchain

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

// Cells visible from 2 elevation units above cell (100, 100), out to 5000 map units
let seen = viewshed(&dem, 100, 100, 2.0, 5000.0);

// Write output
let mut out = std::fs::File::create("slopes.tif").unwrap();
write_geotiff(&slopes, &meta, &mut out).unwrap();
```

## COG sample formats

`CogParams::format` picks the sample type, `SampleFormat::F64` by default. An
8-bit image written as `U8` is an eighth the size of the same image as `F64`.
In terrano-wasm, `writeCog` takes the format as a string, `"u8"` through `"f64"`.

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

On integer formats a block average can round onto the
nodata value, which turns a cell of real data into an absent one in that
overview level. Pick a nodata at the edge of the range, not in the middle of
the data.

## CLI

The `terrano` binary runs on a DEM it generates itself and reads or writes no
raster files. Use `terrano-core` for real work. `stats` prints elevation, slope
and aspect at the centre cell, `hillshade` prints the hillshade value there.

```sh
terrano stats --width 10 --height 10
terrano hillshade --azimuth 315 --altitude 45
```

## License

AGPL-3.0-or-later, see [LICENSE](LICENSE).

Copyright (C) 2026 Grok Image Compression Inc.
