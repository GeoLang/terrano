# Terrano

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

## Usage

```rust
use terrano_core::{
    Raster, slope, hillshade, aspect, contours, flow_direction,
    flow_accumulation, watershed, read_geotiff, write_geotiff,
};

// Read a DEM
let dem = read_geotiff("elevation.tif").unwrap();

// Terrain derivatives
let slopes = slope(&dem);
let hs = hillshade(&dem, 315.0, 45.0);
let asp = aspect(&dem);

// Contour lines at 10m intervals
let lines = contours(&dem, 10.0);

// Hydrology
let flow_dir = flow_direction(&dem);
let accumulation = flow_accumulation(&flow_dir);
let basins = watershed(&flow_dir, &pour_points);

// Write output
write_geotiff("slopes.tif", &slopes, &metadata).unwrap();
```

## CLI

```sh
terrano hillshade --input dem.tif --azimuth 315 --altitude 45 --output hillshade.tif
terrano slope --input dem.tif --output slope.tif
terrano contour --input dem.tif --interval 10 --output contours.geojson
terrano flow --input dem.tif --output flow_acc.tif
terrano watershed --input dem.tif --pour-points points.geojson --output basins.tif
```

## License

AGPL-3.0-or-later
