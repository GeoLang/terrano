# Terrano

Raster algebra and terrain analysis engine for the TileTopia-HQ GIS stack.

## Features

- **Raster type** — 2D grid with nodata handling and cell size
- **Map algebra** — unary (add, multiply, sqrt, abs, log) and binary (add, subtract, multiply, divide, min, max) cell-by-cell operations
- **Reclassification** — value-range-based reclassification
- **Terrain analysis** — hillshade, slope, aspect (Horn's method, 3×3 kernel)

## Usage

```rust
use terrano_core::{Raster, slope, hillshade};

let dem = Raster::from_vec(5, 5, elevation_data, 10.0, -9999.0).unwrap();
let slopes = slope(&dem);
let hs = hillshade(&dem, 315.0, 45.0);
```

## CLI

```sh
terrano stats --width 20 --height 20
terrano hillshade --azimuth 315 --altitude 45
```

## License

AGPL-3.0-or-later
