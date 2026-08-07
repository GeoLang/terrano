use crate::Raster;

/// Compute hillshade from a DEM raster.
///
/// # Parameters
/// - `dem`: Input elevation raster
/// - `azimuth`: Sun azimuth in degrees (0=N, 90=E, 180=S, 270=W)
/// - `altitude`: Sun altitude angle in degrees above horizon
pub fn hillshade(dem: &Raster, azimuth: f64, altitude: f64) -> Raster {
    let az = azimuth.to_radians();
    let alt = altitude.to_radians();
    let mut result = Raster::new(dem.width(), dem.height(), dem.cell_size, dem.nodata);

    for row in 1..dem.height() - 1 {
        for col in 1..dem.width() - 1 {
            let (dzdx, dzdy) = gradient(dem, row, col);
            if dzdx.is_nan() || dzdy.is_nan() {
                continue;
            }

            let slope_rad = (dzdx * dzdx + dzdy * dzdy).sqrt().atan();
            let aspect_rad = dzdy.atan2(-dzdx);

            let hs = 255.0
                * (alt.cos() * slope_rad.cos()
                    + alt.sin() * slope_rad.sin() * (az - aspect_rad).cos());
            result.set(row, col, hs.clamp(0.0, 255.0));
        }
    }
    result
}

/// Compute slope in degrees from a DEM raster.
pub fn slope(dem: &Raster) -> Raster {
    let mut result = Raster::new(dem.width(), dem.height(), dem.cell_size, dem.nodata);

    for row in 1..dem.height() - 1 {
        for col in 1..dem.width() - 1 {
            let (dzdx, dzdy) = gradient(dem, row, col);
            if dzdx.is_nan() || dzdy.is_nan() {
                continue;
            }
            let s = (dzdx * dzdx + dzdy * dzdy).sqrt().atan().to_degrees();
            result.set(row, col, s);
        }
    }
    result
}

/// Compute aspect in degrees from a DEM raster (0=N, clockwise).
pub fn aspect(dem: &Raster) -> Raster {
    let mut result = Raster::new(dem.width(), dem.height(), dem.cell_size, dem.nodata);

    for row in 1..dem.height() - 1 {
        for col in 1..dem.width() - 1 {
            let (dzdx, dzdy) = gradient(dem, row, col);
            if dzdx.is_nan() || dzdy.is_nan() {
                continue;
            }
            if dzdx == 0.0 && dzdy == 0.0 {
                // aspect is undefined on flat cells, leave nodata
                continue;
            }
            let mut a = dzdy.atan2(-dzdx).to_degrees();
            if a < 0.0 {
                a += 360.0;
            }
            result.set(row, col, a);
        }
    }
    result
}

/// Compute gradient (dz/dx, dz/dy) using Horn's method (3x3 neighbourhood).
fn gradient(dem: &Raster, row: usize, col: usize) -> (f64, f64) {
    let get = |r: usize, c: usize| -> f64 { dem.get(r, c).unwrap_or(dem.nodata) };

    let a = get(row - 1, col - 1);
    let b = get(row - 1, col);
    let c = get(row - 1, col + 1);
    let d = get(row, col - 1);
    let f = get(row, col + 1);
    let g = get(row + 1, col - 1);
    let h = get(row + 1, col);
    let i = get(row + 1, col + 1);

    // Check for nodata
    for v in [a, b, c, d, f, g, h, i] {
        if dem.is_nodata(v) {
            return (f64::NAN, f64::NAN);
        }
    }

    let cell = dem.cell_size;
    let dzdx = ((c + 2.0 * f + i) - (a + 2.0 * d + g)) / (8.0 * cell);
    let dzdy = ((g + 2.0 * h + i) - (a + 2.0 * b + c)) / (8.0 * cell);
    (dzdx, dzdy)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_dem() -> Raster {
        // 5x5 flat DEM at elevation 100
        Raster::from_vec(5, 5, vec![100.0; 25], 10.0, -9999.0).unwrap()
    }

    fn sloped_dem() -> Raster {
        // 5x5 DEM with linear slope in x-direction
        let mut data = vec![0.0; 25];
        for row in 0..5 {
            for col in 0..5 {
                data[row * 5 + col] = col as f64 * 10.0;
            }
        }
        Raster::from_vec(5, 5, data, 1.0, -9999.0).unwrap()
    }

    #[test]
    fn test_slope_flat() {
        let result = slope(&flat_dem());
        // Interior cells should have zero slope
        assert!((result.get(2, 2).unwrap()).abs() < 1e-10);
    }

    #[test]
    fn test_slope_linear() {
        let result = slope(&sloped_dem());
        let s = result.get(2, 2).unwrap();
        // dz/dx = 10/1 = 10 per cell, slope = atan(10) ≈ 84.3°
        assert!((s - 84.289).abs() < 0.1);
    }

    #[test]
    fn test_hillshade_flat() {
        let result = hillshade(&flat_dem(), 315.0, 45.0);
        // Flat surface at 45° sun altitude → hillshade ≈ 255 * cos(0) * ... ≈ 180
        let hs = result.get(2, 2).unwrap();
        assert!(hs > 170.0 && hs < 190.0);
    }
}
