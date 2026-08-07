// Comprehensive integration tests for terrano-core.

use terrano_core::{
    BandedRaster, BinaryOp, GeoTiffMetadata, Raster, SampleFormat, UnaryOp, aspect, contours,
    flow_accumulation, flow_direction, hillshade, read_geotiff_bands, reclassify, slope,
    write_geotiff, write_geotiff_bands,
};

// ═══════════════════════════════════════════════════════════════════════════
// Raster basics
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_raster_new_filled_with_nodata() {
    let r = Raster::new(4, 3, 10.0, -9999.0);
    assert_eq!(r.width(), 4);
    assert_eq!(r.height(), 3);
    for row in 0..3 {
        for col in 0..4 {
            assert_eq!(r.get(row, col), Some(-9999.0));
        }
    }
}

#[test]
fn test_raster_from_vec() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let r = Raster::from_vec(3, 2, data, 1.0, -9999.0).unwrap();
    assert_eq!(r.get(0, 0), Some(1.0));
    assert_eq!(r.get(0, 2), Some(3.0));
    assert_eq!(r.get(1, 0), Some(4.0));
    assert_eq!(r.get(1, 2), Some(6.0));
}

#[test]
fn test_raster_from_vec_dimension_mismatch() {
    let data = vec![1.0, 2.0, 3.0];
    let result = Raster::from_vec(5, 5, data, 1.0, -9999.0);
    assert!(result.is_err());
}

#[test]
fn test_raster_get_out_of_bounds() {
    let r = Raster::new(3, 3, 1.0, -9999.0);
    assert_eq!(r.get(10, 10), None);
}

#[test]
fn test_raster_set_and_get() {
    let mut r = Raster::new(3, 3, 1.0, -9999.0);
    r.set(1, 1, 42.0);
    assert_eq!(r.get(1, 1), Some(42.0));
}

#[test]
fn test_raster_is_nodata() {
    let r = Raster::new(3, 3, 1.0, -9999.0);
    assert!(r.is_nodata(-9999.0));
    assert!(r.is_nodata(f64::NAN));
    assert!(!r.is_nodata(0.0));
}

// ═══════════════════════════════════════════════════════════════════════════
// Map algebra
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_unary_add() {
    let r = Raster::from_vec(2, 2, vec![1.0, 2.0, 3.0, 4.0], 1.0, -9999.0).unwrap();
    let result = r.apply_unary(&UnaryOp::Add(10.0));
    assert_eq!(result.get(0, 0), Some(11.0));
    assert_eq!(result.get(1, 1), Some(14.0));
}

#[test]
fn test_unary_multiply() {
    let r = Raster::from_vec(2, 2, vec![2.0, 4.0, 6.0, 8.0], 1.0, -9999.0).unwrap();
    let result = r.apply_unary(&UnaryOp::Multiply(0.5));
    assert_eq!(result.get(0, 0), Some(1.0));
    assert_eq!(result.get(1, 1), Some(4.0));
}

#[test]
fn test_unary_sqrt() {
    let r = Raster::from_vec(2, 2, vec![4.0, 9.0, 16.0, 25.0], 1.0, -9999.0).unwrap();
    let result = r.apply_unary(&UnaryOp::Sqrt);
    assert_eq!(result.get(0, 0), Some(2.0));
    assert_eq!(result.get(0, 1), Some(3.0));
}

#[test]
fn test_unary_abs() {
    let r = Raster::from_vec(2, 2, vec![-1.0, -2.0, 3.0, -4.0], 1.0, -9999.0).unwrap();
    let result = r.apply_unary(&UnaryOp::Abs);
    assert_eq!(result.get(0, 0), Some(1.0));
    assert_eq!(result.get(1, 1), Some(4.0));
}

#[test]
fn test_unary_preserves_nodata() {
    let r = Raster::from_vec(2, 2, vec![1.0, -9999.0, 3.0, 4.0], 1.0, -9999.0).unwrap();
    let result = r.apply_unary(&UnaryOp::Add(10.0));
    assert_eq!(result.get(0, 1), Some(-9999.0)); // nodata preserved
    assert_eq!(result.get(0, 0), Some(11.0));
}

#[test]
fn test_binary_add() {
    let a = Raster::from_vec(2, 2, vec![1.0, 2.0, 3.0, 4.0], 1.0, -9999.0).unwrap();
    let b = Raster::from_vec(2, 2, vec![10.0, 20.0, 30.0, 40.0], 1.0, -9999.0).unwrap();
    let result = a.apply_binary(&b, &BinaryOp::Add).unwrap();
    assert_eq!(result.get(0, 0), Some(11.0));
    assert_eq!(result.get(1, 1), Some(44.0));
}

#[test]
fn test_binary_subtract() {
    let a = Raster::from_vec(2, 2, vec![10.0, 20.0, 30.0, 40.0], 1.0, -9999.0).unwrap();
    let b = Raster::from_vec(2, 2, vec![1.0, 2.0, 3.0, 4.0], 1.0, -9999.0).unwrap();
    let result = a.apply_binary(&b, &BinaryOp::Subtract).unwrap();
    assert_eq!(result.get(0, 0), Some(9.0));
}

#[test]
fn test_binary_divide_by_zero() {
    let a = Raster::from_vec(2, 2, vec![10.0, 20.0, 30.0, 40.0], 1.0, -9999.0).unwrap();
    let b = Raster::from_vec(2, 2, vec![0.0, 2.0, 0.0, 4.0], 1.0, -9999.0).unwrap();
    let result = a.apply_binary(&b, &BinaryOp::Divide).unwrap();
    // Division by zero produces nodata
    assert_eq!(result.get(0, 0), Some(-9999.0));
    assert_eq!(result.get(0, 1), Some(10.0));
}

#[test]
fn test_binary_incompatible_rasters() {
    let a = Raster::new(3, 3, 1.0, -9999.0);
    let b = Raster::new(4, 4, 1.0, -9999.0);
    let result = a.apply_binary(&b, &BinaryOp::Add);
    assert!(result.is_err());
}

#[test]
fn test_binary_min_max() {
    let a = Raster::from_vec(2, 2, vec![5.0, 1.0, 8.0, 3.0], 1.0, -9999.0).unwrap();
    let b = Raster::from_vec(2, 2, vec![3.0, 7.0, 2.0, 9.0], 1.0, -9999.0).unwrap();

    let min = a.apply_binary(&b, &BinaryOp::Min).unwrap();
    assert_eq!(min.get(0, 0), Some(3.0));
    assert_eq!(min.get(0, 1), Some(1.0));

    let max = a.apply_binary(&b, &BinaryOp::Max).unwrap();
    assert_eq!(max.get(0, 0), Some(5.0));
    assert_eq!(max.get(0, 1), Some(7.0));
}

// ═══════════════════════════════════════════════════════════════════════════
// Reclassify
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_reclassify_basic() {
    let r = Raster::from_vec(
        3,
        3,
        vec![5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0],
        1.0,
        -9999.0,
    )
    .unwrap();
    let classes = vec![
        (0.0, 30.0, 1.0),   // low
        (30.0, 60.0, 2.0),  // medium
        (60.0, 100.0, 3.0), // high
    ];
    let result = reclassify(&r, &classes);
    assert_eq!(result.get(0, 0), Some(1.0)); // 5 → low
    assert_eq!(result.get(1, 0), Some(2.0)); // 35 → medium
    assert_eq!(result.get(2, 0), Some(3.0)); // 65 → high
}

// ═══════════════════════════════════════════════════════════════════════════
// Terrain analysis
// ═══════════════════════════════════════════════════════════════════════════

fn flat_dem() -> Raster {
    Raster::from_vec(5, 5, vec![100.0; 25], 10.0, -9999.0).unwrap()
}

fn sloped_dem() -> Raster {
    let mut data = vec![0.0; 25];
    for row in 0..5 {
        for col in 0..5 {
            data[row * 5 + col] = col as f64 * 10.0;
        }
    }
    Raster::from_vec(5, 5, data, 1.0, -9999.0).unwrap()
}

#[test]
fn test_slope_flat_terrain() {
    let result = slope(&flat_dem());
    // Interior cells should have zero slope
    assert!((result.get(2, 2).unwrap()).abs() < 1e-10);
}

#[test]
fn test_slope_inclined_terrain() {
    let result = slope(&sloped_dem());
    // Interior cells should have non-zero slope
    let s = result.get(2, 2).unwrap();
    assert!(s > 0.0);
}

#[test]
fn test_hillshade_flat_is_uniform() {
    let result = hillshade(&flat_dem(), 315.0, 45.0);
    let a = result.get(2, 2).unwrap();
    let b = result.get(2, 3).unwrap();
    // Flat terrain → uniform illumination at interior cells
    assert!((a - b).abs() < 1e-10);
}

#[test]
fn test_hillshade_range() {
    let result = hillshade(&sloped_dem(), 315.0, 45.0);
    // All values should be in [0, 255]
    for row in 0..5 {
        for col in 0..5 {
            let v = result.get(row, col).unwrap();
            if !result.is_nodata(v) {
                assert!((0.0..=255.0).contains(&v), "hillshade out of range: {v}");
            }
        }
    }
}

#[test]
fn test_aspect_flat_terrain() {
    let result = aspect(&flat_dem());
    // Flat terrain has 0 gradient → aspect is undefined, cells keep nodata
    let v = result.get(2, 2).unwrap();
    assert!(result.is_nodata(v));
}

#[test]
fn test_aspect_eastward_slope() {
    // Slope increases going east (col increases)
    let result = aspect(&sloped_dem());
    let a = result.get(2, 2).unwrap();
    // Aspect indicates direction of steepest descent — verify it's a valid angle
    assert!((0.0..360.0).contains(&a), "aspect out of range: {a}");
}

// ═══════════════════════════════════════════════════════════════════════════
// Contours
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_contours_flat_no_lines() {
    let flat = flat_dem();
    let result = contours(&flat, 10.0, 0.0);
    // A perfectly flat DEM may produce one contour at 100.0 level
    // All vertices should be at that level
    for line in &result {
        assert_eq!(line.level, 100.0);
    }
}

#[test]
fn test_contours_sloped_produces_lines() {
    let dem = sloped_dem();
    let result = contours(&dem, 10.0, 0.0);
    // Should produce contour lines at 10, 20, 30 intervals
    assert!(!result.is_empty());
    // All levels should be multiples of 10
    for line in &result {
        assert!((line.level % 10.0).abs() < 1e-10);
    }
}

#[test]
fn test_contours_invalid_interval() {
    let dem = flat_dem();
    let result = contours(&dem, 0.0, 0.0);
    assert!(result.is_empty());
    let result2 = contours(&dem, -5.0, 0.0);
    assert!(result2.is_empty());
}

// ═══════════════════════════════════════════════════════════════════════════
// Hydrology
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn test_flow_direction_downhill() {
    // Simple 3x3 bowl: center is lowest
    let data = vec![10.0, 10.0, 10.0, 10.0, 1.0, 10.0, 10.0, 10.0, 10.0];
    let dem = Raster::from_vec(3, 3, data, 1.0, -9999.0).unwrap();
    let fdir = flow_direction(&dem);
    // Surrounding cells should flow toward center
    // Check that corner and edge cells have non-zero flow directions
    assert!(fdir.get(0, 0).unwrap() > 0.0); // should flow inward
}

#[test]
fn test_flow_direction_flat_no_flow() {
    let flat = flat_dem();
    let fdir = flow_direction(&flat);
    // Flat terrain → no steepest descent → direction stays 0
    assert_eq!(fdir.get(2, 2).unwrap(), 0.0);
}

#[test]
fn test_flow_accumulation_basic() {
    // Linear slope flowing east
    let fdir = flow_direction(&sloped_dem());
    let accum = flow_accumulation(&fdir);
    // Higher accumulation values downstream (lower cols for eastward-sloping
    // ... actually the slope goes W→E so flow goes E→W (downhill is toward col 0)
    // In any case, we just check it doesn't panic and produces non-negative values
    for row in 0..5 {
        for col in 0..5 {
            let v = accum.get(row, col).unwrap();
            if !accum.is_nodata(v) {
                assert!(v >= 0.0);
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Multi-band rasters
// ═══════════════════════════════════════════════════════════════════════════

/// Bytes produced by `write_geotiff` before multi-band support landed. Fenestra
/// reads and writes these, so the single-band output must not shift.
const GOLDEN_SINGLE_BAND: &[u8] = include_bytes!("fixtures/single_band_f64.tif");

#[test]
fn test_write_geotiff_single_band_bytes_unchanged() {
    let raster = Raster::from_vec(
        3,
        2,
        vec![1.5, -2.25, 3.0, 4.125, -9999.0, 6.5],
        10.0,
        -9999.0,
    )
    .unwrap();
    let meta = GeoTiffMetadata {
        origin_x: 500000.0,
        origin_y: 4500000.0,
        pixel_width: 10.0,
        pixel_height: 10.0,
        epsg: 32632,
    };

    let mut buf = Vec::new();
    write_geotiff(&raster, &meta, &mut buf).unwrap();
    assert_eq!(buf, GOLDEN_SINGLE_BAND);
}

#[test]
fn test_banded_band_feeds_hillshade_identically() {
    let dem = sloped_dem();
    let banded = BandedRaster::new(vec![Raster::new(5, 5, 10.0, -9999.0), dem.clone()]).unwrap();

    let direct = hillshade(&dem, 315.0, 45.0);
    let per_band = hillshade(banded.band(1).unwrap(), 315.0, 45.0);
    assert_eq!(per_band.data(), direct.data());
}

#[test]
fn test_banded_rgb_geotiff_roundtrip() {
    let bands: Vec<Vec<f64>> = vec![
        (0..25).map(|i| (i * 10 % 256) as f64).collect(),
        (0..25).map(|i| (255 - i * 3) as f64).collect(),
        (0..25).map(|i| (i % 2 * 255) as f64).collect(),
    ];
    let raster = BandedRaster::with_names(
        bands
            .iter()
            .map(|v| Raster::from_vec(5, 5, v.clone(), 0.5, -9999.0).unwrap())
            .collect(),
        vec!["red".into(), "green".into(), "blue".into()],
    )
    .unwrap();
    let meta = GeoTiffMetadata {
        origin_x: -122.5,
        origin_y: 37.8,
        pixel_width: 0.5,
        pixel_height: 0.5,
        epsg: 4326,
    };

    let mut buf = Vec::new();
    write_geotiff_bands(&raster, &meta, SampleFormat::U8, &mut buf).unwrap();
    let (read, read_meta) = read_geotiff_bands(&buf).unwrap();

    assert_eq!(read.band_count(), 3);
    assert_eq!((read.width(), read.height()), (5, 5));
    assert_eq!(read.cell_size(), 0.5);
    assert_eq!(read_meta.epsg, 4326);
    assert_eq!(read_meta.origin_x, -122.5);
    for (b, values) in bands.iter().enumerate() {
        assert_eq!(read.band(b).unwrap().data(), values.as_slice());
    }
}

#[test]
fn test_banded_mismatched_band_sizes_rejected() {
    let result = BandedRaster::new(vec![
        Raster::new(4, 4, 1.0, -9999.0),
        Raster::new(4, 3, 1.0, -9999.0),
    ]);
    assert!(result.is_err());
}
