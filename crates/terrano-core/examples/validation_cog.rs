//! Writes the COGs that CI hands to GDAL's `validate_cloud_optimized_geotiff.py`.
//!
//! Both images are larger than 512 pixels on a side, so the validator fails
//! them unless the overviews are really present.

use std::{env, fs, path::PathBuf};
use terrano_core::{BandedRaster, CogParams, Raster, SampleFormat, write_cog, write_cog_bands};

fn gradient(width: usize, height: usize, scale: f64) -> Raster {
    let data = (0..width * height)
        .map(|i| (i % 251) as f64 * scale)
        .collect();
    Raster::from_vec(width, height, data, 1.0, -9999.0).unwrap()
}

fn main() {
    let dir = PathBuf::from(
        env::args()
            .nth(1)
            .expect("usage: validation_cog <output dir>"),
    );
    fs::create_dir_all(&dir).unwrap();

    let single = CogParams {
        overview_levels: 3,
        deflate: true,
        nodata: Some(-9999.0),
        format: SampleFormat::F32,
        ..CogParams::default()
    };
    let mut file = fs::File::create(dir.join("single.tif")).unwrap();
    write_cog(&gradient(1030, 700, 1.0), &single, &mut file).unwrap();

    let rgb = CogParams {
        overview_levels: 2,
        nodata: None,
        format: SampleFormat::U8,
        ..CogParams::default()
    };
    let bands = BandedRaster::new(vec![
        gradient(800, 600, 1.0),
        gradient(800, 600, 0.5),
        gradient(800, 600, 0.25),
    ])
    .unwrap();
    let mut file = fs::File::create(dir.join("rgb.tif")).unwrap();
    write_cog_bands(&bands, &rgb, &mut file).unwrap();
}
