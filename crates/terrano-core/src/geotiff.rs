use crate::{Error, Raster};
use std::io::{self, Write};

/// GeoTIFF metadata for georeferencing.
#[derive(Debug, Clone)]
pub struct GeoTiffMetadata {
    /// X coordinate of the top-left corner
    pub origin_x: f64,
    /// Y coordinate of the top-left corner
    pub origin_y: f64,
    /// Pixel width in map units
    pub pixel_width: f64,
    /// Pixel height in map units (positive value; image grows downward)
    pub pixel_height: f64,
    /// EPSG code for the CRS (e.g., 4326 for WGS84, 32632 for UTM zone 32N)
    pub epsg: u16,
}

/// Write a raster as a minimal GeoTIFF file.
///
/// Produces a valid TIFF with:
/// - 64-bit float samples
/// - GeoKey directory (ModelTiepointTag, ModelPixelScaleTag, GeoKeyDirectoryTag)
/// - Single strip
///
/// # Arguments
/// * `raster` — the raster data to write
/// * `meta` — georeferencing metadata
/// * `writer` — output destination
pub fn write_geotiff<W: Write>(
    raster: &Raster,
    meta: &GeoTiffMetadata,
    writer: &mut W,
) -> Result<(), Error> {
    let width = raster.width() as u32;
    let height = raster.height() as u32;
    let sample_size: u32 = 8; // f64
    let strip_byte_count = width * height * sample_size;

    // We'll write a little-endian TIFF
    let mut buf: Vec<u8> = Vec::new();

    // TIFF Header (8 bytes)
    buf.write_all(b"II")?; // Little-endian
    write_u16(&mut buf, 42)?; // Magic
    write_u32(&mut buf, 8)?; // Offset to first IFD (immediately after header)

    // Wait - IFD comes after header. Let's plan the layout:
    // [0..8] = header
    // [8..] = IFD
    // After IFD: tag data (tiepoint, pixel scale, geokeys), then strip data

    // GeoTIFF-specific tag data we'll need to write after the IFD:
    let tiepoint: [f64; 6] = [0.0, 0.0, 0.0, meta.origin_x, meta.origin_y, 0.0];
    let pixel_scale: [f64; 3] = [meta.pixel_width, meta.pixel_height, 0.0];

    // GeoKeyDirectoryTag: version, revision, minor, numkeys, then key entries
    // Key 1024 = GTModelTypeGeoKey: 1=Projected, 2=Geographic
    // Key 2048 = GeographicTypeGeoKey (for geographic CRS)
    // Key 3072 = ProjectedCSTypeGeoKey (for projected CRS)
    let model_type: u16 = if meta.epsg < 5000 { 2 } else { 1 };
    let geo_keys: Vec<u16> = if model_type == 2 {
        vec![
            1, 1, 0, 2, // version 1.1.0, 2 keys
            1024, 0, 1, model_type, // GTModelTypeGeoKey = Geographic
            2048, 0, 1, meta.epsg, // GeographicTypeGeoKey
        ]
    } else {
        vec![
            1, 1, 0, 2, // version 1.1.0, 2 keys
            1024, 0, 1, model_type, // GTModelTypeGeoKey = Projected
            3072, 0, 1, meta.epsg, // ProjectedCSTypeGeoKey
        ]
    };

    // Number of IFD entries
    let num_tags: u16 = 11;
    let ifd_start: u32 = 8;
    let ifd_size: u32 = 2 + (num_tags as u32) * 12 + 4; // count + entries + next_ifd
    let data_start: u32 = ifd_start + ifd_size;

    // Layout of extra data after IFD:
    let tiepoint_offset = data_start;
    let tiepoint_size = 48u32; // 6 * 8
    let pixel_scale_offset = tiepoint_offset + tiepoint_size;
    let pixel_scale_size = 24u32; // 3 * 8
    let geo_keys_offset = pixel_scale_offset + pixel_scale_size;
    let geo_keys_size = (geo_keys.len() as u32) * 2;
    let strip_offset = geo_keys_offset + geo_keys_size;

    // Re-build buffer from scratch
    buf.clear();

    // TIFF header
    buf.write_all(b"II")?;
    write_u16(&mut buf, 42)?;
    write_u32(&mut buf, ifd_start)?;

    // IFD entry count
    write_u16(&mut buf, num_tags)?;

    // Tag entries (must be sorted by tag number)
    // 256: ImageWidth
    write_ifd_entry(&mut buf, 256, 3, 1, width)?;
    // 257: ImageLength
    write_ifd_entry(&mut buf, 257, 3, 1, height)?;
    // 258: BitsPerSample
    write_ifd_entry(&mut buf, 258, 3, 1, 64)?;
    // 259: Compression (1 = None)
    write_ifd_entry(&mut buf, 259, 3, 1, 1)?;
    // 262: PhotometricInterpretation (1 = MinIsBlack)
    write_ifd_entry(&mut buf, 262, 3, 1, 1)?;
    // 273: StripOffsets
    write_ifd_entry(&mut buf, 273, 4, 1, strip_offset)?;
    // 277: SamplesPerPixel
    write_ifd_entry(&mut buf, 277, 3, 1, 1)?;
    // 278: RowsPerStrip
    write_ifd_entry(&mut buf, 278, 3, 1, height)?;
    // 279: StripByteCounts
    write_ifd_entry(&mut buf, 279, 4, 1, strip_byte_count)?;
    // 339: SampleFormat (3 = IEEE floating point)
    write_ifd_entry(&mut buf, 339, 3, 1, 3)?;
    // 33550: ModelPixelScaleTag (DOUBLE, count=3)
    write_ifd_entry(&mut buf, 33550, 12, 3, pixel_scale_offset)?;

    // We need to also fit ModelTiepointTag (33922) and GeoKeyDirectoryTag (34735)
    // But we only have 11 tags. Let me add them.
    // Actually let me recalculate with 13 tags total.

    // Let me redo this properly with all needed tags.
    buf.clear();

    let num_tags: u16 = 13;
    let ifd_size: u32 = 2 + (num_tags as u32) * 12 + 4;
    let data_start: u32 = ifd_start + ifd_size;

    let tiepoint_offset = data_start;
    let pixel_scale_offset = tiepoint_offset + tiepoint_size;
    let geo_keys_offset = pixel_scale_offset + pixel_scale_size;
    let geo_keys_size = (geo_keys.len() as u32) * 2;
    let strip_offset = geo_keys_offset + geo_keys_size;

    // TIFF header
    buf.write_all(b"II")?;
    write_u16(&mut buf, 42)?;
    write_u32(&mut buf, ifd_start)?;

    // IFD
    write_u16(&mut buf, num_tags)?;

    // Sorted by tag number
    write_ifd_entry(&mut buf, 256, 3, 1, width)?; // ImageWidth
    write_ifd_entry(&mut buf, 257, 3, 1, height)?; // ImageLength
    write_ifd_entry(&mut buf, 258, 3, 1, 64)?; // BitsPerSample
    write_ifd_entry(&mut buf, 259, 3, 1, 1)?; // Compression=None
    write_ifd_entry(&mut buf, 262, 3, 1, 1)?; // PhotometricInterpretation=MinIsBlack
    write_ifd_entry(&mut buf, 273, 4, 1, strip_offset)?; // StripOffsets
    write_ifd_entry(&mut buf, 277, 3, 1, 1)?; // SamplesPerPixel
    write_ifd_entry(&mut buf, 278, 3, 1, height)?; // RowsPerStrip
    write_ifd_entry(&mut buf, 279, 4, 1, strip_byte_count)?; // StripByteCounts
    write_ifd_entry(&mut buf, 339, 3, 1, 3)?; // SampleFormat=Float
    write_ifd_entry(&mut buf, 33550, 12, 3, pixel_scale_offset)?; // ModelPixelScaleTag
    write_ifd_entry(&mut buf, 33922, 12, 6, tiepoint_offset)?; // ModelTiepointTag
    write_ifd_entry(&mut buf, 34735, 3, geo_keys.len() as u32, geo_keys_offset)?;
    // GeoKeyDirectoryTag

    // Next IFD offset (0 = no more IFDs)
    write_u32(&mut buf, 0)?;

    // Extra data: tiepoint
    for &val in &tiepoint {
        write_f64(&mut buf, val)?;
    }

    // pixel scale
    for &val in &pixel_scale {
        write_f64(&mut buf, val)?;
    }

    // geo keys
    for &val in &geo_keys {
        write_u16(&mut buf, val)?;
    }

    // Strip data (raster values as f64, row-major)
    for row in 0..raster.height() {
        for col in 0..raster.width() {
            let val = raster.get(row, col).unwrap_or(raster.nodata);
            write_f64(&mut buf, val)?;
        }
    }

    writer.write_all(&buf)?;
    Ok(())
}

/// Read a minimal GeoTIFF file (uncompressed, single band, f64).
///
/// This is a basic reader supporting the subset written by `write_geotiff`.
pub fn read_geotiff(data: &[u8]) -> Result<(Raster, GeoTiffMetadata), Error> {
    if data.len() < 8 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "file too small").into());
    }

    let le = data[0] == b'I' && data[1] == b'I';
    if !le {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "only little-endian TIFF supported",
        )
        .into());
    }

    let magic = read_u16_at(data, 2);
    if magic != 42 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "not a TIFF file").into());
    }

    let ifd_offset = read_u32_at(data, 4) as usize;
    if ifd_offset + 2 > data.len() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid IFD offset").into());
    }

    let num_entries = read_u16_at(data, ifd_offset) as usize;
    let mut width = 0u32;
    let mut height = 0u32;
    let mut strip_offset = 0u32;
    let mut tiepoint_offset = 0u32;
    let mut pixel_scale_offset = 0u32;
    let mut geo_keys_offset = 0u32;
    let mut geo_keys_count = 0u32;

    for i in 0..num_entries {
        let entry_off = ifd_offset + 2 + i * 12;
        let tag = read_u16_at(data, entry_off);
        let _typ = read_u16_at(data, entry_off + 2);
        let count = read_u32_at(data, entry_off + 4);
        let value = read_u32_at(data, entry_off + 8);

        match tag {
            256 => width = value,
            257 => height = value,
            273 => strip_offset = value,
            33550 => pixel_scale_offset = value,
            33922 => tiepoint_offset = value,
            34735 => {
                geo_keys_count = count;
                geo_keys_offset = value;
            }
            _ => {}
        }
    }

    if width == 0 || height == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "missing width/height").into());
    }

    // Read pixel scale
    let px_w = read_f64_at(data, pixel_scale_offset as usize);
    let px_h = read_f64_at(data, pixel_scale_offset as usize + 8);

    // Read tiepoint
    let origin_x = read_f64_at(data, tiepoint_offset as usize + 24);
    let origin_y = read_f64_at(data, tiepoint_offset as usize + 32);

    // Read EPSG from geokeys
    let mut epsg = 0u16;
    if geo_keys_count >= 8 {
        // Skip header (4 shorts), then read key entries
        let num_keys = read_u16_at(data, geo_keys_offset as usize + 6) as usize;
        for k in 0..num_keys {
            let base = geo_keys_offset as usize + 8 + k * 8;
            let key_id = read_u16_at(data, base);
            let key_val = read_u16_at(data, base + 6);
            if key_id == 2048 || key_id == 3072 {
                epsg = key_val;
            }
        }
    }

    // Read raster data
    let nodata = -9999.0;
    let mut raster = Raster::new(width as usize, height as usize, px_w, nodata);
    let base = strip_offset as usize;
    for row in 0..height as usize {
        for col in 0..width as usize {
            let off = base + (row * width as usize + col) * 8;
            let val = read_f64_at(data, off);
            raster.set(row, col, val);
        }
    }

    let meta = GeoTiffMetadata {
        origin_x,
        origin_y,
        pixel_width: px_w,
        pixel_height: px_h,
        epsg,
    };

    Ok((raster, meta))
}

// --- Helper functions ---

fn write_u16(buf: &mut Vec<u8>, val: u16) -> Result<(), Error> {
    buf.write_all(&val.to_le_bytes())?;
    Ok(())
}

fn write_u32(buf: &mut Vec<u8>, val: u32) -> Result<(), Error> {
    buf.write_all(&val.to_le_bytes())?;
    Ok(())
}

fn write_f64(buf: &mut Vec<u8>, val: f64) -> Result<(), Error> {
    buf.write_all(&val.to_le_bytes())?;
    Ok(())
}

fn write_ifd_entry(
    buf: &mut Vec<u8>,
    tag: u16,
    typ: u16,
    count: u32,
    value: u32,
) -> Result<(), Error> {
    write_u16(buf, tag)?;
    write_u16(buf, typ)?;
    write_u32(buf, count)?;
    write_u32(buf, value)?;
    Ok(())
}

fn read_u16_at(data: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes([data[offset], data[offset + 1]])
}

fn read_u32_at(data: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes([
        data[offset],
        data[offset + 1],
        data[offset + 2],
        data[offset + 3],
    ])
}

fn read_f64_at(data: &[u8], offset: usize) -> f64 {
    let bytes: [u8; 8] = data[offset..offset + 8].try_into().unwrap();
    f64::from_le_bytes(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geotiff_roundtrip() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let raster = Raster::from_vec(3, 2, data.clone(), 10.0, -9999.0).unwrap();

        let meta = GeoTiffMetadata {
            origin_x: 500000.0,
            origin_y: 4500000.0,
            pixel_width: 10.0,
            pixel_height: 10.0,
            epsg: 32632,
        };

        let mut buf = Vec::new();
        write_geotiff(&raster, &meta, &mut buf).unwrap();

        let (read_raster, read_meta) = read_geotiff(&buf).unwrap();
        assert_eq!(read_raster.width(), 3);
        assert_eq!(read_raster.height(), 2);
        assert!((read_meta.origin_x - 500000.0).abs() < 1e-6);
        assert!((read_meta.origin_y - 4500000.0).abs() < 1e-6);
        assert!((read_meta.pixel_width - 10.0).abs() < 1e-6);
        assert_eq!(read_meta.epsg, 32632);

        for row in 0..2 {
            for col in 0..3 {
                let expected = data[row * 3 + col];
                let actual = read_raster.get(row, col).unwrap();
                assert!(
                    (actual - expected).abs() < 1e-10,
                    "mismatch at ({},{}): expected {}, got {}",
                    row,
                    col,
                    expected,
                    actual
                );
            }
        }
    }

    #[test]
    fn test_geotiff_geographic_crs() {
        let raster = Raster::from_vec(2, 2, vec![0.0; 4], 1.0, -9999.0).unwrap();
        let meta = GeoTiffMetadata {
            origin_x: -180.0,
            origin_y: 90.0,
            pixel_width: 1.0,
            pixel_height: 1.0,
            epsg: 4326,
        };

        let mut buf = Vec::new();
        write_geotiff(&raster, &meta, &mut buf).unwrap();
        let (_, read_meta) = read_geotiff(&buf).unwrap();
        assert_eq!(read_meta.epsg, 4326);
    }
}
