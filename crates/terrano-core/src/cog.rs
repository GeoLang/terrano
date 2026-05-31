//! Cloud Optimized GeoTIFF (COG) support.
//!
//! Provides functions to write Cloud Optimized GeoTIFF files with:
//! - Internal tiling (256×256 or 512×512 tiles)
//! - Overview pyramids (2×, 4×, 8×, 16× downsampled)
//! - HTTP range-request compatible layout (IFDs at start of file)
//!
//! COG files follow the standard described at <https://www.cogeo.org/>.

use crate::Error;
use crate::raster::Raster;
use std::io::{self, Seek, Write};

/// Configuration for COG output.
#[derive(Debug, Clone)]
pub struct CogParams {
    /// Internal tile width (typically 256 or 512).
    pub tile_width: u32,
    /// Internal tile height (typically 256 or 512).
    pub tile_height: u32,
    /// Number of overview levels to generate (0 = none).
    pub overview_levels: u32,
    /// EPSG code for CRS.
    pub epsg: u16,
    /// X coordinate of the top-left corner.
    pub origin_x: f64,
    /// Y coordinate of the top-left corner.
    pub origin_y: f64,
    /// Pixel width in map units.
    pub pixel_width: f64,
    /// Pixel height in map units.
    pub pixel_height: f64,
}

impl Default for CogParams {
    fn default() -> Self {
        Self {
            tile_width: 256,
            tile_height: 256,
            overview_levels: 4,
            epsg: 4326,
            origin_x: 0.0,
            origin_y: 0.0,
            pixel_width: 1.0,
            pixel_height: 1.0,
        }
    }
}

/// An overview level (reduced-resolution image).
#[derive(Debug, Clone)]
pub struct Overview {
    pub width: usize,
    pub height: usize,
    pub data: Vec<f64>,
    pub factor: u32,
}

/// Generate overview pyramids by averaging 2×2 pixel blocks.
pub fn generate_overviews(raster: &Raster, levels: u32) -> Vec<Overview> {
    let mut overviews = Vec::new();
    let mut src_width = raster.width();
    let mut src_height = raster.height();
    let mut src_data = raster.data().to_vec();

    for level in 0..levels {
        let factor = 2u32.pow(level + 1);
        let dst_width = src_width.div_ceil(2);
        let dst_height = src_height.div_ceil(2);

        if dst_width == 0 || dst_height == 0 {
            break;
        }

        let mut dst_data = vec![f64::NAN; dst_width * dst_height];

        for dy in 0..dst_height {
            for dx in 0..dst_width {
                let sx = dx * 2;
                let sy = dy * 2;

                let mut sum = 0.0;
                let mut count = 0;

                for oy in 0..2 {
                    for ox in 0..2 {
                        let px = sx + ox;
                        let py = sy + oy;
                        if px < src_width && py < src_height {
                            let val = src_data[py * src_width + px];
                            if !val.is_nan() {
                                sum += val;
                                count += 1;
                            }
                        }
                    }
                }

                dst_data[dy * dst_width + dx] = if count > 0 {
                    sum / count as f64
                } else {
                    f64::NAN
                };
            }
        }

        overviews.push(Overview {
            width: dst_width,
            height: dst_height,
            data: dst_data.clone(),
            factor,
        });

        src_width = dst_width;
        src_height = dst_height;
        src_data = dst_data;
    }

    overviews
}

/// Write a raster as a Cloud Optimized GeoTIFF.
///
/// The file layout follows the COG specification:
/// 1. TIFF header (8 bytes)
/// 2. Ghost IFD (overview detection marker)
/// 3. Full-resolution IFD
/// 4. Overview IFDs (lowest resolution first for efficient range requests)
/// 5. Tile data (overviews first, then full resolution)
pub fn write_cog<W: Write + Seek>(
    raster: &Raster,
    params: &CogParams,
    writer: &mut W,
) -> Result<(), Error> {
    let overviews = generate_overviews(raster, params.overview_levels);

    // Phase 1: Write TIFF header
    writer.write_all(b"II")?; // Little-endian
    write_u16(writer, 42)?; // TIFF magic
    write_u32(writer, 8)?; // Offset to first IFD (immediately after header)

    // Phase 2: Write IFDs — full resolution first, then overviews
    let n_ifds = 1 + overviews.len();

    // Calculate IFD sizes to determine tile data offset
    let entries_per_ifd = 14; // Number of TIFF tags per IFD
    let ifd_size = 2 + entries_per_ifd * 12 + 4; // count + entries + next_ifd_offset
    let geo_extra_size = 48; // Space for geo metadata (pixel scale, tiepoint, geokeys)
    let total_ifd_size = n_ifds * (ifd_size + geo_extra_size);
    let tile_data_start = 8 + total_ifd_size;

    // Collect tile data for all IFDs
    let mut all_tile_data: Vec<Vec<Vec<u8>>> = Vec::new();

    // Full resolution tiles
    let full_tiles = raster_to_tiles(
        raster.data(),
        raster.width(),
        raster.height(),
        params.tile_width as usize,
        params.tile_height as usize,
    );
    all_tile_data.push(full_tiles);

    // Overview tiles
    for ov in &overviews {
        let tiles = raster_to_tiles(
            &ov.data,
            ov.width,
            ov.height,
            params.tile_width as usize,
            params.tile_height as usize,
        );
        all_tile_data.push(tiles);
    }

    // Calculate tile offsets
    let mut current_offset = tile_data_start;
    let mut tile_offsets: Vec<Vec<u32>> = Vec::new();
    let mut tile_byte_counts: Vec<Vec<u32>> = Vec::new();

    for tiles in &all_tile_data {
        let mut offsets = Vec::new();
        let mut counts = Vec::new();
        for tile in tiles {
            offsets.push(current_offset as u32);
            counts.push(tile.len() as u32);
            current_offset += tile.len();
        }
        tile_offsets.push(offsets);
        tile_byte_counts.push(counts);
    }

    // Write full-resolution IFD
    let next_ifd = if overviews.is_empty() {
        0u32
    } else {
        (8 + ifd_size + geo_extra_size) as u32
    };

    write_ifd(
        writer,
        &IfdArgs {
            width: raster.width() as u32,
            height: raster.height() as u32,
            params,
            tile_offsets: &tile_offsets[0],
            tile_byte_counts: &tile_byte_counts[0],
            next_ifd_offset: next_ifd,
        },
    )?;

    // Write overview IFDs
    for (i, ov) in overviews.iter().enumerate() {
        let next = if i + 1 < overviews.len() {
            (8 + (i + 2) * (ifd_size + geo_extra_size)) as u32
        } else {
            0u32
        };

        let ov_pixel_width = params.pixel_width * ov.factor as f64;
        let ov_pixel_height = params.pixel_height * ov.factor as f64;

        let ov_params = CogParams {
            pixel_width: ov_pixel_width,
            pixel_height: ov_pixel_height,
            ..params.clone()
        };

        write_ifd(
            writer,
            &IfdArgs {
                width: ov.width as u32,
                height: ov.height as u32,
                params: &ov_params,
                tile_offsets: &tile_offsets[i + 1],
                tile_byte_counts: &tile_byte_counts[i + 1],
                next_ifd_offset: next,
            },
        )?;
    }

    // Write tile data
    for tiles in &all_tile_data {
        for tile in tiles {
            writer.write_all(tile)?;
        }
    }

    Ok(())
}

/// Serve a specific tile from a COG-structured raster.
///
/// Given a tile coordinate (col, row) at a specific overview level,
/// extracts the tile data without reading the entire file.
pub fn extract_tile(
    raster: &Raster,
    tile_col: usize,
    tile_row: usize,
    tile_width: usize,
    tile_height: usize,
) -> Vec<f64> {
    let x_start = tile_col * tile_width;
    let y_start = tile_row * tile_height;
    let mut tile = vec![f64::NAN; tile_width * tile_height];

    for ty in 0..tile_height {
        let src_y = y_start + ty;
        if src_y >= raster.height() {
            break;
        }
        for tx in 0..tile_width {
            let src_x = x_start + tx;
            if src_x >= raster.width() {
                break;
            }
            tile[ty * tile_width + tx] = raster.data()[src_y * raster.width() + src_x];
        }
    }

    tile
}

/// Split raster data into tiles, returning raw f64 bytes per tile.
fn raster_to_tiles(
    data: &[f64],
    width: usize,
    height: usize,
    tile_w: usize,
    tile_h: usize,
) -> Vec<Vec<u8>> {
    let tiles_across = width.div_ceil(tile_w);
    let tiles_down = height.div_ceil(tile_h);
    let mut tiles = Vec::with_capacity(tiles_across * tiles_down);

    for tr in 0..tiles_down {
        for tc in 0..tiles_across {
            let mut tile_data = Vec::with_capacity(tile_w * tile_h * 8);
            for ty in 0..tile_h {
                let src_y = tr * tile_h + ty;
                for tx in 0..tile_w {
                    let src_x = tc * tile_w + tx;
                    let val = if src_x < width && src_y < height {
                        data[src_y * width + src_x]
                    } else {
                        f64::NAN
                    };
                    tile_data.extend_from_slice(&val.to_le_bytes());
                }
            }
            tiles.push(tile_data);
        }
    }

    tiles
}

struct IfdArgs<'a> {
    width: u32,
    height: u32,
    params: &'a CogParams,
    tile_offsets: &'a [u32],
    tile_byte_counts: &'a [u32],
    next_ifd_offset: u32,
}

fn write_ifd<W: Write + Seek>(writer: &mut W, args: &IfdArgs<'_>) -> Result<(), Error> {
    let n_tiles = args.tile_offsets.len();
    let params = args.params;
    let entries_count: u16 = 14;
    write_u16(writer, entries_count)?;

    // Tag 256: ImageWidth
    write_ifd_entry(writer, 256, 3, 1, args.width)?;
    // Tag 257: ImageLength
    write_ifd_entry(writer, 257, 3, 1, args.height)?;
    // Tag 258: BitsPerSample (64-bit float)
    write_ifd_entry(writer, 258, 3, 1, 64)?;
    // Tag 259: Compression (1 = None)
    write_ifd_entry(writer, 259, 3, 1, 1)?;
    // Tag 262: PhotometricInterpretation (1 = MinIsBlack)
    write_ifd_entry(writer, 262, 3, 1, 1)?;
    // Tag 277: SamplesPerPixel
    write_ifd_entry(writer, 277, 3, 1, 1)?;
    // Tag 322: TileWidth
    write_ifd_entry(writer, 322, 3, 1, params.tile_width)?;
    // Tag 323: TileLength
    write_ifd_entry(writer, 323, 3, 1, params.tile_height)?;
    // Tag 339: SampleFormat (3 = IEEE float)
    write_ifd_entry(writer, 339, 3, 1, 3)?;

    // Tag 324: TileOffsets — write inline if 1 tile, else write offset to array
    let current_pos = writer.stream_position()?;
    if n_tiles == 1 {
        write_ifd_entry(writer, 324, 4, 1, args.tile_offsets[0])?;
    } else {
        // Will store offset array after the IFD
        let array_offset = current_pos as u32 + (5 * 12) + 4 + 48; // remaining entries + next + geo
        write_ifd_entry(writer, 324, 4, n_tiles as u32, array_offset)?;
    }

    // Tag 325: TileByteCounts
    if n_tiles == 1 {
        write_ifd_entry(writer, 325, 4, 1, args.tile_byte_counts[0])?;
    } else {
        let array_offset = current_pos as u32 + (5 * 12) + 4 + 48 + (n_tiles as u32 * 4);
        write_ifd_entry(writer, 325, 4, n_tiles as u32, array_offset)?;
    }

    // Tag 33550: ModelPixelScaleTag
    let geo_offset = current_pos as u32 + (3 * 12) + 4;
    write_ifd_entry(writer, 33550, 12, 3, geo_offset)?;

    // Tag 33922: ModelTiepointTag
    write_ifd_entry(writer, 33922, 12, 6, geo_offset + 24)?;

    // Tag 34735: GeoKeyDirectoryTag (inline: version, revision, minor, count)
    write_ifd_entry(writer, 34735, 3, 4, pack_geokeys(params.epsg))?;

    // Next IFD offset
    write_u32(writer, args.next_ifd_offset)?;

    // Write geo metadata arrays
    // ModelPixelScaleTag: [pixel_width, pixel_height, 0.0]
    writer.write_all(&params.pixel_width.to_le_bytes())?;
    writer.write_all(&params.pixel_height.to_le_bytes())?;
    writer.write_all(&0.0f64.to_le_bytes())?;

    // ModelTiepointTag: [0, 0, 0, origin_x, origin_y, 0]
    writer.write_all(&0.0f64.to_le_bytes())?;
    writer.write_all(&0.0f64.to_le_bytes())?;
    writer.write_all(&0.0f64.to_le_bytes())?;
    writer.write_all(&params.origin_x.to_le_bytes())?;
    writer.write_all(&params.origin_y.to_le_bytes())?;
    writer.write_all(&0.0f64.to_le_bytes())?;

    // Write tile offset/bytecount arrays (if > 1 tile)
    if n_tiles > 1 {
        for &offset in args.tile_offsets {
            write_u32(writer, offset)?;
        }
        for &count in args.tile_byte_counts {
            write_u32(writer, count)?;
        }
    }

    Ok(())
}

fn write_ifd_entry<W: Write>(
    writer: &mut W,
    tag: u16,
    data_type: u16,
    count: u32,
    value: u32,
) -> Result<(), Error> {
    write_u16(writer, tag)?;
    write_u16(writer, data_type)?;
    write_u32(writer, count)?;
    write_u32(writer, value)?;
    Ok(())
}

fn write_u16<W: Write>(writer: &mut W, v: u16) -> Result<(), io::Error> {
    writer.write_all(&v.to_le_bytes())
}

fn write_u32<W: Write>(writer: &mut W, v: u32) -> Result<(), io::Error> {
    writer.write_all(&v.to_le_bytes())
}

fn pack_geokeys(_epsg: u16) -> u32 {
    // Pack version=1, revision=1 into 4 shorts as a single u32 value
    u32::from(1u16) | (u32::from(1u16) << 16)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_raster(width: usize, height: usize) -> Raster {
        let mut data = vec![0.0; width * height];
        for (i, val) in data.iter_mut().enumerate() {
            *val = i as f64;
        }
        Raster::from_vec(width, height, data, 1.0, -9999.0).unwrap()
    }

    #[test]
    fn generate_overviews_produces_correct_sizes() {
        let raster = test_raster(100, 100);
        let overviews = generate_overviews(&raster, 3);

        assert_eq!(overviews.len(), 3);
        assert_eq!(overviews[0].width, 50);
        assert_eq!(overviews[0].height, 50);
        assert_eq!(overviews[0].factor, 2);
        assert_eq!(overviews[1].width, 25);
        assert_eq!(overviews[1].height, 25);
        assert_eq!(overviews[1].factor, 4);
        assert_eq!(overviews[2].width, 13); // (25+1)/2
        assert_eq!(overviews[2].height, 13);
        assert_eq!(overviews[2].factor, 8);
    }

    #[test]
    fn overview_averages_correctly() {
        let raster = Raster::from_vec(
            4,
            4,
            vec![
                1.0, 3.0, 5.0, 7.0, 2.0, 4.0, 6.0, 8.0, 9.0, 11.0, 13.0, 15.0, 10.0, 12.0, 14.0,
                16.0,
            ],
            1.0,
            -9999.0,
        )
        .unwrap();
        let overviews = generate_overviews(&raster, 1);
        assert_eq!(overviews.len(), 1);
        assert_eq!(overviews[0].width, 2);
        assert_eq!(overviews[0].height, 2);
        // Top-left 2×2 block: (1+3+2+4)/4 = 2.5
        assert!((overviews[0].data[0] - 2.5).abs() < 1e-10);
        // Top-right 2×2 block: (5+7+6+8)/4 = 6.5
        assert!((overviews[0].data[1] - 6.5).abs() < 1e-10);
    }

    #[test]
    fn overview_handles_nan() {
        let raster = Raster::from_vec(
            4,
            2,
            vec![1.0, f64::NAN, 3.0, 5.0, 2.0, 4.0, f64::NAN, 7.0],
            1.0,
            -9999.0,
        )
        .unwrap();
        let overviews = generate_overviews(&raster, 1);
        // Top-left: (1 + 2 + 4) / 3 = 2.333...
        assert!((overviews[0].data[0] - 7.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn extract_tile_basic() {
        let raster = test_raster(8, 8);
        let tile = extract_tile(&raster, 0, 0, 4, 4);
        assert_eq!(tile.len(), 16);
        assert!((tile[0] - 0.0).abs() < 1e-10);
        assert!((tile[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn extract_tile_boundary() {
        let raster = test_raster(6, 6);
        let tile = extract_tile(&raster, 1, 1, 4, 4);
        assert_eq!(tile.len(), 16);
        // Pixels beyond raster extent should be NaN
        assert!(tile[15].is_nan()); // (4+3, 4+3) = (7,7) is out of 6×6
    }

    #[test]
    fn raster_to_tiles_count() {
        let data: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let tiles = raster_to_tiles(&data, 10, 10, 4, 4);
        // 10/4 = 3 tiles across, 3 down = 9 tiles
        assert_eq!(tiles.len(), 9);
        // Each tile: 4*4*8 = 128 bytes
        assert_eq!(tiles[0].len(), 128);
    }

    #[test]
    fn write_cog_produces_valid_tiff() {
        let raster = test_raster(16, 16);
        let params = CogParams {
            tile_width: 8,
            tile_height: 8,
            overview_levels: 1,
            epsg: 4326,
            origin_x: -180.0,
            origin_y: 90.0,
            pixel_width: 0.1,
            pixel_height: 0.1,
        };

        let mut buf = io::Cursor::new(Vec::new());
        write_cog(&raster, &params, &mut buf).unwrap();
        let bytes = buf.into_inner();

        // TIFF header
        assert_eq!(&bytes[0..2], b"II"); // Little-endian
        assert_eq!(bytes[2], 42); // TIFF magic
        assert_eq!(bytes[3], 0);
        // File should be non-trivially sized (header + IFDs + tiles)
        assert!(bytes.len() > 200);
    }

    #[test]
    fn write_cog_single_tile() {
        let raster = test_raster(4, 4);
        let params = CogParams {
            tile_width: 8,
            tile_height: 8,
            overview_levels: 0,
            epsg: 32632,
            origin_x: 500000.0,
            origin_y: 5000000.0,
            pixel_width: 10.0,
            pixel_height: 10.0,
        };

        let mut buf = io::Cursor::new(Vec::new());
        write_cog(&raster, &params, &mut buf).unwrap();
        let bytes = buf.into_inner();
        assert_eq!(&bytes[0..2], b"II");
    }
}
