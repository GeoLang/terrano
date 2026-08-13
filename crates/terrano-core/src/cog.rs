//! Cloud Optimized GeoTIFF (COG) support.
//!
//! Writes COGs with internal tiling, overview pyramids, and all IFDs at
//! the start of the file, and reads them back windowed: [`CogReader`] over
//! a [`RangeRead`] source fetches only the tiles a window touches, so a
//! `Range`-request http transport streams remote files without downloading
//! them. The read side covers real-world cogs: uncompressed or deflate
//! tiles, horizontal and floating-point predictors, integer and float
//! sample types (decoded to f64, nodata mapped to NaN). The writer produces
//! 64-bit float tiles, raw or deflate, and its output passes GDAL's
//! `validate_cloud_optimized_geotiff.py`.
//!
//! Multi-band files are pixel-interleaved (PlanarConfiguration 1) and go
//! through [`write_cog_bands`] and [`CogReader::read_window_bands`]. The
//! single-band [`write_cog`] and [`CogReader::read_window`] are the
//! one-band case of the same pipeline.
//!
//! COG files follow the standard described at <https://www.cogeo.org/>.

use crate::Error;
use crate::banded::BandedRaster;
use crate::geotiff::GeoTiffMetadata;
use crate::raster::Raster;
use std::io::{self, Write};

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
    /// Compress tiles with zlib deflate.
    pub deflate: bool,
    /// Value declared as GDAL_NODATA and substituted for NaN samples.
    ///
    /// `None` writes no nodata tag and leaves NaN in the file, which only
    /// readers that treat NaN as absent will understand.
    pub nodata: Option<f64>,
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
            deflate: false,
            nodata: Some(f64::NAN),
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
/// 2. Full-resolution IFD, then overview IFDs by decreasing resolution,
///    each followed by its geo arrays, GeoKey directory, and tile
///    offset/count arrays
/// 3. Tile data, smallest overview first and full resolution last
///
/// All IFD metadata sits at the start of the file, so a reader learns the
/// full layout from one small prefix read, and a viewer that only wants a
/// zoomed-out view stops reading after the overviews.
///
/// Overview IFDs carry NewSubfileType, which is what marks them a pyramid
/// rather than extra pages.
pub fn write_cog<W: Write>(
    raster: &Raster,
    params: &CogParams,
    writer: &mut W,
) -> Result<(), Error> {
    write_bands(&[raster], params, writer)
}

/// Write a [`BandedRaster`] as a pixel-interleaved multi-band COG.
///
/// Same layout as [`write_cog`], with SamplesPerPixel set to the band count
/// and PlanarConfiguration chunky. Overviews are block-averaged per band.
pub fn write_cog_bands<W: Write>(
    bands: &BandedRaster,
    params: &CogParams,
    writer: &mut W,
) -> Result<(), Error> {
    let refs: Vec<&Raster> = bands.bands().iter().collect();
    write_bands(&refs, params, writer)
}

fn write_bands<W: Write>(
    bands: &[&Raster],
    params: &CogParams,
    writer: &mut W,
) -> Result<(), Error> {
    let samples = bands.len();
    let first = bands
        .first()
        .ok_or_else(|| Error::InvalidInput("a cog needs at least one band".into()))?;
    let (width, height) = (first.width(), first.height());
    let per_band: Vec<Vec<Overview>> = bands
        .iter()
        .map(|b| generate_overviews(b, params.overview_levels))
        .collect();
    let overviews = &per_band[0];

    let fill = params.nodata.unwrap_or(f64::NAN);
    let mut all_tile_data: Vec<Vec<Vec<u8>>> = Vec::new();
    let planes: Vec<&[f64]> = bands.iter().map(|b| b.data()).collect();
    all_tile_data.push(bands_to_tiles(
        &planes,
        width,
        height,
        params.tile_width as usize,
        params.tile_height as usize,
        fill,
    ));
    for (i, ov) in overviews.iter().enumerate() {
        let planes: Vec<&[f64]> = per_band.iter().map(|b| b[i].data.as_slice()).collect();
        all_tile_data.push(bands_to_tiles(
            &planes,
            ov.width,
            ov.height,
            params.tile_width as usize,
            params.tile_height as usize,
            fill,
        ));
    }
    if params.deflate {
        for tiles in &mut all_tile_data {
            for tile in tiles.iter_mut() {
                let mut enc =
                    flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
                enc.write_all(tile)?;
                *tile = enc.finish()?;
            }
        }
    }

    // per-ifd byte footprint decides every offset up front
    let nodata_text = params.nodata.map(nodata_text);
    let layouts: Vec<IfdLayout> = all_tile_data
        .iter()
        .enumerate()
        .map(|(level, tiles)| IfdLayout {
            samples,
            n_tiles: tiles.len(),
            overview: level > 0,
            nodata: nodata_text.clone(),
        })
        .collect();
    let mut ifd_starts = Vec::with_capacity(layouts.len());
    let mut at = 8usize;
    for layout in &layouts {
        ifd_starts.push(at);
        at += layout.total();
    }

    // a viewer reading coarse-to-fine wants the small levels first, and the
    // cog spec requires that order
    let mut tile_offsets: Vec<Vec<u32>> = vec![Vec::new(); all_tile_data.len()];
    let mut tile_byte_counts: Vec<Vec<u32>> = vec![Vec::new(); all_tile_data.len()];
    for level in (0..all_tile_data.len()).rev() {
        for tile in &all_tile_data[level] {
            tile_offsets[level].push(at as u32);
            tile_byte_counts[level].push(tile.len() as u32);
            at += tile.len();
        }
    }
    if at > u32::MAX as usize {
        return Err(Error::InvalidInput(
            "image too large for a classic tiff, offsets are 32-bit".into(),
        ));
    }

    writer.write_all(b"II")?;
    write_u16(writer, 42)?;
    write_u32(writer, 8)?;

    for i in 0..all_tile_data.len() {
        let (level_width, level_height, pixel_width, pixel_height) = if i == 0 {
            (
                width as u32,
                height as u32,
                params.pixel_width,
                params.pixel_height,
            )
        } else {
            let ov = &overviews[i - 1];
            (
                ov.width as u32,
                ov.height as u32,
                params.pixel_width * ov.factor as f64,
                params.pixel_height * ov.factor as f64,
            )
        };
        let next = if i + 1 < ifd_starts.len() {
            ifd_starts[i + 1] as u32
        } else {
            0
        };
        let ifd = build_ifd(&IfdArgs {
            base: ifd_starts[i] as u32,
            width: level_width,
            height: level_height,
            layout: &layouts[i],
            pixel_width,
            pixel_height,
            params,
            tile_offsets: &tile_offsets[i],
            tile_byte_counts: &tile_byte_counts[i],
            next_ifd_offset: next,
        });
        debug_assert_eq!(ifd.len(), layouts[i].total());
        writer.write_all(&ifd)?;
    }

    for tiles in all_tile_data.iter().rev() {
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
#[cfg(test)]
fn raster_to_tiles(
    data: &[f64],
    width: usize,
    height: usize,
    tile_w: usize,
    tile_h: usize,
) -> Vec<Vec<u8>> {
    bands_to_tiles(&[data], width, height, tile_w, tile_h, f64::NAN)
}

/// Split per-band planes into pixel-interleaved tiles of raw f64 bytes.
///
/// `nodata` fills the padding past the image edge and replaces NaN samples,
/// so the value the file declares absent is the one actually stored.
fn bands_to_tiles(
    planes: &[&[f64]],
    width: usize,
    height: usize,
    tile_w: usize,
    tile_h: usize,
    nodata: f64,
) -> Vec<Vec<u8>> {
    let tiles_across = width.div_ceil(tile_w);
    let tiles_down = height.div_ceil(tile_h);
    let mut tiles = Vec::with_capacity(tiles_across * tiles_down);

    for tr in 0..tiles_down {
        for tc in 0..tiles_across {
            let mut tile_data = Vec::with_capacity(tile_w * tile_h * planes.len() * 8);
            for ty in 0..tile_h {
                let src_y = tr * tile_h + ty;
                for tx in 0..tile_w {
                    let src_x = tc * tile_w + tx;
                    let inside = src_x < width && src_y < height;
                    for plane in planes {
                        let mut val = if inside {
                            plane[src_y * width + src_x]
                        } else {
                            nodata
                        };
                        if val.is_nan() {
                            val = nodata;
                        }
                        tile_data.extend_from_slice(&val.to_le_bytes());
                    }
                }
            }
            tiles.push(tile_data);
        }
    }

    tiles
}

/// entries every level carries: 256, 257, 258, 259, 262, 277, 322, 323,
/// 324, 325, 339, 33550, 33922, 34735
const BASE_IFD_ENTRIES: usize = 14;
/// 4-short header plus GTModelType, GTRasterType and the crs key
const GEO_KEY_SHORTS: usize = 4 + 3 * 4;
// pixel scale (24) + tiepoint (48) + geokey directory
const GEO_BYTES: usize = 24 + 48 + GEO_KEY_SHORTS * 2;
/// an ascii value of four bytes or fewer sits in the entry value field
const INLINE_BYTES: usize = 4;

/// BitsPerSample and SampleFormat carry one short per sample, and three or
/// more no longer fit in the entry value field
fn sample_array_bytes(samples: usize) -> usize {
    if samples > 2 { samples * 4 } else { 0 }
}

/// GDAL_NODATA is ascii, and GDAL spells a non-finite value in lowercase
fn nodata_text(value: f64) -> String {
    if value.is_nan() {
        "nan".into()
    } else if value.is_infinite() {
        if value > 0.0 {
            "inf".into()
        } else {
            "-inf".into()
        }
    } else {
        format!("{value}")
    }
}

/// Byte layout of one level's IFD.
///
/// The size used to place the IFD and the bytes later written for it both
/// come from here, so they cannot drift apart as optional entries change.
struct IfdLayout {
    samples: usize,
    n_tiles: usize,
    overview: bool,
    nodata: Option<String>,
}

impl IfdLayout {
    /// bands past the first are unspecified data, and libtiff wants them
    /// declared: colour channels plus extra samples has to reach
    /// SamplesPerPixel
    fn extra_samples(&self) -> usize {
        self.samples - 1
    }

    fn extra_sample_bytes(&self) -> usize {
        if self.extra_samples() > 2 {
            self.extra_samples() * 2
        } else {
            0
        }
    }

    fn entries(&self) -> usize {
        BASE_IFD_ENTRIES
            + usize::from(self.overview)
            + 2 * usize::from(self.samples > 1)
            + usize::from(self.nodata.is_some())
    }

    fn directory_bytes(&self) -> usize {
        2 + self.entries() * 12 + 4
    }

    /// ascii bytes including the terminating nul, zero when it fits inline
    fn nodata_bytes(&self) -> usize {
        match &self.nodata {
            Some(text) if text.len() + 1 > INLINE_BYTES => text.len() + 1,
            _ => 0,
        }
    }

    fn tile_array_bytes(&self) -> usize {
        if self.n_tiles > 1 {
            self.n_tiles * 8
        } else {
            0
        }
    }

    fn total(&self) -> usize {
        self.directory_bytes()
            + GEO_BYTES
            + sample_array_bytes(self.samples)
            + self.extra_sample_bytes()
            + self.nodata_bytes()
            + self.tile_array_bytes()
    }
}

struct IfdArgs<'a> {
    base: u32,
    width: u32,
    height: u32,
    layout: &'a IfdLayout,
    pixel_width: f64,
    pixel_height: f64,
    params: &'a CogParams,
    tile_offsets: &'a [u32],
    tile_byte_counts: &'a [u32],
    next_ifd_offset: u32,
}

fn build_ifd(args: &IfdArgs<'_>) -> Vec<u8> {
    let layout = args.layout;
    let n_tiles = args.tile_offsets.len();
    let samples = layout.samples;
    let aux = args.base + layout.directory_bytes() as u32;
    let scale_off = aux;
    let tiepoint_off = aux + 24;
    let geokeys_off = aux + 72;
    let bits_off = aux + GEO_BYTES as u32;
    let formats_off = bits_off + if samples > 2 { samples as u32 * 2 } else { 0 };
    let extra_off = aux + GEO_BYTES as u32 + sample_array_bytes(samples) as u32;
    let nodata_off = extra_off + layout.extra_sample_bytes() as u32;
    let arrays_off = nodata_off + layout.nodata_bytes() as u32;

    // one short per sample: inline while two fit in the value field
    let packed = |value: u32| -> u32 {
        if samples == 1 {
            value
        } else {
            value | (value << 16)
        }
    };

    let mut buf = Vec::with_capacity(layout.total());
    push_u16(&mut buf, layout.entries() as u16);
    // NewSubfileType: reduced-resolution, the bit that makes readers treat
    // these ifds as a pyramid instead of separate pages
    if layout.overview {
        push_entry(&mut buf, 254, 4, 1, 1);
    }
    push_entry(&mut buf, 256, 3, 1, args.width);
    push_entry(&mut buf, 257, 3, 1, args.height);
    let bits_value = if samples > 2 { bits_off } else { packed(64) };
    push_entry(&mut buf, 258, 3, samples as u32, bits_value);
    push_entry(&mut buf, 259, 3, 1, if args.params.deflate { 8 } else { 1 });
    push_entry(&mut buf, 262, 3, 1, 1);
    push_entry(&mut buf, 277, 3, 1, samples as u32);
    if samples > 1 {
        push_entry(&mut buf, 284, 3, 1, 1);
    }
    push_entry(&mut buf, 322, 3, 1, args.params.tile_width);
    push_entry(&mut buf, 323, 3, 1, args.params.tile_height);
    if n_tiles == 1 {
        push_entry(&mut buf, 324, 4, 1, args.tile_offsets[0]);
        push_entry(&mut buf, 325, 4, 1, args.tile_byte_counts[0]);
    } else {
        push_entry(&mut buf, 324, 4, n_tiles as u32, arrays_off);
        push_entry(
            &mut buf,
            325,
            4,
            n_tiles as u32,
            arrays_off + n_tiles as u32 * 4,
        );
    }
    if samples > 1 {
        let value = if layout.extra_sample_bytes() > 0 {
            extra_off
        } else {
            0
        };
        push_entry(&mut buf, 338, 3, layout.extra_samples() as u32, value);
    }
    let format_value = if samples > 2 { formats_off } else { packed(3) };
    push_entry(&mut buf, 339, 3, samples as u32, format_value);
    push_entry(&mut buf, 33550, 12, 3, scale_off);
    push_entry(&mut buf, 33922, 12, 6, tiepoint_off);
    push_entry(&mut buf, 34735, 3, GEO_KEY_SHORTS as u32, geokeys_off);
    if let Some(text) = &layout.nodata {
        let count = text.len() as u32 + 1;
        let value = if layout.nodata_bytes() == 0 {
            let mut inline = [0u8; INLINE_BYTES];
            inline[..text.len()].copy_from_slice(text.as_bytes());
            u32::from_le_bytes(inline)
        } else {
            nodata_off
        };
        push_entry(&mut buf, 42113, 2, count, value);
    }
    push_u32(&mut buf, args.next_ifd_offset);

    buf.extend_from_slice(&args.pixel_width.to_le_bytes());
    buf.extend_from_slice(&args.pixel_height.to_le_bytes());
    buf.extend_from_slice(&0.0f64.to_le_bytes());

    for v in [
        0.0,
        0.0,
        0.0,
        args.params.origin_x,
        args.params.origin_y,
        0.0,
    ] {
        buf.extend_from_slice(&v.to_le_bytes());
    }

    // GeoKeyDirectory: header, then keys sorted by id
    let geographic = (4000..=4999).contains(&args.params.epsg);
    let model_type: u16 = if geographic { 2 } else { 1 };
    let crs_key: u16 = if geographic { 2048 } else { 3072 };
    for v in [
        1u16,
        1,
        0,
        3,
        1024,
        0,
        1,
        model_type,
        1025,
        0,
        1,
        1,
        crs_key,
        0,
        1,
        args.params.epsg,
    ] {
        push_u16(&mut buf, v);
    }

    if samples > 2 {
        for _ in 0..samples {
            push_u16(&mut buf, 64);
        }
        for _ in 0..samples {
            push_u16(&mut buf, 3);
        }
    }

    if layout.extra_sample_bytes() > 0 {
        for _ in 0..layout.extra_samples() {
            push_u16(&mut buf, 0);
        }
    }

    if layout.nodata_bytes() > 0 {
        let text = layout.nodata.as_ref().expect("nodata bytes imply a value");
        buf.extend_from_slice(text.as_bytes());
        buf.push(0);
    }

    if n_tiles > 1 {
        for &offset in args.tile_offsets {
            push_u32(&mut buf, offset);
        }
        for &count in args.tile_byte_counts {
            push_u32(&mut buf, count);
        }
    }

    buf
}

fn push_u16(buf: &mut Vec<u8>, v: u16) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn push_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn push_entry(buf: &mut Vec<u8>, tag: u16, data_type: u16, count: u32, value: u32) {
    push_u16(buf, tag);
    push_u16(buf, data_type);
    push_u32(buf, count);
    push_u32(buf, value);
}

fn write_u16<W: Write>(writer: &mut W, v: u16) -> Result<(), io::Error> {
    writer.write_all(&v.to_le_bytes())
}

fn write_u32<W: Write>(writer: &mut W, v: u32) -> Result<(), io::Error> {
    writer.write_all(&v.to_le_bytes())
}

/// byte-range access to a COG, the seam between the reader and its
/// transport. implement it over an http client with `Range` requests to
/// stream remote files, the reader never asks for more than it needs
pub trait RangeRead {
    fn read_range(&mut self, offset: u64, len: u64) -> Result<Vec<u8>, Error>;

    /// several ranges in one call, results in request order. transports
    /// that can multiplex override this, the default reads sequentially
    fn read_ranges(&mut self, ranges: &[(u64, u64)]) -> Result<Vec<Vec<u8>>, Error> {
        ranges
            .iter()
            .map(|&(offset, len)| self.read_range(offset, len))
            .collect()
    }
}

impl RangeRead for std::fs::File {
    fn read_range(&mut self, offset: u64, len: u64) -> Result<Vec<u8>, Error> {
        use std::io::{Read, Seek, SeekFrom};
        self.seek(SeekFrom::Start(offset))?;
        let mut buf = vec![0u8; len as usize];
        self.read_exact(&mut buf)?;
        Ok(buf)
    }
}

impl RangeRead for &[u8] {
    fn read_range(&mut self, offset: u64, len: u64) -> Result<Vec<u8>, Error> {
        let start = offset as usize;
        let end = start + len as usize;
        if end > self.len() {
            return Err(Error::Format(format!(
                "range {start}..{end} beyond end of data ({} bytes)",
                self.len()
            )));
        }
        Ok(self[start..end].to_vec())
    }
}

/// one resolution level of a cog: the full image or an overview
#[derive(Debug, Clone)]
pub struct CogLevel {
    pub width: usize,
    pub height: usize,
    pub tile_width: usize,
    pub tile_height: usize,
    pub pixel_width: f64,
    pub pixel_height: f64,
    /// samples per pixel, one per band
    pub samples: usize,
    tile_offsets: Vec<u64>,
    tile_byte_counts: Vec<u64>,
    compression: u16,
    predictor: u16,
    bits: u32,
    format: u32,
    nodata: Option<f64>,
}

impl CogLevel {
    fn tiles_across(&self) -> usize {
        self.width.div_ceil(self.tile_width)
    }
}

/// windowed cog reader over a [`RangeRead`] source.
///
/// `open` learns the layout from small header reads, `read_window` fetches
/// only the tiles a window touches. multi-band files read through
/// [`CogReader::read_window_bands`]
pub struct CogReader<R: RangeRead> {
    source: R,
    levels: Vec<CogLevel>,
    meta: GeoTiffMetadata,
}

impl<R: RangeRead> CogReader<R> {
    pub fn open(mut source: R) -> Result<Self, Error> {
        let header = source.read_range(0, 8)?;
        if &header[0..2] != b"II" || u16::from_le_bytes([header[2], header[3]]) != 42 {
            return Err(Error::Format("not a little-endian tiff".into()));
        }
        let mut next = u64::from(u32::from_le_bytes([
            header[4], header[5], header[6], header[7],
        ]));
        let mut levels = Vec::new();
        let mut meta: Option<GeoTiffMetadata> = None;
        while next != 0 {
            let count = u16::from_le_bytes(source.read_range(next, 2)?.try_into().unwrap()) as u64;
            let block = source.read_range(next + 2, count * 12 + 4)?;
            let (level, level_meta, following) = parse_ifd(&mut source, &block, count as usize)?;
            levels.push(level);
            if meta.is_none() {
                meta = Some(level_meta);
            }
            next = following;
        }
        if levels.is_empty() {
            return Err(Error::Format("tiff has no ifd".into()));
        }
        Ok(CogReader {
            source,
            levels,
            meta: meta.expect("meta set with first level"),
        })
    }

    pub fn levels(&self) -> &[CogLevel] {
        &self.levels
    }

    pub fn meta(&self) -> &GeoTiffMetadata {
        &self.meta
    }

    /// coarsest level still at least as fine as the target pixel size,
    /// level 0 when the target outresolves the file
    pub fn select_level(&self, target_pixel_width: f64) -> usize {
        let mut best = 0;
        for (i, level) in self.levels.iter().enumerate() {
            if level.pixel_width <= target_pixel_width * (1.0 + 1e-9)
                && level.pixel_width > self.levels[best].pixel_width
            {
                best = i;
            }
        }
        best
    }

    /// read a pixel window from one level, fetching only intersecting
    /// tiles. pixels outside the image come back NaN
    pub fn read_window(
        &mut self,
        level: usize,
        col0: usize,
        row0: usize,
        cols: usize,
        rows: usize,
    ) -> Result<Raster, Error> {
        let samples = self.level_at(level)?.samples;
        if samples != 1 {
            return Err(Error::Format(format!(
                "{samples} samples per pixel, read multi-band cogs with read_window_bands"
            )));
        }
        let bands = self.read_window_bands(level, col0, row0, cols, rows)?;
        Ok(bands.into_bands().remove(0))
    }

    /// read a pixel window from one level as one [`Raster`] per band,
    /// decoding each touched tile once. single-band files come back with
    /// one band. pixels outside the image come back NaN
    pub fn read_window_bands(
        &mut self,
        level: usize,
        col0: usize,
        row0: usize,
        cols: usize,
        rows: usize,
    ) -> Result<BandedRaster, Error> {
        let l = self.level_at(level)?.clone();
        let mut out = vec![vec![f64::NAN; cols * rows]; l.samples];
        if col0 < l.width && row0 < l.height && cols > 0 && rows > 0 {
            let last_col = (col0 + cols - 1).min(l.width - 1);
            let last_row = (row0 + rows - 1).min(l.height - 1);
            let mut tiles = Vec::new();
            let mut ranges = Vec::new();
            for tr in row0 / l.tile_height..=last_row / l.tile_height {
                for tc in col0 / l.tile_width..=last_col / l.tile_width {
                    let idx = tr * l.tiles_across() + tc;
                    tiles.push((tc, tr));
                    ranges.push((l.tile_offsets[idx], l.tile_byte_counts[idx]));
                }
            }
            // one call for every touched tile, so a multiplexing
            // transport fetches them concurrently
            let fetched = self.source.read_ranges(&ranges)?;
            for ((tc, tr), bytes) in tiles.into_iter().zip(fetched) {
                let values = decode_tile(&l, bytes)?;
                for (band, plane) in out.iter_mut().enumerate() {
                    copy_tile(&values, &l, band, tc, tr, col0, row0, cols, rows, plane);
                }
            }
        }
        let bands = out
            .into_iter()
            .map(|plane| Raster::from_vec(cols, rows, plane, l.pixel_width, f64::NAN))
            .collect::<Result<Vec<_>, _>>()?;
        BandedRaster::new(bands)
    }

    fn level_at(&self, level: usize) -> Result<&CogLevel, Error> {
        self.levels
            .get(level)
            .ok_or_else(|| Error::Format(format!("no overview level {level}")))
    }
}

/// raw tile bytes to f64 samples: decompress, undo the predictor, widen
/// the sample type, and map the declared nodata to NaN. multi-band samples
/// stay pixel-interleaved
fn decode_tile(l: &CogLevel, bytes: Vec<u8>) -> Result<Vec<f64>, Error> {
    let sample = (l.bits / 8) as usize;
    let expected = l.tile_width * l.tile_height * l.samples * sample;
    let mut raw = match l.compression {
        1 => bytes,
        8 => {
            use std::io::Read;
            let mut out = Vec::with_capacity(expected);
            flate2::read::ZlibDecoder::new(bytes.as_slice()).read_to_end(&mut out)?;
            out
        }
        other => return Err(Error::Format(format!("compression {other} unsupported"))),
    };
    if raw.len() != expected {
        return Err(Error::Format(format!(
            "tile decoded to {} bytes, expected {expected}",
            raw.len()
        )));
    }
    match l.predictor {
        1 => {}
        2 => undo_horizontal(&mut raw, l.tile_width * l.samples, l.samples, l.bits),
        3 => raw = undo_fp(&raw, l.tile_width * l.samples, sample),
        other => return Err(Error::Format(format!("predictor {other} unsupported"))),
    }
    let mut vals = widen(&raw, l.bits, l.format);
    if let Some(nd) = l.nodata {
        for v in &mut vals {
            if *v == nd || (nd.is_nan() && v.is_nan()) {
                *v = f64::NAN;
            }
        }
    }
    Ok(vals)
}

/// undo per-row horizontal differencing of integer samples (predictor 2).
/// `count` is samples per row across all bands, each difference taken
/// against the same band of the previous pixel
fn undo_horizontal(raw: &mut [u8], count: usize, spp: usize, bits: u32) {
    let sample = (bits / 8) as usize;
    for row in raw.chunks_exact_mut(count * sample) {
        match bits {
            8 => {
                for i in spp..count {
                    row[i] = row[i].wrapping_add(row[i - spp]);
                }
            }
            16 => {
                for i in spp..count {
                    let p = (i - spp) * 2;
                    let prev = u16::from_le_bytes(row[p..p + 2].try_into().unwrap());
                    let cur = u16::from_le_bytes(row[i * 2..i * 2 + 2].try_into().unwrap());
                    row[i * 2..i * 2 + 2].copy_from_slice(&prev.wrapping_add(cur).to_le_bytes());
                }
            }
            _ => {
                for i in spp..count {
                    let p = (i - spp) * 4;
                    let prev = u32::from_le_bytes(row[p..p + 4].try_into().unwrap());
                    let cur = u32::from_le_bytes(row[i * 4..i * 4 + 4].try_into().unwrap());
                    row[i * 4..i * 4 + 4].copy_from_slice(&prev.wrapping_add(cur).to_le_bytes());
                }
            }
        }
    }
}

/// undo the floating-point predictor (3): per row, cumulative byte sum,
/// then reassemble from byte planes stored most significant first
fn undo_fp(raw: &[u8], count: usize, sample: usize) -> Vec<u8> {
    let mut out = vec![0u8; raw.len()];
    let stride = count * sample;
    for (r, row) in raw.chunks_exact(stride).enumerate() {
        let mut planes = row.to_vec();
        for i in 1..stride {
            planes[i] = planes[i].wrapping_add(planes[i - 1]);
        }
        for i in 0..count {
            for k in 0..sample {
                // byte k of the big-endian value sits in plane k
                out[r * stride + i * sample + (sample - 1 - k)] = planes[k * count + i];
            }
        }
    }
    out
}

/// widen little-endian samples of the declared format to f64
fn widen(raw: &[u8], bits: u32, format: u32) -> Vec<f64> {
    let le2 = |c: &[u8]| [c[0], c[1]];
    let le4 = |c: &[u8]| [c[0], c[1], c[2], c[3]];
    match (format, bits) {
        (1, 8) => raw.iter().map(|&b| f64::from(b)).collect(),
        (2, 8) => raw.iter().map(|&b| f64::from(b as i8)).collect(),
        (1, 16) => raw
            .chunks_exact(2)
            .map(|c| f64::from(u16::from_le_bytes(le2(c))))
            .collect(),
        (2, 16) => raw
            .chunks_exact(2)
            .map(|c| f64::from(i16::from_le_bytes(le2(c))))
            .collect(),
        (1, 32) => raw
            .chunks_exact(4)
            .map(|c| f64::from(u32::from_le_bytes(le4(c))))
            .collect(),
        (2, 32) => raw
            .chunks_exact(4)
            .map(|c| f64::from(i32::from_le_bytes(le4(c))))
            .collect(),
        (3, 32) => raw
            .chunks_exact(4)
            .map(|c| f64::from(f32::from_le_bytes(le4(c))))
            .collect(),
        _ => raw
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().expect("8-byte samples")))
            .collect(),
    }
}

/// copy one band of a decoded tile into the window buffer
#[allow(clippy::too_many_arguments)]
fn copy_tile(
    values: &[f64],
    l: &CogLevel,
    band: usize,
    tc: usize,
    tr: usize,
    col0: usize,
    row0: usize,
    cols: usize,
    rows: usize,
    out: &mut [f64],
) {
    let tile_x = tc * l.tile_width;
    let tile_y = tr * l.tile_height;
    for ty in 0..l.tile_height {
        let img_y = tile_y + ty;
        if img_y < row0 || img_y > row0 + rows - 1 || img_y >= l.height {
            continue;
        }
        for tx in 0..l.tile_width {
            let img_x = tile_x + tx;
            if img_x < col0 || img_x > col0 + cols - 1 || img_x >= l.width {
                continue;
            }
            out[(img_y - row0) * cols + (img_x - col0)] =
                values[(ty * l.tile_width + tx) * l.samples + band];
        }
    }
}

/// a SHORT-valued ifd entry that may hold one value per sample: up to two
/// sit in the entry itself, more live at an offset
#[derive(Default)]
struct ShortList {
    count: u32,
    inline: [u8; 4],
}

impl ShortList {
    fn at(block: &[u8], entry: usize, count: u32) -> Self {
        let mut inline = [0u8; 4];
        inline.copy_from_slice(&block[entry + 8..entry + 12]);
        Self { count, inline }
    }

    fn resolve<R: RangeRead>(&self, source: &mut R) -> Result<Vec<u16>, Error> {
        if self.count == 0 {
            return Ok(Vec::new());
        }
        let bytes = if self.count <= 2 {
            self.inline.to_vec()
        } else {
            let at = u32::from_le_bytes(self.inline);
            source.read_range(u64::from(at), u64::from(self.count) * 2)?
        };
        Ok(bytes
            .chunks_exact(2)
            .take(self.count as usize)
            .map(|c| u16::from_le_bytes(c.try_into().unwrap()))
            .collect())
    }
}

fn parse_ifd<R: RangeRead>(
    source: &mut R,
    block: &[u8],
    count: usize,
) -> Result<(CogLevel, GeoTiffMetadata, u64), Error> {
    let entry_u16 = |off: usize| u16::from_le_bytes(block[off..off + 2].try_into().unwrap());
    let entry_u32 = |off: usize| u32::from_le_bytes(block[off..off + 4].try_into().unwrap());

    let (mut width, mut height, mut tile_w, mut tile_h) = (0usize, 0usize, 0usize, 0usize);
    let (mut compression, mut bits, mut format, mut samples) = (1u32, 64u32, 3u32, 1u32);
    let (mut predictor, mut planar) = (1u32, 1u32);
    let mut bits_entry = ShortList::default();
    let mut format_entry = ShortList::default();
    let (mut offsets_at, mut counts_at, mut n_tiles) = (0u32, 0u32, 0usize);
    let (mut scale_at, mut tiepoint_at) = (0u32, 0u32);
    let (mut geokeys_at, mut geokeys_count) = (0u32, 0u32);
    let (mut nodata_at, mut nodata_count) = (0u32, 0u32);
    let mut nodata_inline = [0u8; 4];

    for i in 0..count {
        let e = i * 12;
        let tag = entry_u16(e);
        let entry_count = entry_u32(e + 4);
        let value = entry_u32(e + 8);
        match tag {
            256 => width = value as usize,
            257 => height = value as usize,
            258 => bits_entry = ShortList::at(block, e, entry_count),
            259 => compression = value & 0xffff,
            277 => samples = value & 0xffff,
            284 => planar = value & 0xffff,
            317 => predictor = value & 0xffff,
            322 => tile_w = value as usize,
            323 => tile_h = value as usize,
            324 => {
                offsets_at = value;
                n_tiles = entry_count as usize;
            }
            325 => counts_at = value,
            339 => format_entry = ShortList::at(block, e, entry_count),
            33550 => scale_at = value,
            33922 => tiepoint_at = value,
            34735 => {
                geokeys_at = value;
                geokeys_count = entry_count;
            }
            // GDAL_NODATA, an ascii number
            42113 => {
                nodata_at = value;
                nodata_count = entry_count;
                nodata_inline.copy_from_slice(&block[e + 8..e + 12]);
            }
            _ => {}
        }
    }

    let bits_list = bits_entry.resolve(source)?;
    let format_list = format_entry.resolve(source)?;
    if let Some(&first) = bits_list.first() {
        if bits_list.iter().any(|&b| b != first) {
            return Err(Error::Format(
                "BitsPerSample differs across samples, every band must share one bit depth".into(),
            ));
        }
        bits = u32::from(first);
    }
    if let Some(&first) = format_list.first() {
        if format_list.iter().any(|&f| f != first) {
            return Err(Error::Format(
                "SampleFormat differs across samples, every band must share one sample format"
                    .into(),
            ));
        }
        format = u32::from(first);
    }

    if width == 0 || height == 0 {
        return Err(Error::Format("ifd missing width/height".into()));
    }
    if tile_w == 0 || tile_h == 0 {
        return Err(Error::Format(
            "not a tiled tiff, stripped files are read whole via read_geotiff".into(),
        ));
    }
    if compression != 1 && compression != 8 {
        return Err(Error::Format(format!(
            "compression {compression} not supported, only uncompressed (1) and deflate (8)"
        )));
    }
    if samples == 0 {
        return Err(Error::Format("SamplesPerPixel is 0".into()));
    }
    if samples > 1 && planar != 1 {
        return Err(Error::Format(
            "only pixel-interleaved (PlanarConfiguration 1) cogs supported".into(),
        ));
    }
    let known = matches!(
        (format, bits),
        (1, 8) | (1, 16) | (1, 32) | (2, 8) | (2, 16) | (2, 32) | (3, 32) | (3, 64)
    );
    if !known {
        return Err(Error::Format(format!(
            "unsupported sample layout: {bits} bits format {format}"
        )));
    }
    let predictor_ok = match predictor {
        1 => true,
        2 => format != 3,
        3 => format == 3,
        _ => false,
    };
    if !predictor_ok {
        return Err(Error::Format(format!(
            "predictor {predictor} not supported for sample format {format}"
        )));
    }

    let nodata = if nodata_count > 0 {
        let bytes = if nodata_count <= 4 {
            nodata_inline[..nodata_count as usize].to_vec()
        } else {
            source.read_range(u64::from(nodata_at), u64::from(nodata_count))?
        };
        let text: String = bytes
            .iter()
            .take_while(|&&b| b != 0)
            .map(|&b| b as char)
            .collect();
        text.trim().parse::<f64>().ok()
    } else {
        None
    };

    let read_u32s = |source: &mut R, at: u32, n: usize| -> Result<Vec<u64>, Error> {
        let bytes = source.read_range(u64::from(at), n as u64 * 4)?;
        Ok(bytes
            .chunks_exact(4)
            .map(|c| u64::from(u32::from_le_bytes(c.try_into().unwrap())))
            .collect())
    };
    let (tile_offsets, tile_byte_counts) = if n_tiles <= 1 {
        (vec![u64::from(offsets_at)], vec![u64::from(counts_at)])
    } else {
        (
            read_u32s(source, offsets_at, n_tiles)?,
            read_u32s(source, counts_at, n_tiles)?,
        )
    };
    let expected_tiles = width.div_ceil(tile_w) * height.div_ceil(tile_h);
    if tile_offsets.len() != expected_tiles {
        return Err(Error::Format(format!(
            "tile count {} does not match image geometry ({expected_tiles})",
            tile_offsets.len()
        )));
    }

    let f64_at =
        |bytes: &[u8], off: usize| f64::from_le_bytes(bytes[off..off + 8].try_into().unwrap());
    let scale = source.read_range(u64::from(scale_at), 24)?;
    let tiepoint = source.read_range(u64::from(tiepoint_at), 48)?;
    let pixel_width = f64_at(&scale, 0);
    let pixel_height = f64_at(&scale, 8);

    let mut epsg = 0u16;
    if geokeys_count >= 8 {
        let keys = source.read_range(u64::from(geokeys_at), u64::from(geokeys_count) * 2)?;
        let key_u16 = |off: usize| u16::from_le_bytes(keys[off..off + 2].try_into().unwrap());
        let num_keys = key_u16(6) as usize;
        for k in 0..num_keys {
            let base = 8 + k * 8;
            if base + 8 > keys.len() {
                break;
            }
            let key_id = key_u16(base);
            if key_id == 2048 || key_id == 3072 {
                epsg = key_u16(base + 6);
            }
        }
    }

    let meta = GeoTiffMetadata {
        origin_x: f64_at(&tiepoint, 24),
        origin_y: f64_at(&tiepoint, 32),
        pixel_width,
        pixel_height,
        epsg,
    };
    let next = u64::from(u32::from_le_bytes(
        block[count * 12..count * 12 + 4].try_into().unwrap(),
    ));
    Ok((
        CogLevel {
            width,
            height,
            tile_width: tile_w,
            tile_height: tile_h,
            pixel_width,
            pixel_height,
            samples: samples as usize,
            tile_offsets,
            tile_byte_counts,
            compression: compression as u16,
            predictor: predictor as u16,
            bits,
            format,
            nodata,
        },
        meta,
        next,
    ))
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
            deflate: false,
            ..CogParams::default()
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
            deflate: false,
            ..CogParams::default()
        };

        let mut buf = io::Cursor::new(Vec::new());
        write_cog(&raster, &params, &mut buf).unwrap();
        let bytes = buf.into_inner();
        assert_eq!(&bytes[0..2], b"II");
    }

    struct RecordingReader<'a> {
        data: &'a [u8],
        fetched: usize,
    }

    impl RangeRead for RecordingReader<'_> {
        fn read_range(&mut self, offset: u64, len: u64) -> Result<Vec<u8>, Error> {
            self.fetched += len as usize;
            let mut slice = self.data;
            slice.read_range(offset, len)
        }
    }

    fn write_test_cog(width: usize, height: usize, levels: u32, epsg: u16) -> (Raster, Vec<u8>) {
        let raster = test_raster(width, height);
        let params = CogParams {
            tile_width: 16,
            tile_height: 16,
            overview_levels: levels,
            epsg,
            origin_x: 10.0,
            origin_y: 50.0,
            pixel_width: 0.5,
            pixel_height: 0.5,
            deflate: false,
            ..CogParams::default()
        };
        let mut buf = io::Cursor::new(Vec::new());
        write_cog(&raster, &params, &mut buf).unwrap();
        (raster, buf.into_inner())
    }

    #[test]
    fn read_window_roundtrips_full_resolution() {
        let (raster, bytes) = write_test_cog(100, 80, 2, 4326);
        let mut reader = CogReader::open(bytes.as_slice()).unwrap();
        let window = reader.read_window(0, 7, 5, 40, 30).unwrap();
        for row in 0..30 {
            for col in 0..40 {
                let expected = raster.data()[(row + 5) * 100 + (col + 7)];
                let got = window.data()[row * 40 + col];
                assert!(
                    (got - expected).abs() < 1e-12,
                    "mismatch at ({col},{row}): {got} vs {expected}"
                );
            }
        }
    }

    #[test]
    fn read_window_overview_matches_generated() {
        let (raster, bytes) = write_test_cog(64, 64, 1, 4326);
        let overviews = generate_overviews(&raster, 1);
        let mut reader = CogReader::open(bytes.as_slice()).unwrap();
        assert_eq!(reader.levels().len(), 2);
        let window = reader.read_window(1, 0, 0, 32, 32).unwrap();
        for i in 0..32 * 32 {
            assert!((window.data()[i] - overviews[0].data[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn read_window_pads_outside_with_nan() {
        let (_, bytes) = write_test_cog(20, 20, 0, 4326);
        let mut reader = CogReader::open(bytes.as_slice()).unwrap();
        let window = reader.read_window(0, 10, 10, 20, 20).unwrap();
        assert!(!window.data()[0].is_nan());
        assert!(window.data()[20 * 10 + 15].is_nan());
    }

    #[test]
    fn epsg_roundtrips_through_geokeys() {
        let (_, bytes) = write_test_cog(20, 20, 0, 32632);
        let reader = CogReader::open(bytes.as_slice()).unwrap();
        assert_eq!(reader.meta().epsg, 32632);
        assert!((reader.meta().origin_x - 10.0).abs() < 1e-12);
        assert!((reader.meta().pixel_width - 0.5).abs() < 1e-12);

        let (_, bytes) = write_test_cog(20, 20, 0, 4326);
        let reader = CogReader::open(bytes.as_slice()).unwrap();
        assert_eq!(reader.meta().epsg, 4326);
    }

    #[test]
    fn select_level_picks_coarsest_sufficient() {
        let (_, bytes) = write_test_cog(128, 128, 3, 4326);
        let reader = CogReader::open(bytes.as_slice()).unwrap();
        // levels at 0.5, 1.0, 2.0, 4.0
        assert_eq!(reader.select_level(0.5), 0);
        assert_eq!(reader.select_level(0.1), 0);
        assert_eq!(reader.select_level(1.4), 1);
        assert_eq!(reader.select_level(2.0), 2);
        assert_eq!(reader.select_level(100.0), 3);
    }

    #[test]
    fn window_read_fetches_a_fraction_of_the_file() {
        let (_, bytes) = write_test_cog(512, 512, 0, 4326);
        let total = bytes.len();
        let mut reader = CogReader::open(RecordingReader {
            data: &bytes,
            fetched: 0,
        })
        .unwrap();
        let window = reader.read_window(0, 100, 100, 16, 16).unwrap();
        assert!(!window.data()[0].is_nan());
        let fetched = reader.source.fetched;
        assert!(
            fetched * 10 < total,
            "fetched {fetched} of {total} bytes for one small window"
        );
    }

    #[test]
    fn overview_ifd_declares_scaled_pixel_size() {
        let (_, bytes) = write_test_cog(64, 64, 2, 4326);
        let reader = CogReader::open(bytes.as_slice()).unwrap();
        let widths: Vec<f64> = reader.levels().iter().map(|l| l.pixel_width).collect();
        assert_eq!(widths, vec![0.5, 1.0, 2.0]);
    }

    #[test]
    fn deflate_roundtrips_and_shrinks_the_file() {
        let raster = test_raster(100, 80);
        let mut params = CogParams {
            tile_width: 16,
            tile_height: 16,
            overview_levels: 2,
            ..CogParams::default()
        };
        let mut plain = io::Cursor::new(Vec::new());
        write_cog(&raster, &params, &mut plain).unwrap();
        params.deflate = true;
        let mut packed = io::Cursor::new(Vec::new());
        write_cog(&raster, &params, &mut packed).unwrap();
        let (plain, packed) = (plain.into_inner(), packed.into_inner());
        assert!(
            packed.len() < plain.len(),
            "deflate did not shrink the file"
        );

        let mut reader = CogReader::open(packed.as_slice()).unwrap();
        let window = reader.read_window(0, 0, 0, 100, 80).unwrap();
        for (a, b) in window.data().iter().zip(raster.data()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    /// minimal single-tile tiff with an arbitrary sample layout, for
    /// exercising decode paths the writer does not produce
    #[allow(clippy::too_many_arguments)]
    fn build_test_tiff(
        width: u32,
        height: u32,
        bits: u32,
        format: u32,
        compression: u16,
        predictor: u16,
        nodata: Option<&str>,
        tile: &[u8],
    ) -> Vec<u8> {
        let nod: Option<Vec<u8>> = nodata.map(|s| {
            let mut v = s.as_bytes().to_vec();
            v.push(0);
            v
        });
        let n: u32 = if nod.is_some() { 16 } else { 15 };
        let aux = 8 + 2 + n * 12 + 4;
        let (scale_off, tp_off, gk_off, nod_off) = (aux, aux + 24, aux + 72, aux + 88);
        let nod_stored = nod.as_ref().filter(|v| v.len() > 4);
        let tile_off = nod_off + nod_stored.map_or(0, |v| v.len() as u32);

        let mut buf = Vec::new();
        buf.extend_from_slice(b"II");
        push_u16(&mut buf, 42);
        push_u32(&mut buf, 8);
        push_u16(&mut buf, n as u16);
        push_entry(&mut buf, 256, 3, 1, width);
        push_entry(&mut buf, 257, 3, 1, height);
        push_entry(&mut buf, 258, 3, 1, bits);
        push_entry(&mut buf, 259, 3, 1, u32::from(compression));
        push_entry(&mut buf, 262, 3, 1, 1);
        push_entry(&mut buf, 277, 3, 1, 1);
        push_entry(&mut buf, 317, 3, 1, u32::from(predictor));
        push_entry(&mut buf, 322, 3, 1, width);
        push_entry(&mut buf, 323, 3, 1, height);
        push_entry(&mut buf, 324, 4, 1, tile_off);
        push_entry(&mut buf, 325, 4, 1, tile.len() as u32);
        push_entry(&mut buf, 339, 3, 1, format);
        push_entry(&mut buf, 33550, 12, 3, scale_off);
        push_entry(&mut buf, 33922, 12, 6, tp_off);
        push_entry(&mut buf, 34735, 3, 8, gk_off);
        if let Some(v) = &nod {
            let value = if v.len() <= 4 {
                let mut b = [0u8; 4];
                b[..v.len()].copy_from_slice(v);
                u32::from_le_bytes(b)
            } else {
                nod_off
            };
            push_entry(&mut buf, 42113, 2, v.len() as u32, value);
        }
        push_u32(&mut buf, 0);

        for v in [1.0f64, 1.0, 0.0] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for v in [0.0f64, 0.0, 0.0, 5.0, 55.0, 0.0] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for v in [1u16, 1, 0, 1, 2048, 0, 1, 4326] {
            push_u16(&mut buf, v);
        }
        if let Some(v) = nod_stored {
            buf.extend_from_slice(v);
        }
        buf.extend_from_slice(tile);
        buf
    }

    fn zlib(bytes: &[u8]) -> Vec<u8> {
        let mut enc = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(bytes).unwrap();
        enc.finish().unwrap()
    }

    fn read_all(bytes: &[u8], w: usize, h: usize) -> Vec<f64> {
        let mut reader = CogReader::open(bytes).unwrap();
        reader.read_window(0, 0, 0, w, h).unwrap().data().to_vec()
    }

    #[test]
    fn uint16_and_int16_tiles_decode() {
        let vals: Vec<u16> = (0..256u16).map(|i| i * 17 % 991).collect();
        let raw: Vec<u8> = vals.iter().flat_map(|v| v.to_le_bytes()).collect();
        let tiff = build_test_tiff(16, 16, 16, 1, 1, 1, None, &raw);
        let got = read_all(&tiff, 16, 16);
        for (g, v) in got.iter().zip(&vals) {
            assert_eq!(*g, f64::from(*v));
        }

        let vals: Vec<i16> = (0..256i32).map(|i| (i * 89 % 1000 - 500) as i16).collect();
        let raw: Vec<u8> = vals.iter().flat_map(|v| v.to_le_bytes()).collect();
        let tiff = build_test_tiff(16, 16, 16, 2, 1, 1, None, &raw);
        let got = read_all(&tiff, 16, 16);
        for (g, v) in got.iter().zip(&vals) {
            assert_eq!(*g, f64::from(*v));
        }
    }

    #[test]
    fn float32_tile_decodes() {
        let vals: Vec<f32> = (0..256).map(|i| i as f32 * 1.25 - 100.5).collect();
        let raw: Vec<u8> = vals.iter().flat_map(|v| v.to_le_bytes()).collect();
        let tiff = build_test_tiff(16, 16, 32, 3, 1, 1, None, &raw);
        let got = read_all(&tiff, 16, 16);
        for (g, v) in got.iter().zip(&vals) {
            assert_eq!(*g, f64::from(*v));
        }
    }

    #[test]
    fn horizontal_predictor_deflate_uint16_decodes() {
        let vals: Vec<u16> = (0..256u16).map(|i| 1000 + (i % 16) * 3 + i / 16).collect();
        let mut diffed = Vec::with_capacity(256);
        for row in vals.chunks(16) {
            diffed.push(row[0]);
            for i in 1..16 {
                diffed.push(row[i].wrapping_sub(row[i - 1]));
            }
        }
        let raw: Vec<u8> = diffed.iter().flat_map(|v| v.to_le_bytes()).collect();
        let tiff = build_test_tiff(16, 16, 16, 1, 8, 2, None, &zlib(&raw));
        let got = read_all(&tiff, 16, 16);
        for (g, v) in got.iter().zip(&vals) {
            assert_eq!(*g, f64::from(*v));
        }
    }

    #[test]
    fn floating_point_predictor_decodes() {
        let vals: Vec<f32> = (0..256).map(|i| (i as f32).sin() * 500.0).collect();
        let width = 16usize;
        let mut encoded = Vec::new();
        for row in vals.chunks(width) {
            // byte planes, most significant first, then forward differencing
            let mut planes = vec![0u8; width * 4];
            for (i, v) in row.iter().enumerate() {
                for (k, b) in v.to_be_bytes().iter().enumerate() {
                    planes[k * width + i] = *b;
                }
            }
            for i in (1..planes.len()).rev() {
                planes[i] = planes[i].wrapping_sub(planes[i - 1]);
            }
            encoded.extend_from_slice(&planes);
        }
        let tiff = build_test_tiff(16, 16, 32, 3, 8, 3, None, &zlib(&encoded));
        let got = read_all(&tiff, 16, 16);
        for (g, v) in got.iter().zip(&vals) {
            assert_eq!(*g, f64::from(*v));
        }
    }

    #[test]
    fn nodata_maps_to_nan_inline_and_stored() {
        // short ascii fits in the value field
        let vals: Vec<u16> = (0..256u16)
            .map(|i| if i % 7 == 0 { 0 } else { i })
            .collect();
        let raw: Vec<u8> = vals.iter().flat_map(|v| v.to_le_bytes()).collect();
        let tiff = build_test_tiff(16, 16, 16, 1, 1, 1, Some("0"), &raw);
        let got = read_all(&tiff, 16, 16);
        for (g, v) in got.iter().zip(&vals) {
            if *v == 0 {
                assert!(g.is_nan(), "nodata value survived: {g}");
            } else {
                assert_eq!(*g, f64::from(*v));
            }
        }

        // long ascii goes through the offset path
        let vals: Vec<f32> = (0..256)
            .map(|i| if i % 5 == 0 { -9999.25 } else { i as f32 })
            .collect();
        let raw: Vec<u8> = vals.iter().flat_map(|v| v.to_le_bytes()).collect();
        let tiff = build_test_tiff(16, 16, 32, 3, 1, 1, Some("-9999.25"), &raw);
        let got = read_all(&tiff, 16, 16);
        for (i, (g, v)) in got.iter().zip(&vals).enumerate() {
            if i % 5 == 0 {
                assert!(g.is_nan());
            } else {
                assert_eq!(*g, f64::from(*v));
            }
        }
    }

    fn test_bands(width: usize, height: usize, count: usize) -> BandedRaster {
        let bands = (0..count)
            .map(|b| {
                let data = (0..width * height)
                    .map(|i| i as f64 * 0.25 + b as f64 * 1000.0)
                    .collect();
                Raster::from_vec(width, height, data, 1.0, -9999.0).unwrap()
            })
            .collect();
        BandedRaster::new(bands).unwrap()
    }

    fn banded_params(deflate: bool) -> CogParams {
        CogParams {
            tile_width: 16,
            tile_height: 16,
            overview_levels: 1,
            epsg: 4326,
            origin_x: 10.0,
            origin_y: 50.0,
            pixel_width: 0.5,
            pixel_height: 0.5,
            deflate,
            ..CogParams::default()
        }
    }

    fn write_test_cog_bands(bands: &BandedRaster, deflate: bool) -> Vec<u8> {
        let mut buf = io::Cursor::new(Vec::new());
        write_cog_bands(bands, &banded_params(deflate), &mut buf).unwrap();
        buf.into_inner()
    }

    fn assert_close(got: &[f64], want: &[f64], what: &str) {
        assert_eq!(got.len(), want.len(), "{what}: length");
        for (i, (g, w)) in got.iter().zip(want).enumerate() {
            assert!((g - w).abs() < 1e-12, "{what}: {g} vs {w} at {i}");
        }
    }

    #[test]
    fn banded_roundtrips_at_base_and_overview() {
        for count in [2, 3, 4] {
            for deflate in [false, true] {
                let bands = test_bands(40, 30, count);
                let bytes = write_test_cog_bands(&bands, deflate);
                let mut reader = CogReader::open(bytes.as_slice()).unwrap();
                assert_eq!(reader.levels().len(), 2);
                assert_eq!(reader.levels()[0].samples, count);

                let base = reader.read_window_bands(0, 0, 0, 40, 30).unwrap();
                assert_eq!(base.band_count(), count);
                for b in 0..count {
                    assert_close(
                        base.band(b).unwrap().data(),
                        bands.band(b).unwrap().data(),
                        &format!("band {b} of {count}, deflate {deflate}"),
                    );
                }

                let ov = reader.read_window_bands(1, 0, 0, 20, 15).unwrap();
                for b in 0..count {
                    let expected = generate_overviews(bands.band(b).unwrap(), 1);
                    assert_close(
                        ov.band(b).unwrap().data(),
                        &expected[0].data,
                        &format!("overview band {b} of {count}"),
                    );
                }
            }
        }
    }

    #[test]
    fn banded_window_offset_from_the_origin() {
        let bands = test_bands(40, 30, 3);
        let bytes = write_test_cog_bands(&bands, true);
        let mut reader = CogReader::open(bytes.as_slice()).unwrap();
        let window = reader.read_window_bands(0, 7, 5, 20, 18).unwrap();
        for b in 0..3 {
            let src = bands.band(b).unwrap().data();
            let got = window.band(b).unwrap().data();
            for row in 0..18 {
                for col in 0..20 {
                    let want = src[(row + 5) * 40 + (col + 7)];
                    assert!((got[row * 20 + col] - want).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn read_window_on_a_multi_band_cog_points_at_read_window_bands() {
        let bytes = write_test_cog_bands(&test_bands(20, 20, 3), false);
        let mut reader = CogReader::open(bytes.as_slice()).unwrap();
        let err = reader.read_window(0, 0, 0, 4, 4).unwrap_err().to_string();
        assert!(err.contains("read_window_bands"), "unhelpful error: {err}");
        assert!(err.contains('3'), "error omits the sample count: {err}");
    }

    #[test]
    fn single_band_agrees_through_both_paths() {
        let (raster, bytes) = write_test_cog(64, 48, 1, 4326);
        let mut reader = CogReader::open(bytes.as_slice()).unwrap();
        let plain = reader.read_window(0, 3, 4, 20, 20).unwrap();
        let banded = reader.read_window_bands(0, 3, 4, 20, 20).unwrap();
        assert_eq!(banded.band_count(), 1);
        assert_eq!(plain.data(), banded.band(0).unwrap().data());
        for row in 0..20 {
            for col in 0..20 {
                let want = raster.data()[(row + 4) * 64 + (col + 3)];
                assert!((plain.data()[row * 20 + col] - want).abs() < 1e-12);
            }
        }

        // one band through write_cog_bands matches write_cog
        let one = BandedRaster::new(vec![raster.clone()]).unwrap();
        let mut via_bands = io::Cursor::new(Vec::new());
        write_cog_bands(&one, &banded_params(false), &mut via_bands).unwrap();
        let mut via_plain = io::Cursor::new(Vec::new());
        write_cog(&raster, &banded_params(false), &mut via_plain).unwrap();
        assert_eq!(via_bands.into_inner(), via_plain.into_inner());
    }

    /// single-tile multi-band tiff with explicit per-sample bits and
    /// formats, for layouts the writer does not produce
    fn build_banded_tiff(
        width: u32,
        height: u32,
        bits: &[u16],
        formats: &[u16],
        planar: u16,
        tile: &[u8],
    ) -> Vec<u8> {
        let samples = bits.len() as u32;
        let n: u32 = 15;
        let aux = 8 + 2 + n * 12 + 4;
        let (scale_off, tp_off, gk_off) = (aux, aux + 24, aux + 72);
        let arrays_off = aux + 88;
        let tile_off = arrays_off + if samples > 2 { samples * 4 } else { 0 };
        let pack = |v: &[u16]| -> u32 {
            let hi = v.get(1).copied().unwrap_or(0);
            u32::from(v[0]) | (u32::from(hi) << 16)
        };

        let mut buf = Vec::new();
        buf.extend_from_slice(b"II");
        push_u16(&mut buf, 42);
        push_u32(&mut buf, 8);
        push_u16(&mut buf, n as u16);
        push_entry(&mut buf, 256, 3, 1, width);
        push_entry(&mut buf, 257, 3, 1, height);
        let bits_value = if samples > 2 { arrays_off } else { pack(bits) };
        push_entry(&mut buf, 258, 3, samples, bits_value);
        push_entry(&mut buf, 259, 3, 1, 1);
        push_entry(&mut buf, 262, 3, 1, 1);
        push_entry(&mut buf, 277, 3, 1, samples);
        push_entry(&mut buf, 284, 3, 1, u32::from(planar));
        push_entry(&mut buf, 322, 3, 1, width);
        push_entry(&mut buf, 323, 3, 1, height);
        push_entry(&mut buf, 324, 4, 1, tile_off);
        push_entry(&mut buf, 325, 4, 1, tile.len() as u32);
        let formats_value = if samples > 2 {
            arrays_off + samples * 2
        } else {
            pack(formats)
        };
        push_entry(&mut buf, 339, 3, samples, formats_value);
        push_entry(&mut buf, 33550, 12, 3, scale_off);
        push_entry(&mut buf, 33922, 12, 6, tp_off);
        push_entry(&mut buf, 34735, 3, 8, gk_off);
        push_u32(&mut buf, 0);

        for v in [1.0f64, 1.0, 0.0] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for v in [0.0f64, 0.0, 0.0, 5.0, 55.0, 0.0] {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for v in [1u16, 1, 0, 1, 2048, 0, 1, 4326] {
            push_u16(&mut buf, v);
        }
        if samples > 2 {
            for v in bits.iter().chain(formats) {
                push_u16(&mut buf, *v);
            }
        }
        buf.extend_from_slice(tile);
        buf
    }

    #[test]
    fn foreign_multi_band_tile_deinterleaves() {
        // three uint16 bands, pixel-interleaved, not written by write_cog
        let pixels = 16 * 16;
        let bands: Vec<Vec<u16>> = (0..3)
            .map(|b| (0..pixels).map(|i| (i * 7 + b * 3000) as u16).collect())
            .collect();
        let mut raw = Vec::new();
        for p in 0..pixels {
            for band in &bands {
                raw.extend_from_slice(&band[p].to_le_bytes());
            }
        }
        let tiff = build_banded_tiff(16, 16, &[16; 3], &[1; 3], 1, &raw);
        let mut reader = CogReader::open(tiff.as_slice()).unwrap();
        let got = reader.read_window_bands(0, 0, 0, 16, 16).unwrap();
        assert_eq!(got.band_count(), 3);
        for (b, band) in bands.iter().enumerate() {
            let want: Vec<f64> = band.iter().map(|&v| f64::from(v)).collect();
            assert_close(got.band(b).unwrap().data(), &want, &format!("band {b}"));
        }
    }

    fn open_err(bytes: &[u8]) -> String {
        match CogReader::open(bytes) {
            Err(e) => e.to_string(),
            Ok(_) => panic!("expected open to reject this file"),
        }
    }

    #[test]
    fn mismatched_sample_layouts_error_at_open() {
        let raw = vec![0u8; 16 * 16 * 2 * 8];
        // two samples, one 64-bit and one 32-bit, packed in the value field
        let tiff = build_banded_tiff(16, 16, &[64, 32], &[3, 3], 1, &raw);
        let err = open_err(&tiff);
        assert!(err.contains("BitsPerSample"), "{err}");

        // three samples, formats read from the offset array
        let raw = vec![0u8; 16 * 16 * 3 * 8];
        let tiff = build_banded_tiff(16, 16, &[64; 3], &[3, 3, 1], 1, &raw);
        let err = open_err(&tiff);
        assert!(err.contains("SampleFormat"), "{err}");
    }

    #[test]
    fn band_interleaved_planar_config_is_rejected() {
        let raw = vec![0u8; 16 * 16 * 2 * 8];
        let tiff = build_banded_tiff(16, 16, &[64; 2], &[3; 2], 2, &raw);
        let err = open_err(&tiff);
        assert!(err.contains("PlanarConfiguration"), "{err}");
    }

    // --- cloud optimized layout ---
    //
    // these assert the structure that separates a cog from a merely tiled
    // tiff, the parts a reader relies on to fetch one zoom level without
    // pulling the whole file

    /// one ifd: its own offset and its entries as (tag, type, count, value)
    struct Ifd {
        offset: u32,
        entries: Vec<(u16, u16, u32, u32)>,
    }

    impl Ifd {
        fn entry(&self, tag: u16) -> Option<(u16, u32, u32)> {
            self.entries
                .iter()
                .find(|e| e.0 == tag)
                .map(|&(_, ty, count, value)| (ty, count, value))
        }

        fn value(&self, tag: u16) -> Option<u32> {
            self.entry(tag).map(|(_, _, value)| value)
        }

        /// tile offsets, inline for a single tile and out of line otherwise
        fn tile_offsets(&self, data: &[u8]) -> Vec<u32> {
            let (_, count, value) = self.entry(324).expect("TileOffsets");
            if count == 1 {
                return vec![value];
            }
            (0..count as usize)
                .map(|i| u32_at(data, value as usize + i * 4))
                .collect()
        }
    }

    fn u16_at(data: &[u8], at: usize) -> u16 {
        u16::from_le_bytes(data[at..at + 2].try_into().unwrap())
    }

    fn u32_at(data: &[u8], at: usize) -> u32 {
        u32::from_le_bytes(data[at..at + 4].try_into().unwrap())
    }

    fn f64_at(data: &[u8], at: usize) -> f64 {
        f64::from_le_bytes(data[at..at + 8].try_into().unwrap())
    }

    /// walk the next-ifd chain from the header
    fn ifd_chain(data: &[u8]) -> Vec<Ifd> {
        let mut chain = Vec::new();
        let mut next = u32_at(data, 4);
        while next != 0 {
            let at = next as usize;
            let count = u16_at(data, at) as usize;
            let entries = (0..count)
                .map(|i| {
                    let e = at + 2 + i * 12;
                    (
                        u16_at(data, e),
                        u16_at(data, e + 2),
                        u32_at(data, e + 4),
                        u32_at(data, e + 8),
                    )
                })
                .collect();
            chain.push(Ifd {
                offset: next,
                entries,
            });
            next = u32_at(data, at + 2 + count * 12);
        }
        chain
    }

    /// the ascii text of a GDAL_NODATA entry, inline or out of line
    fn nodata_text_of(data: &[u8], ifd: &Ifd) -> Option<String> {
        let (ty, count, value) = ifd.entry(42113)?;
        assert_eq!(ty, 2, "GDAL_NODATA is ascii");
        let bytes = if count as usize <= INLINE_BYTES {
            value.to_le_bytes()[..count as usize].to_vec()
        } else {
            data[value as usize..value as usize + count as usize].to_vec()
        };
        assert_eq!(bytes.last(), Some(&0), "ascii values are nul terminated");
        Some(String::from_utf8(bytes[..bytes.len() - 1].to_vec()).unwrap())
    }

    #[test]
    fn only_overview_ifds_are_marked_reduced_resolution() {
        let (_, bytes) = write_test_cog(100, 80, 3, 4326);
        let chain = ifd_chain(&bytes);
        assert_eq!(chain.len(), 4, "full resolution plus three overviews");

        assert_eq!(
            chain[0].entry(254),
            None,
            "the full resolution ifd is not a reduced-resolution image"
        );
        for (i, ifd) in chain[1..].iter().enumerate() {
            let (ty, count, value) = ifd
                .entry(254)
                .unwrap_or_else(|| panic!("overview {i} needs NewSubfileType"));
            assert_eq!((ty, count, value), (4, 1, 1), "overview {i}");
        }
    }

    #[test]
    fn ifds_precede_all_tile_data_and_ascend() {
        let (_, bytes) = write_test_cog(100, 80, 3, 4326);
        let chain = ifd_chain(&bytes);

        assert_eq!(chain[0].offset, 8, "the main ifd follows the header");
        for pair in chain.windows(2) {
            assert!(
                pair[1].offset > pair[0].offset,
                "ifds ascend: {} then {}",
                pair[0].offset,
                pair[1].offset
            );
        }

        let last_ifd = chain.last().unwrap().offset;
        let first_tile = chain
            .iter()
            .flat_map(|ifd| ifd.tile_offsets(&bytes))
            .min()
            .unwrap();
        assert!(
            first_tile > last_ifd,
            "tile data at {first_tile} must follow the last ifd at {last_ifd}"
        );
    }

    #[test]
    fn tile_data_runs_smallest_overview_first() {
        let (_, bytes) = write_test_cog(100, 80, 3, 4326);
        let chain = ifd_chain(&bytes);
        let starts: Vec<u32> = chain
            .iter()
            .map(|ifd| *ifd.tile_offsets(&bytes).iter().min().unwrap())
            .collect();

        // level 0 is full resolution and the last level is the smallest
        // overview, so the offsets descend as resolution rises
        for level in 0..starts.len() - 1 {
            assert!(
                starts[level] > starts[level + 1],
                "level {level} data at {} should follow level {} at {}",
                starts[level],
                level + 1,
                starts[level + 1]
            );
        }
    }

    #[test]
    fn tiles_are_row_major_and_cover_each_level() {
        let (_, bytes) = write_test_cog(100, 80, 2, 4326);
        let chain = ifd_chain(&bytes);
        let sizes = [(100u32, 80u32), (50, 40), (25, 20)];

        for (ifd, (width, height)) in chain.iter().zip(sizes) {
            assert_eq!(ifd.value(256), Some(width));
            assert_eq!(ifd.value(257), Some(height));
            assert_eq!(ifd.value(322), Some(16), "TileWidth");
            assert_eq!(ifd.value(323), Some(16), "TileLength");

            let expected = width.div_ceil(16) * height.div_ceil(16);
            let offsets = ifd.tile_offsets(&bytes);
            assert_eq!(
                offsets.len() as u32,
                expected,
                "tile count for {width}x{height}"
            );
            for pair in offsets.windows(2) {
                assert!(pair[1] > pair[0], "tiles ascend in row-major order");
            }
        }
    }

    #[test]
    fn geo_tags_land_on_every_level() {
        let (_, bytes) = write_test_cog(100, 80, 2, 32632);
        let chain = ifd_chain(&bytes);

        for (level, ifd) in chain.iter().enumerate() {
            let scale = ifd.value(33550).expect("ModelPixelScale") as usize;
            let tiepoint = ifd.value(33922).expect("ModelTiepoint") as usize;
            let keys = ifd.value(34735).expect("GeoKeyDirectory") as usize;

            // each overview halves resolution, so its pixels are wider
            let factor = 2f64.powi(level as i32);
            assert_eq!(f64_at(&bytes, scale), 0.5 * factor, "level {level} x scale");
            assert_eq!(
                f64_at(&bytes, scale + 8),
                0.5 * factor,
                "level {level} y scale"
            );

            // the tiepoint ties raster (0,0) to the map origin
            assert_eq!(f64_at(&bytes, tiepoint + 24), 10.0);
            assert_eq!(f64_at(&bytes, tiepoint + 32), 50.0);

            assert_eq!(u16_at(&bytes, keys + 6), 3, "three geo keys");
            // projected model type, pixel-is-area, then the crs
            assert_eq!(u16_at(&bytes, keys + 8), 1024);
            assert_eq!(u16_at(&bytes, keys + 14), 1);
            assert_eq!(u16_at(&bytes, keys + 16), 1025);
            assert_eq!(u16_at(&bytes, keys + 22), 1);
            assert_eq!(u16_at(&bytes, keys + 24), 3072);
            assert_eq!(u16_at(&bytes, keys + 30), 32632);
        }
    }

    #[test]
    fn geographic_crs_uses_the_geographic_geo_keys() {
        let (_, bytes) = write_test_cog(20, 20, 0, 4326);
        let keys = ifd_chain(&bytes)[0].value(34735).unwrap() as usize;
        assert_eq!(u16_at(&bytes, keys + 14), 2, "geographic model type");
        assert_eq!(u16_at(&bytes, keys + 24), 2048, "GeographicTypeGeoKey");
        assert_eq!(u16_at(&bytes, keys + 30), 4326);
    }

    fn write_cog_with(raster: &Raster, nodata: Option<f64>) -> Vec<u8> {
        let params = CogParams {
            tile_width: 16,
            tile_height: 16,
            overview_levels: 1,
            epsg: 4326,
            nodata,
            ..CogParams::default()
        };
        let mut buf = io::Cursor::new(Vec::new());
        write_cog(raster, &params, &mut buf).unwrap();
        buf.into_inner()
    }

    #[test]
    fn nodata_is_declared_on_every_level() {
        let raster = test_raster(40, 40);
        let bytes = write_cog_with(&raster, Some(-9999.0));
        for (level, ifd) in ifd_chain(&bytes).iter().enumerate() {
            assert_eq!(
                nodata_text_of(&bytes, ifd).as_deref(),
                Some("-9999"),
                "level {level}"
            );
        }
    }

    #[test]
    fn nan_nodata_is_spelled_the_way_gdal_spells_it() {
        let raster = test_raster(20, 20);
        let bytes = write_cog_with(&raster, Some(f64::NAN));
        let chain = ifd_chain(&bytes);
        assert_eq!(nodata_text_of(&bytes, &chain[0]).as_deref(), Some("nan"));
    }

    #[test]
    fn no_nodata_tag_when_none_is_asked_for() {
        let raster = test_raster(20, 20);
        let bytes = write_cog_with(&raster, None);
        assert_eq!(ifd_chain(&bytes)[0].entry(42113), None);
    }

    #[test]
    fn nan_samples_become_the_declared_nodata_and_read_back_as_nan() {
        let mut raster = test_raster(40, 40);
        raster.set(3, 5, f64::NAN);
        raster.set(20, 21, f64::NAN);
        let bytes = write_cog_with(&raster, Some(-9999.0));

        // the file stores the declared value, not NaN
        let ifd = &ifd_chain(&bytes)[0];
        let tile = ifd.tile_offsets(&bytes)[0] as usize;
        assert_eq!(f64_at(&bytes, tile + (3 * 16 + 5) * 8), -9999.0);

        // and terrano's reader maps it back to NaN
        let mut reader = CogReader::open(bytes.as_slice()).unwrap();
        let window = reader.read_window(0, 0, 0, 40, 40).unwrap();
        assert!(window.get(3, 5).unwrap().is_nan());
        assert!(window.get(20, 21).unwrap().is_nan());
        assert_eq!(window.get(0, 1).unwrap(), 1.0);
    }

    #[test]
    fn multi_band_declares_its_extra_samples() {
        let bands = BandedRaster::new(vec![
            test_raster(40, 40),
            test_raster(40, 40),
            test_raster(40, 40),
        ])
        .unwrap();
        let bytes = write_test_cog_bands(&bands, false);
        let ifd = &ifd_chain(&bytes)[0];

        assert_eq!(ifd.value(277), Some(3), "SamplesPerPixel");
        // libtiff wants colour channels plus extra samples to reach
        // SamplesPerPixel, and min-is-black contributes one
        let (ty, count, value) = ifd.entry(338).expect("ExtraSamples");
        assert_eq!((ty, count), (3, 2));
        assert_eq!(value, 0, "both extra bands are unspecified data");
    }

    #[test]
    fn deflate_tiles_are_smaller_and_still_read_back() {
        let raster = test_raster(64, 64);
        let plain = write_cog_with(&raster, Some(-9999.0));
        let params = CogParams {
            tile_width: 16,
            tile_height: 16,
            overview_levels: 1,
            epsg: 4326,
            deflate: true,
            nodata: Some(-9999.0),
            ..CogParams::default()
        };
        let mut buf = io::Cursor::new(Vec::new());
        write_cog(&raster, &params, &mut buf).unwrap();
        let squeezed = buf.into_inner();

        assert_eq!(ifd_chain(&plain)[0].value(259), Some(1), "no compression");
        assert_eq!(ifd_chain(&squeezed)[0].value(259), Some(8), "deflate");
        assert!(
            squeezed.len() < plain.len(),
            "deflate should shrink a smooth ramp: {} vs {}",
            squeezed.len(),
            plain.len()
        );

        let mut reader = CogReader::open(squeezed.as_slice()).unwrap();
        let window = reader.read_window(0, 0, 0, 64, 64).unwrap();
        for i in 0..64 * 64 {
            assert_eq!(window.data()[i], raster.data()[i]);
        }
    }
}
