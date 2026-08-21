use crate::{BandedRaster, Error, Raster};
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

/// Sample layout used when writing raster values.
///
/// Integer formats round to the nearest whole number and clamp to their own
/// range, so an out-of-range value is pinned to the nearest end rather than
/// wrapping, and NaN lands on zero.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SampleFormat {
    U8,
    I8,
    U16,
    I16,
    U32,
    I32,
    F32,
    /// 64-bit IEEE floats, as written by [`write_geotiff`].
    F64,
}

impl SampleFormat {
    /// Every format, in the order used to resolve a name or a TIFF tag pair.
    pub const ALL: [SampleFormat; 8] = [
        SampleFormat::U8,
        SampleFormat::I8,
        SampleFormat::U16,
        SampleFormat::I16,
        SampleFormat::U32,
        SampleFormat::I32,
        SampleFormat::F32,
        SampleFormat::F64,
    ];

    /// TIFF BitsPerSample.
    pub fn bits(self) -> u16 {
        match self {
            SampleFormat::U8 | SampleFormat::I8 => 8,
            SampleFormat::U16 | SampleFormat::I16 => 16,
            SampleFormat::U32 | SampleFormat::I32 | SampleFormat::F32 => 32,
            SampleFormat::F64 => 64,
        }
    }

    pub fn bytes(self) -> u32 {
        u32::from(self.bits()) / 8
    }

    /// TIFF SampleFormat tag value: 1 = unsigned, 2 = signed, 3 = IEEE float.
    pub fn tag_value(self) -> u16 {
        match self {
            SampleFormat::U8 | SampleFormat::U16 | SampleFormat::U32 => 1,
            SampleFormat::I8 | SampleFormat::I16 | SampleFormat::I32 => 2,
            SampleFormat::F32 | SampleFormat::F64 => 3,
        }
    }

    /// Short name accepted by the wasm binding, e.g. `"u8"`.
    pub fn name(self) -> &'static str {
        match self {
            SampleFormat::U8 => "u8",
            SampleFormat::I8 => "i8",
            SampleFormat::U16 => "u16",
            SampleFormat::I16 => "i16",
            SampleFormat::U32 => "u32",
            SampleFormat::I32 => "i32",
            SampleFormat::F32 => "f32",
            SampleFormat::F64 => "f64",
        }
    }

    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|f| f.name() == name)
    }

    /// The format a TIFF declares through BitsPerSample and SampleFormat.
    pub fn from_tiff(bits: u16, tag_value: u16) -> Option<Self> {
        Self::ALL
            .into_iter()
            .find(|f| f.bits() == bits && f.tag_value() == tag_value)
    }

    /// Lowest and highest representable value, `None` for the float formats.
    pub fn integer_range(self) -> Option<(f64, f64)> {
        match self {
            SampleFormat::U8 => Some((0.0, 255.0)),
            SampleFormat::I8 => Some((-128.0, 127.0)),
            SampleFormat::U16 => Some((0.0, 65535.0)),
            SampleFormat::I16 => Some((-32768.0, 32767.0)),
            SampleFormat::U32 => Some((0.0, 4294967295.0)),
            SampleFormat::I32 => Some((-2147483648.0, 2147483647.0)),
            SampleFormat::F32 | SampleFormat::F64 => None,
        }
    }

    /// Append one little-endian sample.
    pub(crate) fn encode(self, value: f64, out: &mut Vec<u8>) {
        let value = match self.integer_range() {
            // a saturating cast, so NaN lands on 0
            Some((low, high)) => value.round().clamp(low, high),
            None => value,
        };
        match self {
            SampleFormat::U8 => out.push(value as u8),
            SampleFormat::I8 => out.extend_from_slice(&(value as i8).to_le_bytes()),
            SampleFormat::U16 => out.extend_from_slice(&(value as u16).to_le_bytes()),
            SampleFormat::I16 => out.extend_from_slice(&(value as i16).to_le_bytes()),
            SampleFormat::U32 => out.extend_from_slice(&(value as u32).to_le_bytes()),
            SampleFormat::I32 => out.extend_from_slice(&(value as i32).to_le_bytes()),
            SampleFormat::F32 => out.extend_from_slice(&(value as f32).to_le_bytes()),
            SampleFormat::F64 => out.extend_from_slice(&value.to_le_bytes()),
        }
    }

    /// Widen one little-endian sample from the start of `bytes`.
    pub(crate) fn decode(self, bytes: &[u8]) -> f64 {
        let two = || u16::from_le_bytes(bytes[..2].try_into().unwrap());
        let four = || u32::from_le_bytes(bytes[..4].try_into().unwrap());
        match self {
            SampleFormat::U8 => f64::from(bytes[0]),
            SampleFormat::I8 => f64::from(bytes[0] as i8),
            SampleFormat::U16 => f64::from(two()),
            SampleFormat::I16 => f64::from(two() as i16),
            SampleFormat::U32 => f64::from(four()),
            SampleFormat::I32 => f64::from(four() as i32),
            SampleFormat::F32 => f64::from(f32::from_bits(four())),
            SampleFormat::F64 => f64::from_le_bytes(bytes[..8].try_into().unwrap()),
        }
    }
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

    let tiepoint: [f64; 6] = [0.0, 0.0, 0.0, meta.origin_x, meta.origin_y, 0.0];
    let pixel_scale: [f64; 3] = [meta.pixel_width, meta.pixel_height, 0.0];
    let geo_keys = geo_keys_for(meta.epsg);

    let num_tags: u16 = 13;
    let ifd_start: u32 = 8;
    let ifd_size: u32 = 2 + (num_tags as u32) * 12 + 4; // count + entries + next_ifd

    // layout after the ifd, in the order these blocks are written below
    let tiepoint_offset = ifd_start + ifd_size;
    let pixel_scale_offset = tiepoint_offset + 48;
    let geo_keys_offset = pixel_scale_offset + 24;
    let strip_offset = geo_keys_offset + geo_keys.len() as u32 * 2;

    let mut buf: Vec<u8> = Vec::with_capacity((strip_offset + strip_byte_count) as usize);

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

/// Write a multi-band raster as a GeoTIFF, e.g. an RGB or RGBA image.
///
/// Samples are pixel-interleaved (PlanarConfiguration 1): every reader supports
/// chunky layout, and RGB(A) imagery arrives that way anyway.
///
/// PhotometricInterpretation is RGB for three or more bands, min-is-black
/// otherwise; a 4th band of an RGB image is tagged as unassociated alpha.
/// Band names are not stored in the file.
///
/// # Arguments
/// * `raster` — the bands to write
/// * `meta` — georeferencing metadata
/// * `format` — sample layout, e.g. [`SampleFormat::U8`] for image bands
/// * `writer` — output destination
pub fn write_geotiff_bands<W: Write>(
    raster: &BandedRaster,
    meta: &GeoTiffMetadata,
    format: SampleFormat,
    writer: &mut W,
) -> Result<(), Error> {
    let width = raster.width() as u32;
    let height = raster.height() as u32;
    let samples = raster.band_count() as u32;

    let strip_bytes =
        u64::from(width) * u64::from(height) * u64::from(samples) * u64::from(format.bytes());
    // classic TIFF offsets are 32-bit, so the whole file has to fit in 4 GiB
    if strip_bytes + u64::from(u16::MAX) > u64::from(u32::MAX) {
        return Err(Error::InvalidInput(
            "image too large for a classic TIFF".into(),
        ));
    }
    let strip_byte_count = strip_bytes as u32;

    let (photometric, named_samples): (u16, u32) = if samples >= 3 { (2, 3) } else { (1, 1) };
    // 2 = unassociated alpha for the band after RGB, 0 = unspecified for anything beyond
    let extra_samples: Vec<u16> = (0..samples - named_samples)
        .map(|i| if photometric == 2 && i == 0 { 2 } else { 0 })
        .collect();

    let bits_per_sample = vec![format.bits(); samples as usize];
    let sample_formats = vec![format.tag_value(); samples as usize];
    let geo_keys = geo_keys_for(meta.epsg);
    let tiepoint: [f64; 6] = [0.0, 0.0, 0.0, meta.origin_x, meta.origin_y, 0.0];
    let pixel_scale: [f64; 3] = [meta.pixel_width, meta.pixel_height, 0.0];

    let num_tags: u16 = if extra_samples.is_empty() { 14 } else { 15 };
    let ifd_start: u32 = 8;
    let mut next = ifd_start + 2 + u32::from(num_tags) * 12 + 4;

    // order of these allocations must match the order the blocks are written below
    let bits_value = alloc_shorts(&mut next, &bits_per_sample);
    let sample_format_value = alloc_shorts(&mut next, &sample_formats);
    let extra_samples_value = alloc_shorts(&mut next, &extra_samples);
    let pixel_scale_offset = next;
    next += 24;
    let tiepoint_offset = next;
    next += 48;
    let geo_keys_offset = next;
    next += geo_keys.len() as u32 * 2;
    let strip_offset = next;

    let mut buf: Vec<u8> = Vec::with_capacity((strip_offset + strip_byte_count) as usize);

    // TIFF header
    buf.write_all(b"II")?;
    write_u16(&mut buf, 42)?;
    write_u32(&mut buf, ifd_start)?;

    // IFD, sorted by tag number
    write_u16(&mut buf, num_tags)?;
    write_ifd_entry(&mut buf, 256, 4, 1, width)?; // ImageWidth
    write_ifd_entry(&mut buf, 257, 4, 1, height)?; // ImageLength
    write_ifd_entry(&mut buf, 258, 3, samples, bits_value)?; // BitsPerSample
    write_ifd_entry(&mut buf, 259, 3, 1, 1)?; // Compression=None
    write_ifd_entry(&mut buf, 262, 3, 1, u32::from(photometric))?; // PhotometricInterpretation
    write_ifd_entry(&mut buf, 273, 4, 1, strip_offset)?; // StripOffsets
    write_ifd_entry(&mut buf, 277, 3, 1, samples)?; // SamplesPerPixel
    write_ifd_entry(&mut buf, 278, 4, 1, height)?; // RowsPerStrip
    write_ifd_entry(&mut buf, 279, 4, 1, strip_byte_count)?; // StripByteCounts
    write_ifd_entry(&mut buf, 284, 3, 1, 1)?; // PlanarConfiguration=chunky
    if !extra_samples.is_empty() {
        write_ifd_entry(
            &mut buf,
            338,
            3,
            extra_samples.len() as u32,
            extra_samples_value,
        )?; // ExtraSamples
    }
    write_ifd_entry(&mut buf, 339, 3, samples, sample_format_value)?; // SampleFormat
    write_ifd_entry(&mut buf, 33550, 12, 3, pixel_scale_offset)?; // ModelPixelScaleTag
    write_ifd_entry(&mut buf, 33922, 12, 6, tiepoint_offset)?; // ModelTiepointTag
    write_ifd_entry(&mut buf, 34735, 3, geo_keys.len() as u32, geo_keys_offset)?; // GeoKeyDirectoryTag
    write_u32(&mut buf, 0)?; // no further IFDs

    write_shorts_block(&mut buf, &bits_per_sample)?;
    write_shorts_block(&mut buf, &sample_formats)?;
    write_shorts_block(&mut buf, &extra_samples)?;
    for &val in &pixel_scale {
        write_f64(&mut buf, val)?;
    }
    for &val in &tiepoint {
        write_f64(&mut buf, val)?;
    }
    for &val in &geo_keys {
        write_u16(&mut buf, val)?;
    }

    // strip data, pixel-interleaved
    for row in 0..raster.height() {
        for col in 0..raster.width() {
            for band in raster.bands() {
                let val = band.get(row, col).unwrap_or(band.nodata);
                format.encode(val, &mut buf);
            }
        }
    }

    writer.write_all(&buf)?;
    Ok(())
}

/// Read a multi-band GeoTIFF written by [`write_geotiff_bands`].
///
/// Supports uncompressed, pixel-interleaved files in any [`SampleFormat`].
/// Strips are concatenated using RowsPerStrip and StripOffsets. A single-band
/// file reads back as a one-band raster; band names are not carried by the
/// file, so bands come back unnamed.
pub fn read_geotiff_bands(data: &[u8]) -> Result<(BandedRaster, GeoTiffMetadata), Error> {
    let ifd_offset = first_ifd_offset(data)?;
    let num_entries = read_u16_at(data, ifd_offset) as usize;
    if ifd_offset + 2 + num_entries * 12 > data.len() {
        return Err(Error::Format("truncated IFD".into()));
    }

    let mut width = 0u32;
    let mut height = 0u32;
    let mut samples = 1u32;
    let mut planar = 1u16;
    let mut compression = 0u16;
    let mut rows_per_strip = 0u32;
    let mut strip_offsets: Vec<u32> = Vec::new();
    let mut bits_per_sample: Vec<u16> = Vec::new();
    let mut sample_formats: Vec<u16> = Vec::new();
    let mut tiepoint_offset = 0u32;
    let mut pixel_scale_offset = 0u32;
    let mut geo_keys_offset = 0u32;
    let mut geo_keys_count = 0u32;

    for i in 0..num_entries {
        let entry_off = ifd_offset + 2 + i * 12;
        let tag = read_u16_at(data, entry_off);
        let count = read_u32_at(data, entry_off + 4);
        let value = read_u32_at(data, entry_off + 8);

        match tag {
            256 => width = value,
            257 => height = value,
            258 => bits_per_sample = read_shorts(data, entry_off)?,
            259 => compression = (value & 0xffff) as u16,
            273 => strip_offsets = read_offset_array(data, entry_off)?,
            277 => samples = value,
            278 => rows_per_strip = value,
            284 => planar = value as u16,
            339 => sample_formats = read_shorts(data, entry_off)?,
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
        return Err(Error::Format("missing width/height".into()));
    }
    if samples == 0 {
        return Err(Error::Format("SamplesPerPixel is 0".into()));
    }
    if compression != 1 {
        return Err(Error::Format(
            "only uncompressed (Compression=1) GeoTIFF is supported".into(),
        ));
    }
    if planar != 1 {
        return Err(Error::Format(
            "only pixel-interleaved (PlanarConfiguration 1) files supported".into(),
        ));
    }
    if pixel_scale_offset == 0 || tiepoint_offset == 0 {
        return Err(Error::Format("missing georeferencing tags".into()));
    }
    if strip_offsets.is_empty() {
        return Err(Error::Format("missing StripOffsets".into()));
    }

    // BitsPerSample and SampleFormat may be absent (defaults 1 and unsigned int)
    let bits = *bits_per_sample.first().unwrap_or(&1);
    let sample_format = *sample_formats.first().unwrap_or(&1);
    if bits_per_sample.iter().any(|&b| b != bits)
        || sample_formats.iter().any(|&f| f != sample_format)
    {
        return Err(Error::Format("mixed sample types not supported".into()));
    }
    let Some(format) = SampleFormat::from_tiff(bits, sample_format) else {
        return Err(Error::Format(format!(
            "unsupported sample layout: {bits} bits, format {sample_format}"
        )));
    };

    let rows_per_strip = if rows_per_strip == 0 || rows_per_strip == u32::MAX {
        height
    } else {
        rows_per_strip
    };
    let n_strips = (height as usize).div_ceil(rows_per_strip as usize);
    if strip_offsets.len() != n_strips {
        return Err(Error::Format(format!(
            "expected {n_strips} strips, found {}",
            strip_offsets.len()
        )));
    }

    let sample_bytes = format.bytes() as usize;
    let pixels = width as usize * height as usize;
    let mut band_values: Vec<Vec<f64>> = (0..samples as usize).map(|_| vec![0.0; pixels]).collect();

    for (strip_i, &strip_off) in strip_offsets.iter().enumerate() {
        let row0 = strip_i * rows_per_strip as usize;
        let rows = (rows_per_strip as usize).min(height as usize - row0);
        let strip_pixels = rows * width as usize;
        let strip_bytes = strip_pixels * samples as usize * sample_bytes;
        let base = strip_off as usize;
        if base + strip_bytes > data.len() {
            return Err(Error::Format("truncated strip data".into()));
        }
        for pixel in 0..strip_pixels {
            let dest = row0 * width as usize + pixel;
            for (band, values) in band_values.iter_mut().enumerate() {
                let off = base + (pixel * samples as usize + band) * sample_bytes;
                values[dest] = format.decode(&data[off..]);
            }
        }
    }

    let meta = geo_meta_from(
        data,
        pixel_scale_offset,
        tiepoint_offset,
        geo_keys_offset,
        geo_keys_count,
    );

    let nodata = -9999.0;
    let mut bands = Vec::with_capacity(samples as usize);
    for values in band_values {
        bands.push(Raster::from_vec(
            width as usize,
            height as usize,
            values,
            meta.pixel_width,
            nodata,
        )?);
    }

    Ok((BandedRaster::new(bands)?, meta))
}

/// Read an uncompressed GeoTIFF as a single-band raster (band 0).
///
/// Honours BitsPerSample, SampleFormat, Compression and RowsPerStrip for the
/// same uncompressed, pixel-interleaved files as [`read_geotiff_bands`].
pub fn read_geotiff(data: &[u8]) -> Result<(Raster, GeoTiffMetadata), Error> {
    let (banded, meta) = read_geotiff_bands(data)?;
    let raster = banded
        .into_bands()
        .into_iter()
        .next()
        .ok_or_else(|| Error::Format("GeoTIFF has no bands".into()))?;
    Ok((raster, meta))
}

// --- Helper functions ---

/// GeoKeyDirectoryTag body: version, revision, minor, key count, then key entries.
///
/// Key 1024 = GTModelTypeGeoKey (1=Projected, 2=Geographic), 2048 =
/// GeographicTypeGeoKey, 3072 = ProjectedCSTypeGeoKey.
fn geo_keys_for(epsg: u16) -> Vec<u16> {
    let model_type: u16 = if epsg < 5000 { 2 } else { 1 };
    if model_type == 2 {
        vec![
            1, 1, 0, 2, // version 1.1.0, 2 keys
            1024, 0, 1, model_type, // GTModelTypeGeoKey = Geographic
            2048, 0, 1, epsg, // GeographicTypeGeoKey
        ]
    } else {
        vec![
            1, 1, 0, 2, // version 1.1.0, 2 keys
            1024, 0, 1, model_type, // GTModelTypeGeoKey = Projected
            3072, 0, 1, epsg, // ProjectedCSTypeGeoKey
        ]
    }
}

fn first_ifd_offset(data: &[u8]) -> Result<usize, Error> {
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
    Ok(ifd_offset)
}

fn geo_meta_from(
    data: &[u8],
    pixel_scale_offset: u32,
    tiepoint_offset: u32,
    geo_keys_offset: u32,
    geo_keys_count: u32,
) -> GeoTiffMetadata {
    let pixel_width = read_f64_at(data, pixel_scale_offset as usize);
    let pixel_height = read_f64_at(data, pixel_scale_offset as usize + 8);
    let origin_x = read_f64_at(data, tiepoint_offset as usize + 24);
    let origin_y = read_f64_at(data, tiepoint_offset as usize + 32);

    let mut epsg = 0u16;
    if geo_keys_count >= 8 {
        // skip the 4-short header, then read key entries
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

    GeoTiffMetadata {
        origin_x,
        origin_y,
        pixel_width,
        pixel_height,
        epsg,
    }
}

/// Reserve space for a SHORT array and return the IFD entry's value field.
///
/// One or two shorts fit in the entry itself, so per the TIFF spec they must be
/// packed there rather than pointed at.
fn alloc_shorts(next: &mut u32, values: &[u16]) -> u32 {
    if values.len() <= 2 {
        let low = u32::from(values.first().copied().unwrap_or(0));
        let high = u32::from(values.get(1).copied().unwrap_or(0));
        return low | (high << 16);
    }
    let offset = *next;
    *next += values.len() as u32 * 2;
    offset
}

/// Write the block reserved by [`alloc_shorts`], if it needed one.
fn write_shorts_block(buf: &mut Vec<u8>, values: &[u16]) -> Result<(), Error> {
    if values.len() > 2 {
        for &val in values {
            write_u16(buf, val)?;
        }
    }
    Ok(())
}

/// Read a SHORT array from an IFD entry, inline or out-of-line.
fn read_shorts(data: &[u8], entry_off: usize) -> Result<Vec<u16>, Error> {
    let count = read_u32_at(data, entry_off + 4) as usize;
    let value = read_u32_at(data, entry_off + 8);
    if count <= 2 {
        return Ok((0..count).map(|i| (value >> (16 * i)) as u16).collect());
    }
    let offset = value as usize;
    if offset + count * 2 > data.len() {
        return Err(Error::Format("SHORT array outside file".into()));
    }
    Ok((0..count)
        .map(|i| read_u16_at(data, offset + i * 2))
        .collect())
}

/// Read StripOffsets / StripByteCounts (SHORT or LONG) from an IFD entry.
fn read_offset_array(data: &[u8], entry_off: usize) -> Result<Vec<u32>, Error> {
    let typ = read_u16_at(data, entry_off + 2);
    let count = read_u32_at(data, entry_off + 4) as usize;
    if count == 0 {
        return Ok(Vec::new());
    }
    let elem_size: usize = match typ {
        3 => 2,
        4 => 4,
        _ => {
            return Err(Error::Format(format!(
                "unsupported offset array type {typ}"
            )));
        }
    };
    let inline_n = 4 / elem_size;
    if count <= inline_n {
        let mut out = Vec::with_capacity(count);
        for i in 0..count {
            let at = entry_off + 8 + i * elem_size;
            let v = if elem_size == 2 {
                u32::from(read_u16_at(data, at))
            } else {
                read_u32_at(data, at)
            };
            out.push(v);
        }
        return Ok(out);
    }
    let offset = read_u32_at(data, entry_off + 8) as usize;
    if offset + count * elem_size > data.len() {
        return Err(Error::Format("offset array outside file".into()));
    }
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let at = offset + i * elem_size;
        let v = if elem_size == 2 {
            u32::from(read_u16_at(data, at))
        } else {
            read_u32_at(data, at)
        };
        out.push(v);
    }
    Ok(out)
}

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

    fn utm_meta() -> GeoTiffMetadata {
        GeoTiffMetadata {
            origin_x: 500000.0,
            origin_y: 4500000.0,
            pixel_width: 10.0,
            pixel_height: 10.0,
            epsg: 32632,
        }
    }

    /// Find an IFD entry and return its (count, value field).
    fn ifd_tag(data: &[u8], tag: u16) -> Option<(u32, u32)> {
        let ifd = read_u32_at(data, 4) as usize;
        let entries = read_u16_at(data, ifd) as usize;
        (0..entries)
            .map(|i| ifd + 2 + i * 12)
            .find(|&off| read_u16_at(data, off) == tag)
            .map(|off| (read_u32_at(data, off + 4), read_u32_at(data, off + 8)))
    }

    fn band_of(values: Vec<f64>) -> Raster {
        Raster::from_vec(3, 2, values, 10.0, -9999.0).unwrap()
    }

    fn assert_bands_eq(expected: &[Vec<f64>], actual: &BandedRaster) {
        assert_eq!(actual.band_count(), expected.len());
        for (b, values) in expected.iter().enumerate() {
            let band = actual.band(b).unwrap();
            for row in 0..actual.height() {
                for col in 0..actual.width() {
                    let want = values[row * actual.width() + col];
                    let got = band.get(row, col).unwrap();
                    assert_eq!(got, want, "band {b} at ({row},{col})");
                }
            }
        }
    }

    #[test]
    fn test_geotiff_bands_rgb_u8_roundtrip() {
        let bands = vec![
            vec![0.0, 255.0, 12.0, 30.0, 200.0, 7.0],
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            vec![250.0, 240.0, 230.0, 220.0, 210.0, 200.0],
        ];
        let raster = BandedRaster::with_names(
            bands.iter().cloned().map(band_of).collect(),
            vec!["red".into(), "green".into(), "blue".into()],
        )
        .unwrap();

        let mut buf = Vec::new();
        write_geotiff_bands(&raster, &utm_meta(), SampleFormat::U8, &mut buf).unwrap();

        assert_eq!(
            ifd_tag(&buf, 262).unwrap().1,
            2,
            "photometric should be RGB"
        );
        assert_eq!(ifd_tag(&buf, 284).unwrap().1, 1, "chunky planar config");
        assert_eq!(ifd_tag(&buf, 258).unwrap().0, 3, "3 BitsPerSample values");
        assert!(ifd_tag(&buf, 338).is_none(), "no extra samples for RGB");

        let (read, meta) = read_geotiff_bands(&buf).unwrap();
        assert_eq!(read.width(), 3);
        assert_eq!(read.height(), 2);
        assert_bands_eq(&bands, &read);
        assert_eq!(meta.origin_x, 500000.0);
        assert_eq!(meta.origin_y, 4500000.0);
        assert_eq!(meta.pixel_width, 10.0);
        assert_eq!(meta.pixel_height, 10.0);
        assert_eq!(meta.epsg, 32632);
        // names live in memory only
        assert_eq!(read.band_name(0), None);
    }

    #[test]
    fn test_geotiff_bands_rgba_u8_roundtrip() {
        let bands = vec![
            vec![10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
            vec![11.0, 21.0, 31.0, 41.0, 51.0, 61.0],
            vec![12.0, 22.0, 32.0, 42.0, 52.0, 62.0],
            vec![0.0, 255.0, 128.0, 64.0, 32.0, 16.0],
        ];
        let raster = BandedRaster::new(bands.iter().cloned().map(band_of).collect()).unwrap();

        let mut buf = Vec::new();
        write_geotiff_bands(&raster, &utm_meta(), SampleFormat::U8, &mut buf).unwrap();

        let (count, value) = ifd_tag(&buf, 338).expect("ExtraSamples tag");
        assert_eq!(count, 1);
        assert_eq!(value, 2, "unassociated alpha");
        assert_eq!(ifd_tag(&buf, 277).unwrap().1, 4);

        let (read, meta) = read_geotiff_bands(&buf).unwrap();
        assert_bands_eq(&bands, &read);
        assert_eq!(meta.epsg, 32632);
    }

    #[test]
    fn test_geotiff_bands_f64_roundtrip() {
        let bands = vec![
            vec![1.5, -2.25, 3.0e10, 4.125, -9999.0, 6.5],
            vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        ];
        let raster = BandedRaster::new(bands.iter().cloned().map(band_of).collect()).unwrap();

        let mut buf = Vec::new();
        write_geotiff_bands(&raster, &utm_meta(), SampleFormat::F64, &mut buf).unwrap();

        assert_eq!(
            ifd_tag(&buf, 262).unwrap().1,
            1,
            "two bands are min-is-black"
        );
        // two shorts pack into the entry itself
        assert_eq!(ifd_tag(&buf, 339).unwrap(), (2, 3 | (3 << 16)));
        let (count, value) = ifd_tag(&buf, 338).expect("ExtraSamples tag");
        assert_eq!((count, value), (1, 0), "second band is unspecified data");

        let (read, _) = read_geotiff_bands(&buf).unwrap();
        assert_bands_eq(&bands, &read);
    }

    #[test]
    fn test_geotiff_bands_u8_clamps_out_of_range() {
        let raster =
            BandedRaster::new(vec![band_of(vec![300.0, -5.0, 12.6, f64::NAN, 255.0, 0.4])])
                .unwrap();

        let mut buf = Vec::new();
        write_geotiff_bands(&raster, &utm_meta(), SampleFormat::U8, &mut buf).unwrap();

        let (read, _) = read_geotiff_bands(&buf).unwrap();
        assert_bands_eq(&[vec![255.0, 0.0, 13.0, 0.0, 255.0, 0.0]], &read);
    }

    #[test]
    fn test_geotiff_bands_reads_single_band_write_geotiff_output() {
        let values = vec![1.5, -2.25, 3.0, 4.125, -9999.0, 6.5];
        let mut buf = Vec::new();
        write_geotiff(&band_of(values.clone()), &utm_meta(), &mut buf).unwrap();

        let (read, meta) = read_geotiff_bands(&buf).unwrap();
        assert_bands_eq(&[values], &read);
        assert_eq!(meta.epsg, 32632);
    }

    #[test]
    fn sample_format_names_and_tiff_pairs_agree() {
        for format in SampleFormat::ALL {
            assert_eq!(SampleFormat::from_name(format.name()), Some(format));
            assert_eq!(
                SampleFormat::from_tiff(format.bits(), format.tag_value()),
                Some(format)
            );
        }
        assert_eq!(SampleFormat::from_name("float64"), None);
        assert_eq!(SampleFormat::from_tiff(24, 1), None);
    }

    #[test]
    fn test_geotiff_bands_roundtrip_every_format() {
        for format in SampleFormat::ALL {
            // inside every format's range, and exact in f32
            let values = vec![0.0, 12.0, 100.0, 5.0, 64.0, 1.0];
            let raster = BandedRaster::new(vec![band_of(values.clone())]).unwrap();

            let mut buf = Vec::new();
            write_geotiff_bands(&raster, &utm_meta(), format, &mut buf).unwrap();
            assert_eq!(
                ifd_tag(&buf, 258).unwrap().1,
                u32::from(format.bits()),
                "{} bits",
                format.name()
            );
            assert_eq!(
                ifd_tag(&buf, 339).unwrap().1,
                u32::from(format.tag_value()),
                "{} sample format",
                format.name()
            );

            let (read, _) = read_geotiff_bands(&buf).unwrap();
            assert_bands_eq(&[values], &read);
        }
    }

    #[test]
    fn test_geotiff_bands_signed_samples_survive() {
        let values = vec![-128.0, -1.0, 0.0, 1.0, 127.0, -64.0];
        let raster = BandedRaster::new(vec![band_of(values.clone())]).unwrap();
        let mut buf = Vec::new();
        write_geotiff_bands(&raster, &utm_meta(), SampleFormat::I8, &mut buf).unwrap();
        let (read, _) = read_geotiff_bands(&buf).unwrap();
        assert_bands_eq(&[values], &read);
    }

    #[test]
    fn test_geotiff_bands_truncated_is_error() {
        let raster = BandedRaster::new(vec![band_of(vec![1.0; 6]), band_of(vec![2.0; 6])]).unwrap();
        let mut buf = Vec::new();
        write_geotiff_bands(&raster, &utm_meta(), SampleFormat::U8, &mut buf).unwrap();
        buf.truncate(buf.len() - 4);
        assert!(read_geotiff_bands(&buf).is_err());
    }

    #[test]
    fn test_read_geotiff_u16_uncompressed_values_match() {
        let values = vec![0.0, 1.0, 255.0, 256.0, 1000.0, 65535.0];
        let raster = BandedRaster::new(vec![band_of(values.clone())]).unwrap();
        let mut buf = Vec::new();
        write_geotiff_bands(&raster, &utm_meta(), SampleFormat::U16, &mut buf).unwrap();

        let (read, meta) = read_geotiff(&buf).unwrap();
        assert_eq!(read.width(), 3);
        assert_eq!(read.height(), 2);
        assert_eq!(meta.epsg, 32632);
        for row in 0..2 {
            for col in 0..3 {
                assert_eq!(
                    read.get(row, col),
                    Some(values[row * 3 + col]),
                    "mismatch at ({row},{col})"
                );
            }
        }
    }

    #[test]
    fn test_read_geotiff_rejects_compressed() {
        let raster = Raster::from_vec(2, 2, vec![1.0; 4], 1.0, -9999.0).unwrap();
        let mut buf = Vec::new();
        write_geotiff(&raster, &utm_meta(), &mut buf).unwrap();
        let ifd = read_u32_at(&buf, 4) as usize;
        let entries = read_u16_at(&buf, ifd) as usize;
        for i in 0..entries {
            let off = ifd + 2 + i * 12;
            if read_u16_at(&buf, off) == 259 {
                buf[off + 8] = 5;
                buf[off + 9] = 0;
                buf[off + 10] = 0;
                buf[off + 11] = 0;
            }
        }
        assert!(read_geotiff(&buf).is_err());
    }

    #[test]
    fn test_read_geotiff_u16_multiple_strips() {
        let values = vec![0.0, 1.0, 255.0, 256.0, 1000.0, 65535.0];
        let raster = BandedRaster::new(vec![band_of(values.clone())]).unwrap();
        let mut buf = Vec::new();
        write_geotiff_bands(&raster, &utm_meta(), SampleFormat::U16, &mut buf).unwrap();

        let strip0 = ifd_tag(&buf, 273).unwrap().1;
        let row_bytes = 3 * 2;
        let strip1 = strip0 + row_bytes;
        let offsets_at = buf.len() as u32;
        buf.extend_from_slice(&strip0.to_le_bytes());
        buf.extend_from_slice(&strip1.to_le_bytes());
        let counts_at = buf.len() as u32;
        buf.extend_from_slice(&row_bytes.to_le_bytes());
        buf.extend_from_slice(&row_bytes.to_le_bytes());

        let ifd = read_u32_at(&buf, 4) as usize;
        let entries = read_u16_at(&buf, ifd) as usize;
        for i in 0..entries {
            let off = ifd + 2 + i * 12;
            match read_u16_at(&buf, off) {
                273 => {
                    buf[off + 4..off + 8].copy_from_slice(&2u32.to_le_bytes());
                    buf[off + 8..off + 12].copy_from_slice(&offsets_at.to_le_bytes());
                }
                278 => {
                    buf[off + 8..off + 12].copy_from_slice(&1u32.to_le_bytes());
                }
                279 => {
                    buf[off + 4..off + 8].copy_from_slice(&2u32.to_le_bytes());
                    buf[off + 8..off + 12].copy_from_slice(&counts_at.to_le_bytes());
                }
                _ => {}
            }
        }

        let (read, _) = read_geotiff(&buf).unwrap();
        for row in 0..2 {
            for col in 0..3 {
                assert_eq!(read.get(row, col), Some(values[row * 3 + col]));
            }
        }
    }
}
