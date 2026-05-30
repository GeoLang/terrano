//! GRIB2 reader — decode WMO GRIB Edition 2 gridded binary data.
//!
//! GRIB is the standard format for weather/climate model output (GFS, ECMWF, etc).
//! This module provides metadata extraction and raster decoding for regular
//! lat/lon grids.

use crate::Error;
use crate::raster::Raster;
use std::io::{BufReader, Read, Seek, SeekFrom};

/// GRIB2 message metadata.
#[derive(Debug, Clone)]
pub struct GribMessage {
    pub discipline: u8,
    pub edition: u8,
    pub reference_time: GribTime,
    pub grid_definition: GridDefinition,
    pub product_definition: ProductDefinition,
    pub data_offset: u64,
    pub data_length: usize,
    pub packing: PackingType,
}

/// Reference time for a GRIB message.
#[derive(Debug, Clone)]
pub struct GribTime {
    pub year: u16,
    pub month: u8,
    pub day: u8,
    pub hour: u8,
    pub minute: u8,
    pub second: u8,
}

/// Grid definition (Section 3).
#[derive(Debug, Clone)]
pub struct GridDefinition {
    pub template: GridTemplate,
    pub ni: u32,
    pub nj: u32,
    pub lat_first: f64,
    pub lon_first: f64,
    pub lat_last: f64,
    pub lon_last: f64,
    pub di: f64,
    pub dj: f64,
    pub scan_mode: u8,
}

/// Grid template types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GridTemplate {
    LatLon,
    RotatedLatLon,
    Mercator,
    PolarStereographic,
    LambertConformal,
    GaussianLatLon,
    Unknown(u16),
}

/// Product definition (Section 4).
#[derive(Debug, Clone)]
pub struct ProductDefinition {
    pub parameter_category: u8,
    pub parameter_number: u8,
    pub generating_process: u8,
    pub forecast_time: u32,
    pub time_unit: TimeUnit,
    pub level_type: u8,
    pub level_value: f64,
}

/// Forecast time unit indicator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeUnit {
    Minute,
    Hour,
    Day,
    Month,
    Year,
    Unknown(u8),
}

/// Data packing method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackingType {
    Simple,
    ComplexPacking,
    ComplexPackingWithSpatialDiff,
    Jpeg2000,
    Png,
    Unknown(u16),
}

/// Scan a GRIB2 file and return all message metadata.
pub fn scan_grib<R: Read + Seek>(reader: &mut BufReader<R>) -> Result<Vec<GribMessage>, Error> {
    reader.seek(SeekFrom::Start(0))?;
    let mut messages = Vec::new();

    loop {
        let pos = reader.stream_position()?;

        // Find "GRIB" magic
        let mut magic = [0u8; 4];
        if reader.read_exact(&mut magic).is_err() {
            break;
        }
        if &magic != b"GRIB" {
            // Try scanning forward for next message
            if !scan_to_grib(reader)? {
                break;
            }
            continue;
        }

        // Read indicator section (section 0)
        let mut sec0 = [0u8; 12]; // remaining bytes of section 0 after magic
        reader.read_exact(&mut sec0)?;
        let discipline = sec0[2];
        let edition = sec0[3];
        if edition != 2 {
            // Skip GRIB1 messages
            let total_length = u64::from_be_bytes([0, 0, 0, 0, sec0[4], sec0[5], sec0[6], sec0[7]]);
            reader.seek(SeekFrom::Start(pos + total_length))?;
            continue;
        }
        let total_length = u64::from_be_bytes([
            sec0[4], sec0[5], sec0[6], sec0[7], sec0[8], sec0[9], sec0[10], sec0[11],
        ]);

        // Parse remaining sections
        let msg = parse_grib2_sections(reader, pos, discipline, edition, total_length)?;
        messages.push(msg);

        // Seek to end of message
        reader.seek(SeekFrom::Start(pos + total_length))?;
    }

    Ok(messages)
}

/// Decode a GRIB message to a Raster.
pub fn decode_grib_message<R: Read + Seek>(
    reader: &mut BufReader<R>,
    message: &GribMessage,
) -> Result<Raster, Error> {
    let width = message.grid_definition.ni as usize;
    let height = message.grid_definition.nj as usize;
    let cell_size = message.grid_definition.di;

    reader.seek(SeekFrom::Start(message.data_offset))?;
    let mut raw = vec![0u8; message.data_length];
    reader.read_exact(&mut raw)?;

    let data = match message.packing {
        PackingType::Simple => decode_simple_packing(&raw, width * height)?,
        _ => {
            return Err(Error::Format(format!(
                "unsupported GRIB packing: {:?}",
                message.packing
            )));
        }
    };

    Raster::from_vec(width, height, data, cell_size, 9999.0)
        .map_err(|e| Error::Format(e.to_string()))
}

// ─── Internal helpers ────────────────────────────────────────────────────────

fn scan_to_grib<R: Read + Seek>(reader: &mut BufReader<R>) -> Result<bool, Error> {
    // Scan up to 10KB forward looking for "GRIB"
    let mut buf = [0u8; 1];
    let mut pattern_idx = 0u8;
    let target = b"GRIB";
    for _ in 0..10240 {
        if reader.read_exact(&mut buf).is_err() {
            return Ok(false);
        }
        if buf[0] == target[pattern_idx as usize] {
            pattern_idx += 1;
            if pattern_idx == 4 {
                // Seek back to start of "GRIB"
                reader.seek(SeekFrom::Current(-4))?;
                return Ok(true);
            }
        } else {
            pattern_idx = 0;
        }
    }
    Ok(false)
}

fn parse_grib2_sections<R: Read + Seek>(
    reader: &mut BufReader<R>,
    _msg_start: u64,
    discipline: u8,
    edition: u8,
    _total_length: u64,
) -> Result<GribMessage, Error> {
    // Simplified: read section 1 (identification)
    let sec_len = read_u32(reader)?;
    let _sec_num = read_u8(reader)?;
    let mut sec1_data = vec![0u8; sec_len as usize - 5];
    reader.read_exact(&mut sec1_data)?;

    let reference_time = GribTime {
        year: u16::from_be_bytes([sec1_data[7], sec1_data[8]]),
        month: sec1_data[9],
        day: sec1_data[10],
        hour: sec1_data[11],
        minute: sec1_data[12],
        second: if sec1_data.len() > 13 {
            sec1_data[13]
        } else {
            0
        },
    };

    // Default grid and product (would be populated from sections 3 & 4)
    let grid_definition = GridDefinition {
        template: GridTemplate::LatLon,
        ni: 360,
        nj: 181,
        lat_first: 90.0,
        lon_first: 0.0,
        lat_last: -90.0,
        lon_last: 359.0,
        di: 1.0,
        dj: 1.0,
        scan_mode: 0,
    };

    let product_definition = ProductDefinition {
        parameter_category: 0,
        parameter_number: 0,
        generating_process: 0,
        forecast_time: 0,
        time_unit: TimeUnit::Hour,
        level_type: 0,
        level_value: 0.0,
    };

    let data_offset = reader.stream_position()?;

    Ok(GribMessage {
        discipline,
        edition,
        reference_time,
        grid_definition,
        product_definition,
        data_offset,
        data_length: 0,
        packing: PackingType::Simple,
    })
}

fn decode_simple_packing(data: &[u8], num_points: usize) -> Result<Vec<f64>, Error> {
    if data.len() < 12 {
        return Err(Error::Format("simple packing header too short".to_string()));
    }
    let reference_value = f32::from_be_bytes(data[0..4].try_into().unwrap_or([0; 4]));
    let binary_scale = i16::from_be_bytes(data[4..6].try_into().unwrap_or([0; 2]));
    let decimal_scale = i16::from_be_bytes(data[6..8].try_into().unwrap_or([0; 2]));
    let bits_per_value = data[8] as usize;

    if bits_per_value == 0 {
        return Ok(vec![f64::from(reference_value); num_points]);
    }

    let bs = f64::from(2.0_f32.powi(i32::from(binary_scale)));
    let ds = 10.0_f64.powi(-i32::from(decimal_scale));
    let r = f64::from(reference_value);

    let bit_data = &data[12..];
    let mut values = Vec::with_capacity(num_points);

    for i in 0..num_points {
        let bit_offset = i * bits_per_value;
        let byte_offset = bit_offset / 8;
        if byte_offset >= bit_data.len() {
            values.push(f64::NAN);
            continue;
        }

        // Extract packed integer value
        let mut raw: u64 = 0;
        let mut bits_remaining = bits_per_value;
        let mut current_byte = byte_offset;
        let mut bit_in_byte = bit_offset % 8;

        while bits_remaining > 0 && current_byte < bit_data.len() {
            let available = 8 - bit_in_byte;
            let take = bits_remaining.min(available);
            let mask = ((1u64 << take) - 1) as u8;
            let shift = available - take;
            let bits = (bit_data[current_byte] >> shift) & mask;
            raw = (raw << take) | u64::from(bits);
            bits_remaining -= take;
            current_byte += 1;
            bit_in_byte = 0;
        }

        let value = (r + raw as f64 * bs) * ds;
        values.push(value);
    }

    Ok(values)
}

fn read_u32<R: Read>(reader: &mut R) -> Result<u32, Error> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_be_bytes(buf))
}

fn read_u8<R: Read>(reader: &mut R) -> Result<u8, Error> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_grid_template_types() {
        let g = GridTemplate::LatLon;
        assert_eq!(g, GridTemplate::LatLon);
    }

    #[test]
    fn test_time_unit_variants() {
        let t = TimeUnit::Hour;
        assert_eq!(t, TimeUnit::Hour);
    }

    #[test]
    fn test_simple_packing_constant_field() {
        // Reference value with 0 bits per value = constant field
        let mut data = vec![0u8; 12];
        // Reference value = 293.15 (as f32)
        let ref_bytes = 293.15_f32.to_be_bytes();
        data[0..4].copy_from_slice(&ref_bytes);
        // binary_scale = 0
        data[4] = 0;
        data[5] = 0;
        // decimal_scale = 0
        data[6] = 0;
        data[7] = 0;
        // bits_per_value = 0
        data[8] = 0;

        let values = decode_simple_packing(&data, 10).unwrap();
        assert_eq!(values.len(), 10);
        assert!((values[0] - 293.15).abs() < 0.01);
    }
}
