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
    pub simple_packing: Option<SimplePacking>,
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

/// Simple packing parameters from Data Representation Template 5.0.
#[derive(Debug, Clone, Copy)]
pub struct SimplePacking {
    pub reference_value: f32,
    pub binary_scale: i16,
    pub decimal_scale: i16,
    pub bits_per_value: u8,
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

    let data = match (message.packing, message.simple_packing) {
        (PackingType::Simple, Some(p)) => decode_simple_packing(
            &raw,
            width * height,
            p.reference_value,
            p.binary_scale,
            p.decimal_scale,
            usize::from(p.bits_per_value),
        )?,
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
    msg_start: u64,
    discipline: u8,
    edition: u8,
    total_length: u64,
) -> Result<GribMessage, Error> {
    let msg_end = msg_start + total_length;
    let mut reference_time = GribTime {
        year: 0,
        month: 0,
        day: 0,
        hour: 0,
        minute: 0,
        second: 0,
    };
    let mut grid_definition = GridDefinition {
        template: GridTemplate::Unknown(0),
        ni: 0,
        nj: 0,
        lat_first: 0.0,
        lon_first: 0.0,
        lat_last: 0.0,
        lon_last: 0.0,
        di: 0.0,
        dj: 0.0,
        scan_mode: 0,
    };
    let mut product_definition = ProductDefinition {
        parameter_category: 0,
        parameter_number: 0,
        generating_process: 0,
        forecast_time: 0,
        time_unit: TimeUnit::Hour,
        level_type: 0,
        level_value: 0.0,
    };
    let mut packing = PackingType::Unknown(0);
    let mut simple_packing = None;
    let mut data_offset = 0u64;
    let mut data_length = 0usize;
    let mut have_grid = false;
    let mut have_data = false;

    while reader.stream_position()? + 4 <= msg_end {
        let mut hdr = [0u8; 4];
        reader.read_exact(&mut hdr)?;
        if &hdr == b"7777" {
            break;
        }
        let sec_len = u32::from_be_bytes(hdr);
        if sec_len < 5 {
            return Err(Error::Format("GRIB2 section length < 5".into()));
        }
        let sec_end = reader.stream_position()? - 4 + u64::from(sec_len);
        if sec_end > msg_end {
            return Err(Error::Format("GRIB2 section overruns message".into()));
        }
        let sec_num = read_u8(reader)?;
        let body_len = sec_len as usize - 5;
        let body_start = reader.stream_position()?;
        let mut body = vec![0u8; body_len];
        reader.read_exact(&mut body)?;

        match sec_num {
            1 => reference_time = parse_identification(&body)?,
            3 => {
                grid_definition = parse_grid_definition(&body)?;
                have_grid = true;
            }
            4 => product_definition = parse_product_definition(&body)?,
            5 => {
                let (p, simple) = parse_data_representation(&body)?;
                packing = p;
                simple_packing = simple;
            }
            6 | 2 => {}
            7 => {
                data_offset = body_start;
                data_length = body_len;
                have_data = true;
            }
            _ => {}
        }
    }

    if !have_grid || grid_definition.ni == 0 || grid_definition.nj == 0 {
        return Err(Error::Format(
            "GRIB2 message missing grid definition".into(),
        ));
    }
    if !have_data {
        return Err(Error::Format("GRIB2 message missing data section".into()));
    }

    Ok(GribMessage {
        discipline,
        edition,
        reference_time,
        grid_definition,
        product_definition,
        data_offset,
        data_length,
        packing,
        simple_packing,
    })
}

fn parse_identification(body: &[u8]) -> Result<GribTime, Error> {
    if body.len() < 14 {
        return Err(Error::Format(
            "GRIB2 identification section too short".into(),
        ));
    }
    Ok(GribTime {
        year: u16_be(body, 7),
        month: body[9],
        day: body[10],
        hour: body[11],
        minute: body[12],
        second: body[13],
    })
}

fn parse_grid_definition(body: &[u8]) -> Result<GridDefinition, Error> {
    if body.len() < 70 {
        return Err(Error::Format(
            "GRIB2 grid definition section too short".into(),
        ));
    }
    let tmpl = u16_be(body, 10);
    let scale = coord_scale(u32_be(body, 36), u32_be(body, 40));
    Ok(GridDefinition {
        template: grid_template_from(tmpl),
        ni: u32_be(body, 28),
        nj: u32_be(body, 32),
        lat_first: f64::from(i32_be(body, 44)) * scale,
        lon_first: f64::from(i32_be(body, 48)) * scale,
        lat_last: f64::from(i32_be(body, 53)) * scale,
        lon_last: f64::from(i32_be(body, 57)) * scale,
        di: f64::from(u32_be(body, 61)) * scale,
        dj: f64::from(u32_be(body, 65)) * scale,
        scan_mode: body[69],
    })
}

fn parse_product_definition(body: &[u8]) -> Result<ProductDefinition, Error> {
    if body.len() < 6 {
        return Err(Error::Format(
            "GRIB2 product definition section too short".into(),
        ));
    }
    let mut product = ProductDefinition {
        parameter_category: body[4],
        parameter_number: body[5],
        generating_process: if body.len() > 6 { body[6] } else { 0 },
        forecast_time: 0,
        time_unit: TimeUnit::Hour,
        level_type: 0,
        level_value: 0.0,
    };
    if body.len() >= 23 {
        product.time_unit = time_unit_from(body[12]);
        product.forecast_time = u32_be(body, 13);
        product.level_type = body[17];
        let scale_factor = body[18] as i8;
        product.level_value = f64::from(i32_be(body, 19)) * 10.0_f64.powi(-i32::from(scale_factor));
    }
    Ok(product)
}

fn parse_data_representation(body: &[u8]) -> Result<(PackingType, Option<SimplePacking>), Error> {
    if body.len() < 15 {
        return Err(Error::Format(
            "GRIB2 data representation section too short".into(),
        ));
    }
    let tmpl = u16_be(body, 4);
    let packing = packing_from(tmpl);
    let simple = match packing {
        PackingType::Simple => Some(SimplePacking {
            reference_value: f32::from_be_bytes(body[6..10].try_into().unwrap_or([0; 4])),
            binary_scale: i16_be(body, 10),
            decimal_scale: i16_be(body, 12),
            bits_per_value: body[14],
        }),
        _ => None,
    };
    Ok((packing, simple))
}

fn grid_template_from(n: u16) -> GridTemplate {
    match n {
        0 => GridTemplate::LatLon,
        1 => GridTemplate::RotatedLatLon,
        10 => GridTemplate::Mercator,
        20 => GridTemplate::PolarStereographic,
        30 => GridTemplate::LambertConformal,
        40 => GridTemplate::GaussianLatLon,
        other => GridTemplate::Unknown(other),
    }
}

fn time_unit_from(n: u8) -> TimeUnit {
    match n {
        0 => TimeUnit::Minute,
        1 => TimeUnit::Hour,
        2 => TimeUnit::Day,
        3 => TimeUnit::Month,
        4 => TimeUnit::Year,
        other => TimeUnit::Unknown(other),
    }
}

fn packing_from(n: u16) -> PackingType {
    match n {
        0 => PackingType::Simple,
        2 => PackingType::ComplexPacking,
        3 => PackingType::ComplexPackingWithSpatialDiff,
        40 => PackingType::Jpeg2000,
        41 => PackingType::Png,
        other => PackingType::Unknown(other),
    }
}

fn coord_scale(basic_angle: u32, subdivisions: u32) -> f64 {
    if basic_angle == 0 || subdivisions == 0 {
        1e-6
    } else {
        f64::from(basic_angle) / f64::from(subdivisions)
    }
}

fn decode_simple_packing(
    packed: &[u8],
    num_points: usize,
    reference_value: f32,
    binary_scale: i16,
    decimal_scale: i16,
    bits_per_value: usize,
) -> Result<Vec<f64>, Error> {
    if bits_per_value == 0 {
        return Ok(vec![f64::from(reference_value); num_points]);
    }

    let bs = f64::from(2.0_f32.powi(i32::from(binary_scale)));
    let ds = 10.0_f64.powi(-i32::from(decimal_scale));
    let r = f64::from(reference_value);

    let mut values = Vec::with_capacity(num_points);

    for i in 0..num_points {
        let bit_offset = i * bits_per_value;
        let byte_offset = bit_offset / 8;
        if byte_offset >= packed.len() {
            values.push(f64::NAN);
            continue;
        }

        let mut raw: u64 = 0;
        let mut bits_remaining = bits_per_value;
        let mut current_byte = byte_offset;
        let mut bit_in_byte = bit_offset % 8;

        while bits_remaining > 0 && current_byte < packed.len() {
            let available = 8 - bit_in_byte;
            let take = bits_remaining.min(available);
            let mask = ((1u64 << take) - 1) as u8;
            let shift = available - take;
            let bits = (packed[current_byte] >> shift) & mask;
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

fn read_u8<R: Read>(reader: &mut R) -> Result<u8, Error> {
    let mut buf = [0u8; 1];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn u16_be(data: &[u8], off: usize) -> u16 {
    u16::from_be_bytes([data[off], data[off + 1]])
}

fn i16_be(data: &[u8], off: usize) -> i16 {
    i16::from_be_bytes([data[off], data[off + 1]])
}

fn u32_be(data: &[u8], off: usize) -> u32 {
    u32::from_be_bytes(data[off..off + 4].try_into().unwrap_or([0; 4]))
}

fn i32_be(data: &[u8], off: usize) -> i32 {
    i32::from_be_bytes(data[off..off + 4].try_into().unwrap_or([0; 4]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

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
        let values = decode_simple_packing(&[], 10, 293.15, 0, 0, 0).unwrap();
        assert_eq!(values.len(), 10);
        assert!((values[0] - 293.15).abs() < 0.01);
    }

    #[test]
    fn test_simple_packing_byte_values() {
        let packed = vec![0u8, 10, 20, 30];
        let values = decode_simple_packing(&packed, 4, 10.0, 0, 0, 8).unwrap();
        assert_eq!(values.len(), 4);
        assert!((values[0] - 10.0).abs() < 1e-6);
        assert!((values[1] - 20.0).abs() < 1e-6);
        assert!((values[2] - 30.0).abs() < 1e-6);
        assert!((values[3] - 40.0).abs() < 1e-6);
    }

    fn grib_section(num: u8, body: &[u8]) -> Vec<u8> {
        let len = 5 + body.len() as u32;
        let mut s = Vec::with_capacity(len as usize);
        s.extend_from_slice(&len.to_be_bytes());
        s.push(num);
        s.extend_from_slice(body);
        s
    }

    fn minimal_grib2() -> Vec<u8> {
        let mut sec1 = Vec::new();
        sec1.extend_from_slice(&0u16.to_be_bytes());
        sec1.extend_from_slice(&0u16.to_be_bytes());
        sec1.push(2);
        sec1.push(0);
        sec1.push(1);
        sec1.extend_from_slice(&2024u16.to_be_bytes());
        sec1.push(1);
        sec1.push(15);
        sec1.push(12);
        sec1.push(0);
        sec1.push(0);
        sec1.push(0);
        sec1.push(1);

        let mut sec3 = Vec::new();
        sec3.extend_from_slice(&0u32.to_be_bytes());
        sec3.extend_from_slice(&4u32.to_be_bytes());
        sec3.push(0);
        sec3.push(0);
        sec3.extend_from_slice(&0u16.to_be_bytes());
        sec3.push(6);
        sec3.push(0);
        sec3.extend_from_slice(&0u32.to_be_bytes());
        sec3.push(0);
        sec3.extend_from_slice(&0u32.to_be_bytes());
        sec3.push(0);
        sec3.extend_from_slice(&0u32.to_be_bytes());
        sec3.extend_from_slice(&2u32.to_be_bytes());
        sec3.extend_from_slice(&2u32.to_be_bytes());
        sec3.extend_from_slice(&0u32.to_be_bytes());
        sec3.extend_from_slice(&0u32.to_be_bytes());
        sec3.extend_from_slice(&1_000_000i32.to_be_bytes());
        sec3.extend_from_slice(&0i32.to_be_bytes());
        sec3.push(0x30);
        sec3.extend_from_slice(&0i32.to_be_bytes());
        sec3.extend_from_slice(&1_000_000i32.to_be_bytes());
        sec3.extend_from_slice(&1_000_000u32.to_be_bytes());
        sec3.extend_from_slice(&1_000_000u32.to_be_bytes());
        sec3.push(0);

        let mut sec4 = Vec::new();
        sec4.extend_from_slice(&0u16.to_be_bytes());
        sec4.extend_from_slice(&0u16.to_be_bytes());
        sec4.push(0);
        sec4.push(0);
        sec4.push(2);
        sec4.push(255);
        sec4.push(255);
        sec4.extend_from_slice(&0u16.to_be_bytes());
        sec4.push(0);
        sec4.push(1);
        sec4.extend_from_slice(&0u32.to_be_bytes());
        sec4.push(1);
        sec4.push(0);
        sec4.extend_from_slice(&0u32.to_be_bytes());
        sec4.push(255);
        sec4.push(255);
        sec4.extend_from_slice(&0u32.to_be_bytes());

        let mut sec5 = Vec::new();
        sec5.extend_from_slice(&4u32.to_be_bytes());
        sec5.extend_from_slice(&0u16.to_be_bytes());
        sec5.extend_from_slice(&10.0f32.to_be_bytes());
        sec5.extend_from_slice(&0i16.to_be_bytes());
        sec5.extend_from_slice(&0i16.to_be_bytes());
        sec5.push(8);
        sec5.push(0);

        let mut msg = Vec::new();
        msg.extend_from_slice(b"GRIB");
        msg.extend_from_slice(&[0, 0]);
        msg.push(0);
        msg.push(2);
        msg.extend_from_slice(&[0u8; 8]);
        msg.extend_from_slice(&grib_section(1, &sec1));
        msg.extend_from_slice(&grib_section(3, &sec3));
        msg.extend_from_slice(&grib_section(4, &sec4));
        msg.extend_from_slice(&grib_section(5, &sec5));
        msg.extend_from_slice(&grib_section(6, &[255]));
        msg.extend_from_slice(&grib_section(7, &[0, 10, 20, 30]));
        msg.extend_from_slice(b"7777");
        let total = msg.len() as u64;
        msg[8..16].copy_from_slice(&total.to_be_bytes());
        msg
    }

    #[test]
    fn test_scan_and_decode_minimal_grib2() {
        let bytes = minimal_grib2();
        let mut reader = BufReader::new(Cursor::new(bytes));
        let messages = scan_grib(&mut reader).unwrap();
        assert_eq!(messages.len(), 1);
        let msg = &messages[0];
        assert_eq!(msg.edition, 2);
        assert_eq!(msg.discipline, 0);
        assert_eq!(msg.reference_time.year, 2024);
        assert_eq!(msg.reference_time.month, 1);
        assert_eq!(msg.reference_time.day, 15);
        assert_eq!(msg.grid_definition.template, GridTemplate::LatLon);
        assert_eq!(msg.grid_definition.ni, 2);
        assert_eq!(msg.grid_definition.nj, 2);
        assert!((msg.grid_definition.lat_first - 1.0).abs() < 1e-9);
        assert!((msg.grid_definition.di - 1.0).abs() < 1e-9);
        assert_eq!(msg.packing, PackingType::Simple);
        assert!(msg.data_offset > 0);
        assert_eq!(msg.data_length, 4);

        let raster = decode_grib_message(&mut reader, msg).unwrap();
        assert_eq!(raster.width(), 2);
        assert_eq!(raster.height(), 2);
        assert!((raster.get(0, 0).unwrap() - 10.0).abs() < 1e-4);
        assert!((raster.get(0, 1).unwrap() - 20.0).abs() < 1e-4);
        assert!((raster.get(1, 0).unwrap() - 30.0).abs() < 1e-4);
        assert!((raster.get(1, 1).unwrap() - 40.0).abs() < 1e-4);
    }
}
