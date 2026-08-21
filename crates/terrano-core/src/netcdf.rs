//! NetCDF reader — extract raster data from NetCDF-4 / NetCDF Classic files.
//!
//! Provides a lightweight reader for structured gridded datasets commonly
//! used in climate science, oceanography, and meteorology.

use crate::Error;
use crate::raster::Raster;
use std::collections::HashMap;
use std::io::{BufReader, Read, Seek, SeekFrom};

/// Metadata from a NetCDF file.
#[derive(Debug, Clone)]
pub struct NetCdfMetadata {
    pub dimensions: Vec<Dimension>,
    pub variables: Vec<Variable>,
    pub global_attributes: HashMap<String, AttributeValue>,
}

/// A NetCDF dimension.
#[derive(Debug, Clone)]
pub struct Dimension {
    pub name: String,
    pub size: usize,
    pub is_unlimited: bool,
}

/// A NetCDF variable.
#[derive(Debug, Clone)]
pub struct Variable {
    pub name: String,
    pub data_type: DataType,
    pub dimensions: Vec<String>,
    pub attributes: HashMap<String, AttributeValue>,
    pub shape: Vec<usize>,
    /// Byte offset of this variable's data from the start of the file (`begin`).
    pub begin: u64,
}

/// NetCDF data types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    Byte,
    Short,
    Int,
    Float,
    Double,
    Char,
}

/// Attribute value variants.
#[derive(Debug, Clone)]
pub enum AttributeValue {
    Text(String),
    Int(Vec<i32>),
    Float(Vec<f32>),
    Double(Vec<f64>),
    Short(Vec<i16>),
    Byte(Vec<u8>),
}

impl DataType {
    fn size_bytes(self) -> usize {
        match self {
            Self::Byte | Self::Char => 1,
            Self::Short => 2,
            Self::Int | Self::Float => 4,
            Self::Double => 8,
        }
    }
}

/// Read NetCDF metadata from a reader.
pub fn read_netcdf_metadata<R: Read + Seek>(
    reader: &mut BufReader<R>,
) -> Result<NetCdfMetadata, Error> {
    reader.seek(SeekFrom::Start(0))?;
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;

    // Check magic: "CDF\x01" (classic) or "CDF\x02" (64-bit offset)
    if &magic[0..3] != b"CDF" || (magic[3] != 1 && magic[3] != 2) {
        return Err(Error::Format("not a valid NetCDF classic file".to_string()));
    }

    let version = magic[3];

    // Read number of records
    let num_recs = read_u32(reader)?;
    let _ = num_recs; // unlimited dimension length

    // Read dimensions
    let dimensions = read_dim_list(reader)?;

    // Read global attributes
    let global_attributes = read_att_list(reader)?;

    // Read variables
    let variables = read_var_list(reader, &dimensions, version)?;

    Ok(NetCdfMetadata {
        dimensions,
        variables,
        global_attributes,
    })
}

/// Read a 2D variable as a Raster.
pub fn read_netcdf_variable<R: Read + Seek>(
    reader: &mut BufReader<R>,
    metadata: &NetCdfMetadata,
    variable_name: &str,
) -> Result<Raster, Error> {
    let var = metadata
        .variables
        .iter()
        .find(|v| v.name == variable_name)
        .ok_or_else(|| Error::Format(format!("variable '{variable_name}' not found")))?;

    if var.shape.len() < 2 {
        return Err(Error::Format(format!(
            "variable '{variable_name}' has fewer than 2 dimensions"
        )));
    }

    let height = var.shape[var.shape.len() - 2];
    let width = var.shape[var.shape.len() - 1];

    // Determine scale/offset from CF conventions
    let scale_factor = get_attr_f64(&var.attributes, "scale_factor").unwrap_or(1.0);
    let add_offset = get_attr_f64(&var.attributes, "add_offset").unwrap_or(0.0);
    let fill_value = get_attr_f64(&var.attributes, "_FillValue").unwrap_or(f64::NAN);

    // Cell size from coordinate variables (best effort)
    let cell_size = infer_cell_size(metadata, var);

    // Read raw data (last 2D slice if >2D)
    let total_2d = width * height;
    let byte_size = var.data_type.size_bytes();
    let slice_offset = if var.shape.len() > 2 {
        // Read last time step: skip (product of leading dims - 1) * 2d_size
        let leading: usize = var.shape[..var.shape.len() - 2].iter().product();
        (leading - 1) * total_2d * byte_size
    } else {
        0
    };

    reader.seek(SeekFrom::Start(var.begin + slice_offset as u64))?;

    let mut raw = vec![0u8; total_2d * byte_size];
    reader.read_exact(&mut raw)?;

    let data: Vec<f64> = (0..total_2d)
        .map(|i| {
            let offset = i * byte_size;
            let raw_val = match var.data_type {
                DataType::Float => f64::from(f32::from_be_bytes(
                    raw[offset..offset + 4].try_into().unwrap_or([0; 4]),
                )),
                DataType::Double => {
                    f64::from_be_bytes(raw[offset..offset + 8].try_into().unwrap_or([0; 8]))
                }
                DataType::Short => f64::from(i16::from_be_bytes(
                    raw[offset..offset + 2].try_into().unwrap_or([0; 2]),
                )),
                DataType::Int => f64::from(i32::from_be_bytes(
                    raw[offset..offset + 4].try_into().unwrap_or([0; 4]),
                )),
                DataType::Byte => f64::from(raw[offset]),
                DataType::Char => f64::from(raw[offset]),
            };
            if (raw_val - fill_value).abs() < f64::EPSILON || raw_val.is_nan() {
                f64::NAN
            } else {
                raw_val * scale_factor + add_offset
            }
        })
        .collect();

    Raster::from_vec(width, height, data, cell_size, f64::NAN)
}

// ─── Internal helpers ────────────────────────────────────────────────────────

fn read_u32<R: Read>(reader: &mut R) -> Result<u32, Error> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_be_bytes(buf))
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64, Error> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_be_bytes(buf))
}

fn read_dim_list<R: Read>(reader: &mut R) -> Result<Vec<Dimension>, Error> {
    let tag = read_u32(reader)?;
    let count = read_u32(reader)?;
    if tag == 0 && count == 0 {
        return Ok(Vec::new());
    }
    // tag == 0x0000000A => NC_DIMENSION
    let mut dims = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let name = read_name(reader)?;
        let size = read_u32(reader)? as usize;
        dims.push(Dimension {
            name,
            size,
            is_unlimited: size == 0,
        });
    }
    Ok(dims)
}

fn read_att_list<R: Read>(reader: &mut R) -> Result<HashMap<String, AttributeValue>, Error> {
    let tag = read_u32(reader)?;
    let count = read_u32(reader)?;
    if tag == 0 && count == 0 {
        return Ok(HashMap::new());
    }
    let mut attrs = HashMap::with_capacity(count as usize);
    for _ in 0..count {
        let name = read_name(reader)?;
        let nc_type = read_u32(reader)?;
        let nelems = read_u32(reader)? as usize;
        let value = read_att_value(reader, nc_type, nelems)?;
        attrs.insert(name, value);
    }
    Ok(attrs)
}

fn read_var_list<R: Read>(
    reader: &mut R,
    dims: &[Dimension],
    version: u8,
) -> Result<Vec<Variable>, Error> {
    let tag = read_u32(reader)?;
    let count = read_u32(reader)?;
    if tag == 0 && count == 0 {
        return Ok(Vec::new());
    }
    let mut vars = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let name = read_name(reader)?;
        let ndims = read_u32(reader)? as usize;
        let mut dim_names = Vec::with_capacity(ndims);
        let mut shape = Vec::with_capacity(ndims);
        for _ in 0..ndims {
            let dim_id = read_u32(reader)? as usize;
            if dim_id < dims.len() {
                dim_names.push(dims[dim_id].name.clone());
                shape.push(dims[dim_id].size);
            }
        }
        let attributes = read_att_list(reader)?;
        let nc_type = read_u32(reader)?;
        let _vsize = read_u32(reader)?;
        let begin = if version == 2 {
            read_u64(reader)?
        } else {
            u64::from(read_u32(reader)?)
        };
        let data_type = match nc_type {
            1 => DataType::Byte,
            2 => DataType::Char,
            3 => DataType::Short,
            4 => DataType::Int,
            5 => DataType::Float,
            6 => DataType::Double,
            _ => DataType::Float,
        };
        vars.push(Variable {
            name,
            data_type,
            dimensions: dim_names,
            attributes,
            shape,
            begin,
        });
    }
    Ok(vars)
}

fn read_name<R: Read>(reader: &mut R) -> Result<String, Error> {
    let len = read_u32(reader)? as usize;
    let padded = (len + 3) & !3; // 4-byte aligned
    let mut buf = vec![0u8; padded];
    reader.read_exact(&mut buf)?;
    buf.truncate(len);
    String::from_utf8(buf).map_err(|e| Error::Format(e.to_string()))
}

fn read_att_value<R: Read>(
    reader: &mut R,
    nc_type: u32,
    nelems: usize,
) -> Result<AttributeValue, Error> {
    let byte_len = match nc_type {
        1 => nelems,
        2 => nelems,
        3 => nelems * 2,
        4 => nelems * 4,
        5 => nelems * 4,
        6 => nelems * 8,
        _ => nelems,
    };
    let padded = (byte_len + 3) & !3;
    let mut buf = vec![0u8; padded];
    reader.read_exact(&mut buf)?;
    buf.truncate(byte_len);

    let val = match nc_type {
        1 => AttributeValue::Byte(buf),
        2 => {
            let s = String::from_utf8_lossy(&buf)
                .trim_end_matches('\0')
                .to_string();
            AttributeValue::Text(s)
        }
        3 => {
            let vals: Vec<i16> = buf
                .chunks_exact(2)
                .map(|c| i16::from_be_bytes(c.try_into().unwrap_or([0; 2])))
                .collect();
            AttributeValue::Short(vals)
        }
        4 => {
            let vals: Vec<i32> = buf
                .chunks_exact(4)
                .map(|c| i32::from_be_bytes(c.try_into().unwrap_or([0; 4])))
                .collect();
            AttributeValue::Int(vals)
        }
        5 => {
            let vals: Vec<f32> = buf
                .chunks_exact(4)
                .map(|c| f32::from_be_bytes(c.try_into().unwrap_or([0; 4])))
                .collect();
            AttributeValue::Float(vals)
        }
        6 => {
            let vals: Vec<f64> = buf
                .chunks_exact(8)
                .map(|c| f64::from_be_bytes(c.try_into().unwrap_or([0; 8])))
                .collect();
            AttributeValue::Double(vals)
        }
        _ => AttributeValue::Byte(buf),
    };
    Ok(val)
}

fn get_attr_f64(attrs: &HashMap<String, AttributeValue>, key: &str) -> Option<f64> {
    match attrs.get(key)? {
        AttributeValue::Double(v) => v.first().copied(),
        AttributeValue::Float(v) => v.first().map(|f| f64::from(*f)),
        AttributeValue::Int(v) => v.first().map(|i| f64::from(*i)),
        AttributeValue::Short(v) => v.first().map(|s| f64::from(*s)),
        _ => None,
    }
}

fn infer_cell_size(metadata: &NetCdfMetadata, var: &Variable) -> f64 {
    // Try to find a coordinate variable for the last X dimension
    if let Some(x_dim) = var.dimensions.last() {
        if let Some(x_var) = metadata
            .variables
            .iter()
            .find(|v| v.name == *x_dim && v.shape.len() == 1)
        {
            if x_var.shape[0] > 1 {
                // Rough estimate: assume uniform spacing from axis metadata
                if let Some(AttributeValue::Text(units)) = x_var.attributes.get("units") {
                    if units.contains("degree") {
                        // Assume global 0.25° or similar
                        return 360.0 / x_var.shape[0] as f64;
                    }
                }
            }
        }
    }
    1.0 // default
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    #[test]
    fn test_data_type_sizes() {
        assert_eq!(DataType::Byte.size_bytes(), 1);
        assert_eq!(DataType::Short.size_bytes(), 2);
        assert_eq!(DataType::Float.size_bytes(), 4);
        assert_eq!(DataType::Double.size_bytes(), 8);
    }

    #[test]
    fn test_get_attr_f64() {
        let mut attrs = HashMap::new();
        attrs.insert(
            "scale_factor".to_string(),
            AttributeValue::Double(vec![0.01]),
        );
        attrs.insert(
            "add_offset".to_string(),
            AttributeValue::Float(vec![273.15]),
        );

        assert_eq!(get_attr_f64(&attrs, "scale_factor"), Some(0.01));
        assert!((get_attr_f64(&attrs, "add_offset").unwrap() - 273.15).abs() < 0.01);
        assert_eq!(get_attr_f64(&attrs, "missing"), None);
    }

    fn put_name(buf: &mut Vec<u8>, s: &str) {
        buf.extend_from_slice(&(s.len() as u32).to_be_bytes());
        buf.extend_from_slice(s.as_bytes());
        let pad = (4 - (s.len() % 4)) % 4;
        buf.extend(std::iter::repeat_n(0u8, pad));
    }

    fn classic_two_var_netcdf() -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(b"CDF\x01");
        buf.extend_from_slice(&0u32.to_be_bytes());

        buf.extend_from_slice(&10u32.to_be_bytes());
        buf.extend_from_slice(&2u32.to_be_bytes());
        put_name(&mut buf, "y");
        buf.extend_from_slice(&2u32.to_be_bytes());
        put_name(&mut buf, "x");
        buf.extend_from_slice(&2u32.to_be_bytes());

        buf.extend_from_slice(&0u32.to_be_bytes());
        buf.extend_from_slice(&0u32.to_be_bytes());

        buf.extend_from_slice(&11u32.to_be_bytes());
        buf.extend_from_slice(&2u32.to_be_bytes());

        put_name(&mut buf, "temp");
        buf.extend_from_slice(&2u32.to_be_bytes());
        buf.extend_from_slice(&0u32.to_be_bytes());
        buf.extend_from_slice(&1u32.to_be_bytes());
        buf.extend_from_slice(&0u32.to_be_bytes());
        buf.extend_from_slice(&0u32.to_be_bytes());
        buf.extend_from_slice(&5u32.to_be_bytes());
        buf.extend_from_slice(&16u32.to_be_bytes());
        let temp_begin_at = buf.len();
        buf.extend_from_slice(&0u32.to_be_bytes());

        put_name(&mut buf, "precip");
        buf.extend_from_slice(&2u32.to_be_bytes());
        buf.extend_from_slice(&0u32.to_be_bytes());
        buf.extend_from_slice(&1u32.to_be_bytes());
        buf.extend_from_slice(&0u32.to_be_bytes());
        buf.extend_from_slice(&0u32.to_be_bytes());
        buf.extend_from_slice(&5u32.to_be_bytes());
        buf.extend_from_slice(&16u32.to_be_bytes());
        let precip_begin_at = buf.len();
        buf.extend_from_slice(&0u32.to_be_bytes());

        let temp_begin = buf.len() as u32;
        for v in [1.0f32, 2.0, 3.0, 4.0] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        let precip_begin = buf.len() as u32;
        for v in [10.0f32, 20.0, 30.0, 40.0] {
            buf.extend_from_slice(&v.to_be_bytes());
        }
        buf[temp_begin_at..temp_begin_at + 4].copy_from_slice(&temp_begin.to_be_bytes());
        buf[precip_begin_at..precip_begin_at + 4].copy_from_slice(&precip_begin.to_be_bytes());
        buf
    }

    #[test]
    fn test_read_second_variable_uses_begin_offset() {
        let bytes = classic_two_var_netcdf();
        let mut reader = BufReader::new(Cursor::new(bytes));
        let meta = read_netcdf_metadata(&mut reader).unwrap();
        assert_eq!(meta.variables.len(), 2);
        assert!(meta.variables[1].begin > meta.variables[0].begin);

        let precip = read_netcdf_variable(&mut reader, &meta, "precip").unwrap();
        assert_eq!(precip.width(), 2);
        assert_eq!(precip.height(), 2);
        assert!((precip.get(0, 0).unwrap() - 10.0).abs() < 1e-5);
        assert!((precip.get(0, 1).unwrap() - 20.0).abs() < 1e-5);
        assert!((precip.get(1, 0).unwrap() - 30.0).abs() < 1e-5);
        assert!((precip.get(1, 1).unwrap() - 40.0).abs() < 1e-5);
    }
}
