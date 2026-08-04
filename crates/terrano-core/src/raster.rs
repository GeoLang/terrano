use crate::Error;

/// A 2D raster grid.
#[derive(Debug, Clone)]
pub struct Raster {
    /// Row-major pixel values.
    data: Vec<f64>,
    /// Number of columns.
    width: usize,
    /// Number of rows.
    height: usize,
    /// Cell size in coordinate units.
    pub cell_size: f64,
    /// No-data sentinel value.
    pub nodata: f64,
}

impl Raster {
    pub fn new(width: usize, height: usize, cell_size: f64, nodata: f64) -> Self {
        Self {
            data: vec![nodata; width * height],
            width,
            height,
            cell_size,
            nodata,
        }
    }

    pub fn from_vec(
        width: usize,
        height: usize,
        data: Vec<f64>,
        cell_size: f64,
        nodata: f64,
    ) -> Result<Self, Error> {
        if data.len() != width * height {
            return Err(Error::DimensionMismatch {
                expected: width * height,
                got: data.len(),
            });
        }
        Ok(Self {
            data,
            width,
            height,
            cell_size,
            nodata,
        })
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn height(&self) -> usize {
        self.height
    }

    pub fn get(&self, row: usize, col: usize) -> Option<f64> {
        if row < self.height && col < self.width {
            Some(self.data[row * self.width + col])
        } else {
            None
        }
    }

    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        if row < self.height && col < self.width {
            self.data[row * self.width + col] = value;
        }
    }

    pub fn is_nodata(&self, value: f64) -> bool {
        value == self.nodata || value.is_nan()
    }

    pub fn data(&self) -> &[f64] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Consume the raster, returning its row-major values.
    pub fn into_data(self) -> Vec<f64> {
        self.data
    }
}
