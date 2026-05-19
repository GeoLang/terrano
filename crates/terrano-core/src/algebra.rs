use crate::{Error, Raster};

/// Unary raster operation (applied cell-by-cell).
pub enum UnaryOp {
    Add(f64),
    Multiply(f64),
    Sqrt,
    Abs,
    Log,
}

/// Binary raster operation (applied cell-by-cell between two rasters).
pub enum BinaryOp {
    Add,
    Subtract,
    Multiply,
    Divide,
    Min,
    Max,
}

impl Raster {
    /// Apply a unary operation to every cell, returning a new raster.
    pub fn apply_unary(&self, op: &UnaryOp) -> Raster {
        let mut result = self.clone();
        for val in result.data_mut() {
            if !self.is_nodata(*val) {
                *val = match op {
                    UnaryOp::Add(c) => *val + c,
                    UnaryOp::Multiply(c) => *val * c,
                    UnaryOp::Sqrt => val.sqrt(),
                    UnaryOp::Abs => val.abs(),
                    UnaryOp::Log => val.ln(),
                };
            }
        }
        result
    }

    /// Combine two rasters cell-by-cell.
    pub fn apply_binary(&self, other: &Raster, op: &BinaryOp) -> Result<Raster, Error> {
        if self.width() != other.width() || self.height() != other.height() {
            return Err(Error::IncompatibleRasters);
        }
        let mut result = self.clone();
        for (i, val) in result.data_mut().iter_mut().enumerate() {
            let a = self.data()[i];
            let b = other.data()[i];
            if self.is_nodata(a) || other.is_nodata(b) {
                *val = self.nodata;
            } else {
                *val = match op {
                    BinaryOp::Add => a + b,
                    BinaryOp::Subtract => a - b,
                    BinaryOp::Multiply => a * b,
                    BinaryOp::Divide => {
                        if b == 0.0 {
                            self.nodata
                        } else {
                            a / b
                        }
                    }
                    BinaryOp::Min => a.min(b),
                    BinaryOp::Max => a.max(b),
                };
            }
        }
        Ok(result)
    }
}

/// Reclassify raster values based on breakpoints.
/// Each entry in `classes` is (min_inclusive, max_exclusive, new_value).
pub fn reclassify(raster: &Raster, classes: &[(f64, f64, f64)]) -> Raster {
    let mut result = Raster::new(
        raster.width(),
        raster.height(),
        raster.cell_size,
        raster.nodata,
    );
    for row in 0..raster.height() {
        for col in 0..raster.width() {
            let val = raster.get(row, col).unwrap();
            if raster.is_nodata(val) {
                continue;
            }
            let classified = classes
                .iter()
                .find(|(min, max, _)| val >= *min && val < *max)
                .map(|(_, _, new)| *new)
                .unwrap_or(raster.nodata);
            result.set(row, col, classified);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unary_add() {
        let r = Raster::from_vec(2, 2, vec![1.0, 2.0, 3.0, 4.0], 1.0, -9999.0).unwrap();
        let result = r.apply_unary(&UnaryOp::Add(10.0));
        assert_eq!(result.data(), &[11.0, 12.0, 13.0, 14.0]);
    }

    #[test]
    fn test_binary_add() {
        let a = Raster::from_vec(2, 2, vec![1.0, 2.0, 3.0, 4.0], 1.0, -9999.0).unwrap();
        let b = Raster::from_vec(2, 2, vec![10.0, 20.0, 30.0, 40.0], 1.0, -9999.0).unwrap();
        let result = a.apply_binary(&b, &BinaryOp::Add).unwrap();
        assert_eq!(result.data(), &[11.0, 22.0, 33.0, 44.0]);
    }

    #[test]
    fn test_reclassify() {
        let r = Raster::from_vec(2, 2, vec![1.0, 5.0, 15.0, 25.0], 1.0, -9999.0).unwrap();
        let classes = vec![(0.0, 10.0, 1.0), (10.0, 20.0, 2.0), (20.0, 30.0, 3.0)];
        let result = reclassify(&r, &classes);
        assert_eq!(result.data(), &[1.0, 1.0, 2.0, 3.0]);
    }
}
