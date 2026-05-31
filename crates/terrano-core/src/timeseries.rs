//! Time-series analysis for raster stacks.
//!
//! Provides Earth Observation workflows including:
//! - Temporal compositing (mean, median, max NDVI)
//! - Change detection (differencing, thresholding)
//! - Trend analysis (linear regression per pixel)
//! - Anomaly detection (z-score deviation from historical mean)
//! - Phenology extraction (season start/end from vegetation indices)

use crate::{Error, Raster};

/// A temporal stack of rasters with associated timestamps.
#[derive(Debug, Clone)]
pub struct RasterStack {
    /// Rasters ordered by time (oldest first).
    pub layers: Vec<Raster>,
    /// Timestamps as days since epoch (or any monotonic sequence).
    pub timestamps: Vec<f64>,
}

/// Composite method for temporal reduction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CompositeMethod {
    Mean,
    Median,
    Max,
    Min,
    StdDev,
}

/// Result of per-pixel linear trend analysis.
#[derive(Debug, Clone)]
pub struct TrendResult {
    /// Slope of the linear fit (units per time unit).
    pub slope: Raster,
    /// Intercept of the linear fit.
    pub intercept: Raster,
    /// R² coefficient of determination.
    pub r_squared: Raster,
}

/// Change detection result.
#[derive(Debug, Clone)]
pub struct ChangeResult {
    /// Magnitude of change (absolute difference).
    pub magnitude: Raster,
    /// Binary mask: 1.0 where change exceeds threshold.
    pub change_mask: Raster,
}

/// Phenology metrics extracted from a vegetation index time series.
#[derive(Debug, Clone)]
pub struct PhenologyMetrics {
    /// Day (timestamp value) of season start (green-up).
    pub start_of_season: Raster,
    /// Day (timestamp value) of peak greenness.
    pub peak: Raster,
    /// Day (timestamp value) of season end (senescence).
    pub end_of_season: Raster,
    /// Maximum value during the growing season.
    pub max_value: Raster,
}

impl RasterStack {
    /// Create a new raster stack. All rasters must have same dimensions.
    pub fn new(layers: Vec<Raster>, timestamps: Vec<f64>) -> Result<Self, Error> {
        if layers.is_empty() {
            return Err(Error::InvalidInput(
                "RasterStack requires at least one layer".into(),
            ));
        }
        if layers.len() != timestamps.len() {
            return Err(Error::DimensionMismatch {
                expected: layers.len(),
                got: timestamps.len(),
            });
        }
        let w = layers[0].width();
        let h = layers[0].height();
        for (i, r) in layers.iter().enumerate().skip(1) {
            if r.width() != w || r.height() != h {
                return Err(Error::InvalidInput(format!(
                    "Layer {i} dimensions {}x{} differ from first layer {w}x{h}",
                    r.width(),
                    r.height()
                )));
            }
        }
        Ok(Self { layers, timestamps })
    }

    /// Temporal composite: reduce the stack to a single raster.
    pub fn composite(&self, method: CompositeMethod) -> Raster {
        let w = self.layers[0].width();
        let h = self.layers[0].height();
        let nodata = self.layers[0].nodata;
        let n = w * h;
        let mut out = vec![nodata; n];

        for (i, cell) in out.iter_mut().enumerate() {
            let values: Vec<f64> = self
                .layers
                .iter()
                .map(|r| r.data()[i])
                .filter(|&v| v != nodata && !v.is_nan())
                .collect();

            if values.is_empty() {
                continue;
            }

            *cell = match method {
                CompositeMethod::Mean => values.iter().sum::<f64>() / values.len() as f64,
                CompositeMethod::Median => {
                    let mut sorted = values.clone();
                    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                    let mid = sorted.len() / 2;
                    if sorted.len() % 2 == 0 {
                        (sorted[mid - 1] + sorted[mid]) / 2.0
                    } else {
                        sorted[mid]
                    }
                }
                CompositeMethod::Max => values.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                CompositeMethod::Min => values.iter().cloned().fold(f64::INFINITY, f64::min),
                CompositeMethod::StdDev => {
                    let mean = values.iter().sum::<f64>() / values.len() as f64;
                    let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>()
                        / values.len() as f64;
                    variance.sqrt()
                }
            };
        }

        Raster::from_vec(w, h, out, self.layers[0].cell_size, nodata).unwrap()
    }

    /// Per-pixel linear trend over time.
    pub fn linear_trend(&self) -> TrendResult {
        let w = self.layers[0].width();
        let h = self.layers[0].height();
        let nodata = self.layers[0].nodata;
        let n = w * h;

        let mut slopes = vec![nodata; n];
        let mut intercepts = vec![nodata; n];
        let mut r_squareds = vec![nodata; n];

        for i in 0..n {
            let pairs: Vec<(f64, f64)> = self
                .timestamps
                .iter()
                .zip(self.layers.iter())
                .map(|(&t, r)| (t, r.data()[i]))
                .filter(|&(_, v)| v != nodata && !v.is_nan())
                .collect();

            if pairs.len() < 2 {
                continue;
            }

            let count = pairs.len() as f64;
            let sum_x: f64 = pairs.iter().map(|(x, _)| x).sum();
            let sum_y: f64 = pairs.iter().map(|(_, y)| y).sum();
            let sum_xy: f64 = pairs.iter().map(|(x, y)| x * y).sum();
            let sum_x2: f64 = pairs.iter().map(|(x, _)| x * x).sum();

            let denom = count * sum_x2 - sum_x * sum_x;
            if denom.abs() < 1e-15 {
                continue;
            }

            let slope = (count * sum_xy - sum_x * sum_y) / denom;
            let intercept = (sum_y - slope * sum_x) / count;

            let mean_y = sum_y / count;
            let ss_tot: f64 = pairs.iter().map(|(_, y)| (y - mean_y).powi(2)).sum();
            let ss_res: f64 = pairs
                .iter()
                .map(|(x, y)| (y - (slope * x + intercept)).powi(2))
                .sum();

            let r2 = if ss_tot > 0.0 {
                1.0 - ss_res / ss_tot
            } else {
                0.0
            };

            slopes[i] = slope;
            intercepts[i] = intercept;
            r_squareds[i] = r2;
        }

        let cs = self.layers[0].cell_size;
        TrendResult {
            slope: Raster::from_vec(w, h, slopes, cs, nodata).unwrap(),
            intercept: Raster::from_vec(w, h, intercepts, cs, nodata).unwrap(),
            r_squared: Raster::from_vec(w, h, r_squareds, cs, nodata).unwrap(),
        }
    }

    /// Change detection between first and last layer in the stack.
    pub fn change_detection(&self, threshold: f64) -> ChangeResult {
        let first = &self.layers[0];
        let last = self.layers.last().unwrap();
        let w = first.width();
        let h = first.height();
        let nodata = first.nodata;
        let n = w * h;

        let mut magnitude = vec![nodata; n];
        let mut mask = vec![0.0; n];

        for i in 0..n {
            let v0 = first.data()[i];
            let v1 = last.data()[i];
            if v0 == nodata || v1 == nodata || v0.is_nan() || v1.is_nan() {
                mask[i] = nodata;
                continue;
            }
            let diff = (v1 - v0).abs();
            magnitude[i] = diff;
            mask[i] = if diff >= threshold { 1.0 } else { 0.0 };
        }

        ChangeResult {
            magnitude: Raster::from_vec(w, h, magnitude, first.cell_size, nodata).unwrap(),
            change_mask: Raster::from_vec(w, h, mask, first.cell_size, nodata).unwrap(),
        }
    }

    /// Anomaly detection: z-score of the last observation vs historical mean/stddev.
    pub fn anomaly_zscore(&self) -> Raster {
        let w = self.layers[0].width();
        let h = self.layers[0].height();
        let nodata = self.layers[0].nodata;
        let n = w * h;
        let mut out = vec![nodata; n];

        if self.layers.len() < 3 {
            return Raster::from_vec(w, h, out, self.layers[0].cell_size, nodata).unwrap();
        }

        // Use all layers except the last as historical baseline
        let baseline = &self.layers[..self.layers.len() - 1];
        let current = self.layers.last().unwrap();

        for (i, cell) in out.iter_mut().enumerate() {
            let values: Vec<f64> = baseline
                .iter()
                .map(|r| r.data()[i])
                .filter(|&v| v != nodata && !v.is_nan())
                .collect();

            let cv = current.data()[i];
            if values.len() < 2 || cv == nodata || cv.is_nan() {
                continue;
            }

            let mean = values.iter().sum::<f64>() / values.len() as f64;
            let std = (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>()
                / values.len() as f64)
                .sqrt();

            if std > 1e-10 {
                *cell = (cv - mean) / std;
            } else {
                *cell = 0.0;
            }
        }

        Raster::from_vec(w, h, out, self.layers[0].cell_size, nodata).unwrap()
    }

    /// Extract phenology metrics from a vegetation index time series.
    ///
    /// Uses a threshold-based approach: season starts when value exceeds
    /// `threshold_fraction` of the amplitude, and ends when it drops below.
    pub fn phenology(&self, threshold_fraction: f64) -> PhenologyMetrics {
        let w = self.layers[0].width();
        let h = self.layers[0].height();
        let nodata = self.layers[0].nodata;
        let n = w * h;

        let mut sos = vec![nodata; n];
        let mut peak = vec![nodata; n];
        let mut eos = vec![nodata; n];
        let mut max_val = vec![nodata; n];

        for i in 0..n {
            let series: Vec<(f64, f64)> = self
                .timestamps
                .iter()
                .zip(self.layers.iter())
                .map(|(&t, r)| (t, r.data()[i]))
                .filter(|&(_, v)| v != nodata && !v.is_nan())
                .collect();

            if series.len() < 3 {
                continue;
            }

            let min_v = series.iter().map(|(_, v)| *v).fold(f64::INFINITY, f64::min);
            let max_v = series
                .iter()
                .map(|(_, v)| *v)
                .fold(f64::NEG_INFINITY, f64::max);
            let amplitude = max_v - min_v;

            if amplitude < 1e-10 {
                continue;
            }

            let threshold = min_v + amplitude * threshold_fraction;

            // Find peak
            let peak_idx = series
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(idx, _)| idx)
                .unwrap();

            peak[i] = series[peak_idx].0;
            max_val[i] = series[peak_idx].1;

            // SOS: first crossing above threshold before peak
            for entry in series.iter().take(peak_idx) {
                if entry.1 >= threshold {
                    sos[i] = entry.0;
                    break;
                }
            }

            // EOS: first crossing below threshold after peak
            for entry in series.iter().skip(peak_idx + 1) {
                if entry.1 < threshold {
                    eos[i] = entry.0;
                    break;
                }
            }
        }

        let cs = self.layers[0].cell_size;
        PhenologyMetrics {
            start_of_season: Raster::from_vec(w, h, sos, cs, nodata).unwrap(),
            peak: Raster::from_vec(w, h, peak, cs, nodata).unwrap(),
            end_of_season: Raster::from_vec(w, h, eos, cs, nodata).unwrap(),
            max_value: Raster::from_vec(w, h, max_val, cs, nodata).unwrap(),
        }
    }

    /// Compute spectral index (e.g., NDVI) from two bands per time step.
    /// Formula: (band_a - band_b) / (band_a + band_b)
    pub fn normalized_difference(
        band_a: &RasterStack,
        band_b: &RasterStack,
    ) -> Result<RasterStack, Error> {
        if band_a.layers.len() != band_b.layers.len() {
            return Err(Error::DimensionMismatch {
                expected: band_a.layers.len(),
                got: band_b.layers.len(),
            });
        }

        let mut result_layers = Vec::with_capacity(band_a.layers.len());
        let nodata = band_a.layers[0].nodata;

        for (a, b) in band_a.layers.iter().zip(band_b.layers.iter()) {
            let w = a.width();
            let h = a.height();
            let mut out = vec![nodata; w * h];

            for (cell, (&va, &vb)) in out.iter_mut().zip(a.data().iter().zip(b.data().iter())) {
                if va == nodata || vb == nodata || va.is_nan() || vb.is_nan() {
                    continue;
                }
                let sum = va + vb;
                if sum.abs() > 1e-10 {
                    *cell = (va - vb) / sum;
                }
            }

            result_layers.push(Raster::from_vec(w, h, out, a.cell_size, nodata)?);
        }

        RasterStack::new(result_layers, band_a.timestamps.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_raster(data: Vec<f64>) -> Raster {
        let n = data.len();
        Raster::from_vec(n, 1, data, 1.0, -9999.0).unwrap()
    }

    #[test]
    fn composite_mean() {
        let stack = RasterStack::new(
            vec![
                make_raster(vec![1.0, 2.0, 3.0]),
                make_raster(vec![3.0, 4.0, 5.0]),
                make_raster(vec![5.0, 6.0, 7.0]),
            ],
            vec![1.0, 2.0, 3.0],
        )
        .unwrap();

        let result = stack.composite(CompositeMethod::Mean);
        assert!((result.data()[0] - 3.0).abs() < 1e-10);
        assert!((result.data()[1] - 4.0).abs() < 1e-10);
        assert!((result.data()[2] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn composite_median() {
        let stack = RasterStack::new(
            vec![
                make_raster(vec![1.0, 10.0]),
                make_raster(vec![5.0, 20.0]),
                make_raster(vec![3.0, 30.0]),
            ],
            vec![1.0, 2.0, 3.0],
        )
        .unwrap();

        let result = stack.composite(CompositeMethod::Median);
        assert!((result.data()[0] - 3.0).abs() < 1e-10);
        assert!((result.data()[1] - 20.0).abs() < 1e-10);
    }

    #[test]
    fn composite_max_min() {
        let stack = RasterStack::new(
            vec![make_raster(vec![1.0, 5.0]), make_raster(vec![3.0, 2.0])],
            vec![1.0, 2.0],
        )
        .unwrap();

        let max_r = stack.composite(CompositeMethod::Max);
        assert!((max_r.data()[0] - 3.0).abs() < 1e-10);
        assert!((max_r.data()[1] - 5.0).abs() < 1e-10);

        let min_r = stack.composite(CompositeMethod::Min);
        assert!((min_r.data()[0] - 1.0).abs() < 1e-10);
        assert!((min_r.data()[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn linear_trend_positive() {
        let stack = RasterStack::new(
            vec![
                make_raster(vec![0.0]),
                make_raster(vec![1.0]),
                make_raster(vec![2.0]),
                make_raster(vec![3.0]),
            ],
            vec![0.0, 1.0, 2.0, 3.0],
        )
        .unwrap();

        let trend = stack.linear_trend();
        assert!((trend.slope.data()[0] - 1.0).abs() < 1e-10);
        assert!((trend.intercept.data()[0]).abs() < 1e-10);
        assert!((trend.r_squared.data()[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn linear_trend_flat() {
        let stack = RasterStack::new(
            vec![
                make_raster(vec![5.0]),
                make_raster(vec![5.0]),
                make_raster(vec![5.0]),
            ],
            vec![0.0, 1.0, 2.0],
        )
        .unwrap();

        let trend = stack.linear_trend();
        assert!((trend.slope.data()[0]).abs() < 1e-10);
    }

    #[test]
    fn change_detection_above_threshold() {
        let stack = RasterStack::new(
            vec![
                make_raster(vec![1.0, 5.0, 10.0]),
                make_raster(vec![4.0, 5.5, 10.2]),
            ],
            vec![0.0, 1.0],
        )
        .unwrap();

        let change = stack.change_detection(2.0);
        assert!((change.magnitude.data()[0] - 3.0).abs() < 1e-10);
        assert!((change.change_mask.data()[0] - 1.0).abs() < 1e-10); // exceeds threshold
        assert!((change.change_mask.data()[1]).abs() < 1e-10); // below threshold
    }

    #[test]
    fn anomaly_zscore_high() {
        // Historical: 1,2,3,4 -> mean=2.5, std≈1.12
        // Current: 10 -> z ≈ (10-2.5)/1.12 ≈ 6.7
        let stack = RasterStack::new(
            vec![
                make_raster(vec![1.0]),
                make_raster(vec![2.0]),
                make_raster(vec![3.0]),
                make_raster(vec![4.0]),
                make_raster(vec![10.0]), // anomaly
            ],
            vec![0.0, 1.0, 2.0, 3.0, 4.0],
        )
        .unwrap();

        let zscore = stack.anomaly_zscore();
        assert!(zscore.data()[0] > 5.0); // clearly anomalous
    }

    #[test]
    fn phenology_basic() {
        // Simulates a vegetation growth cycle
        let stack = RasterStack::new(
            vec![
                make_raster(vec![0.1]), // winter
                make_raster(vec![0.3]), // early spring
                make_raster(vec![0.7]), // late spring
                make_raster(vec![0.9]), // summer peak
                make_raster(vec![0.6]), // early fall
                make_raster(vec![0.2]), // late fall
            ],
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        )
        .unwrap();

        let pheno = stack.phenology(0.3);
        // Peak should be at timestamp 4.0
        assert!((pheno.peak.data()[0] - 4.0).abs() < 1e-10);
        // SOS should be at timestamp 2.0 (first > 0.1 + 0.8*0.3 = 0.34)
        assert!((pheno.start_of_season.data()[0] - 3.0).abs() < 1e-10);
        assert!((pheno.max_value.data()[0] - 0.9).abs() < 1e-10);
    }

    #[test]
    fn normalized_difference_ndvi() {
        let nir = RasterStack::new(vec![make_raster(vec![0.8, 0.6])], vec![1.0]).unwrap();

        let red = RasterStack::new(vec![make_raster(vec![0.2, 0.4])], vec![1.0]).unwrap();

        let ndvi = RasterStack::normalized_difference(&nir, &red).unwrap();
        // (0.8-0.2)/(0.8+0.2) = 0.6
        assert!((ndvi.layers[0].data()[0] - 0.6).abs() < 1e-10);
        // (0.6-0.4)/(0.6+0.4) = 0.2
        assert!((ndvi.layers[0].data()[1] - 0.2).abs() < 1e-10);
    }

    #[test]
    fn stack_dimension_mismatch() {
        let result = RasterStack::new(
            vec![make_raster(vec![1.0])],
            vec![1.0, 2.0], // wrong number of timestamps
        );
        assert!(result.is_err());
    }

    #[test]
    fn nodata_handling() {
        let stack = RasterStack::new(
            vec![
                make_raster(vec![-9999.0, 2.0]),
                make_raster(vec![3.0, -9999.0]),
            ],
            vec![1.0, 2.0],
        )
        .unwrap();

        let result = stack.composite(CompositeMethod::Mean);
        assert!((result.data()[0] - 3.0).abs() < 1e-10); // only one valid value
        assert!((result.data()[1] - 2.0).abs() < 1e-10);
    }
}
