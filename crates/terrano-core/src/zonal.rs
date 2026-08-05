//! Statistics of one raster grouped by the zones of another. The zone raster
//! carries a label per cell, whatever a reclass or a rasterized polygon set
//! put there, and every distinct label becomes one row.

use crate::focal::median_of_sorted;
use crate::{Error, Raster};
use std::collections::HashMap;

/// Summary of the value cells falling in one zone.
#[derive(Debug, Clone, PartialEq)]
pub struct ZoneStats {
    pub zone: f64,
    pub count: usize,
    pub min: f64,
    pub max: f64,
    pub mean: f64,
    pub sum: f64,
    pub std: f64,
    pub median: f64,
}

/// Summarize `values` per zone, one row per distinct zone label, ordered by
/// label. A cell counts only where both rasters carry data.
pub fn zonal_stats(values: &Raster, zones: &Raster) -> Result<Vec<ZoneStats>, Error> {
    if values.width() != zones.width() || values.height() != zones.height() {
        return Err(Error::DimensionMismatch {
            expected: values.width() * values.height(),
            got: zones.width() * zones.height(),
        });
    }

    let mut index: HashMap<u64, usize> = HashMap::new();
    let mut labels: Vec<f64> = Vec::new();
    let mut samples: Vec<Vec<f64>> = Vec::new();

    for row in 0..values.height() {
        for col in 0..values.width() {
            let zone = zones.get(row, col).unwrap();
            let value = values.get(row, col).unwrap();
            if zones.is_nodata(zone) || values.is_nodata(value) {
                continue;
            }
            // f64 has no Hash, and the bit pattern is exact for the finite
            // labels a zone raster holds
            let slot = *index.entry(zone.to_bits()).or_insert_with(|| {
                labels.push(zone);
                samples.push(Vec::new());
                labels.len() - 1
            });
            samples[slot].push(value);
        }
    }

    let mut out: Vec<ZoneStats> = labels
        .into_iter()
        .zip(samples)
        .map(|(zone, mut sample)| {
            sample.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let count = sample.len();
            let sum: f64 = sample.iter().sum();
            let mean = sum / count as f64;
            let var = sample.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / count as f64;
            ZoneStats {
                zone,
                count,
                min: sample[0],
                max: sample[count - 1],
                mean,
                sum,
                std: var.sqrt(),
                median: median_of_sorted(&sample),
            }
        })
        .collect();
    out.sort_by(|a, b| a.zone.partial_cmp(&b.zone).unwrap());
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const NODATA: f64 = -9999.0;

    fn raster(width: usize, height: usize, data: Vec<f64>) -> Raster {
        Raster::from_vec(width, height, data, 1.0, NODATA).unwrap()
    }

    #[test]
    fn each_zone_gets_one_row_ordered_by_label() {
        let values = raster(4, 1, vec![10.0, 20.0, 30.0, 40.0]);
        let zones = raster(4, 1, vec![2.0, 1.0, 2.0, 1.0]);

        let stats = zonal_stats(&values, &zones).unwrap();

        assert_eq!(stats.len(), 2);
        assert_eq!(stats[0].zone, 1.0);
        assert_eq!(stats[1].zone, 2.0);
        assert_eq!(stats[0].mean, 30.0); // 20 and 40
        assert_eq!(stats[1].mean, 20.0); // 10 and 30
    }

    #[test]
    fn the_summary_covers_count_extremes_sum_and_spread() {
        let values = raster(4, 1, vec![1.0, 2.0, 3.0, 100.0]);
        let zones = raster(4, 1, vec![7.0; 4]);

        let stats = zonal_stats(&values, &zones).unwrap();

        let z = &stats[0];
        assert_eq!(z.count, 4);
        assert_eq!(z.min, 1.0);
        assert_eq!(z.max, 100.0);
        assert_eq!(z.sum, 106.0);
        assert_eq!(z.mean, 26.5);
        assert_eq!(z.median, 2.5);
        assert!(z.std > 0.0);
    }

    #[test]
    fn a_cell_counts_only_where_both_rasters_carry_data() {
        let values = raster(4, 1, vec![10.0, NODATA, 30.0, f64::NAN]);
        let zones = raster(4, 1, vec![1.0, 1.0, NODATA, 1.0]);

        let stats = zonal_stats(&values, &zones).unwrap();

        assert_eq!(stats.len(), 1);
        assert_eq!(stats[0].count, 1);
        assert_eq!(stats[0].sum, 10.0);
    }

    #[test]
    fn mismatched_grids_are_an_error() {
        let values = raster(4, 1, vec![1.0; 4]);
        let zones = raster(2, 1, vec![1.0; 2]);

        assert!(zonal_stats(&values, &zones).is_err());
    }
}
