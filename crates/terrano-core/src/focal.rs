//! Moving-window statistics. Every output cell summarizes its neighbourhood,
//! which is how a noisy grid gets smoothed, a classified one gets its majority
//! filter, and local relief gets measured.

use crate::Raster;

/// Which summary a focal pass reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FocalStat {
    Min,
    Max,
    Mean,
    Sum,
    Std,
    Median,
    Majority,
    /// max - min, the local relief
    Range,
}

/// Window shape at a given radius.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Neighborhood {
    Square,
    Circle,
}

/// Summarize each cell's neighbourhood.
///
/// The window clips at the raster edge and skips nodata neighbours, so a
/// border cell reports over the part of the window that exists. A nodata cell
/// stays nodata rather than being filled from its neighbours, which would grow
/// the data area on every pass.
pub fn focal_stats(raster: &Raster, radius: usize, shape: Neighborhood, stat: FocalStat) -> Raster {
    let (w, h) = (raster.width(), raster.height());
    let mut result = Raster::new(w, h, raster.cell_size, raster.nodata);
    let offsets = window_offsets(radius, shape);
    let mut sample: Vec<f64> = Vec::with_capacity(offsets.len());

    for row in 0..h {
        for col in 0..w {
            let center = raster.get(row, col).unwrap();
            if raster.is_nodata(center) {
                continue;
            }
            sample.clear();
            for &(dr, dc) in &offsets {
                let r = row as isize + dr;
                let c = col as isize + dc;
                if r < 0 || c < 0 || r >= h as isize || c >= w as isize {
                    continue;
                }
                let v = raster.get(r as usize, c as usize).unwrap();
                if !raster.is_nodata(v) {
                    sample.push(v);
                }
            }
            if !sample.is_empty() {
                result.set(row, col, summarize(&mut sample, stat));
            }
        }
    }
    result
}

fn window_offsets(radius: usize, shape: Neighborhood) -> Vec<(isize, isize)> {
    let r = radius as isize;
    let mut offsets = Vec::new();
    for dr in -r..=r {
        for dc in -r..=r {
            if shape == Neighborhood::Circle && dr * dr + dc * dc > r * r {
                continue;
            }
            offsets.push((dr, dc));
        }
    }
    offsets
}

/// Sorts `sample` in place for the order statistics.
fn summarize(sample: &mut [f64], stat: FocalStat) -> f64 {
    match stat {
        FocalStat::Min => sample.iter().copied().fold(f64::INFINITY, f64::min),
        FocalStat::Max => sample.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        FocalStat::Sum => sample.iter().sum(),
        FocalStat::Mean => sample.iter().sum::<f64>() / sample.len() as f64,
        FocalStat::Std => {
            let mean = sample.iter().sum::<f64>() / sample.len() as f64;
            let var =
                sample.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / sample.len() as f64;
            var.sqrt()
        }
        FocalStat::Range => {
            let min = sample.iter().copied().fold(f64::INFINITY, f64::min);
            let max = sample.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            max - min
        }
        FocalStat::Median => {
            sample.sort_by(|a, b| a.partial_cmp(b).unwrap());
            median_of_sorted(sample)
        }
        FocalStat::Majority => {
            sample.sort_by(|a, b| a.partial_cmp(b).unwrap());
            majority_of_sorted(sample)
        }
    }
}

pub(crate) fn median_of_sorted(sorted: &[f64]) -> f64 {
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

/// The most frequent value, the smallest of them on a tie so the result does
/// not depend on iteration order.
fn majority_of_sorted(sorted: &[f64]) -> f64 {
    let mut best = sorted[0];
    let mut best_count = 0;
    let mut current = sorted[0];
    let mut count = 0;
    for &v in sorted {
        if v == current {
            count += 1;
        } else {
            current = v;
            count = 1;
        }
        if count > best_count {
            best_count = count;
            best = current;
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    const NODATA: f64 = -9999.0;

    fn raster(width: usize, height: usize, data: Vec<f64>) -> Raster {
        Raster::from_vec(width, height, data, 1.0, NODATA).unwrap()
    }

    #[test]
    fn mean_over_a_full_window_averages_the_neighbourhood() {
        let src = raster(3, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let out = focal_stats(&src, 1, Neighborhood::Square, FocalStat::Mean);

        // the centre sees all nine cells
        assert_eq!(out.get(1, 1).unwrap(), 5.0);
        // the top-left corner sees only the 2x2 it clips to
        assert_eq!(out.get(0, 0).unwrap(), (1.0 + 2.0 + 4.0 + 5.0) / 4.0);
    }

    #[test]
    fn min_max_and_range_read_the_window_extremes() {
        let src = raster(3, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);

        let min = focal_stats(&src, 1, Neighborhood::Square, FocalStat::Min);
        let max = focal_stats(&src, 1, Neighborhood::Square, FocalStat::Max);
        let range = focal_stats(&src, 1, Neighborhood::Square, FocalStat::Range);

        assert_eq!(min.get(1, 1).unwrap(), 1.0);
        assert_eq!(max.get(1, 1).unwrap(), 9.0);
        assert_eq!(range.get(1, 1).unwrap(), 8.0);
    }

    #[test]
    fn median_and_majority_disagree_where_the_mode_is_not_central() {
        // eight 7s around a 1: median 7, majority 7, mean pulled down
        let src = raster(3, 3, vec![7.0, 7.0, 7.0, 7.0, 1.0, 7.0, 7.0, 7.0, 7.0]);

        let median = focal_stats(&src, 1, Neighborhood::Square, FocalStat::Median);
        let majority = focal_stats(&src, 1, Neighborhood::Square, FocalStat::Majority);
        let mean = focal_stats(&src, 1, Neighborhood::Square, FocalStat::Mean);

        assert_eq!(median.get(1, 1).unwrap(), 7.0);
        assert_eq!(majority.get(1, 1).unwrap(), 7.0);
        assert!(mean.get(1, 1).unwrap() < 7.0);
    }

    #[test]
    fn a_circle_window_drops_the_diagonal_corners() {
        let src = raster(3, 3, vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0]);
        let square = focal_stats(&src, 1, Neighborhood::Square, FocalStat::Max);
        let circle = focal_stats(&src, 1, Neighborhood::Circle, FocalStat::Max);

        // the 100 sits diagonally from the centre, inside the square window only
        assert_eq!(square.get(1, 1).unwrap(), 100.0);
        assert_eq!(circle.get(1, 1).unwrap(), 0.0);
    }

    #[test]
    fn nodata_is_left_alone_and_kept_out_of_its_neighbours_windows() {
        let src = raster(3, 1, vec![4.0, NODATA, 6.0]);
        let out = focal_stats(&src, 1, Neighborhood::Square, FocalStat::Mean);

        assert!(out.is_nodata(out.get(0, 1).unwrap()));
        // the left cell averages itself alone, the gap contributing nothing
        assert_eq!(out.get(0, 0).unwrap(), 4.0);
    }

    #[test]
    fn a_zero_radius_window_returns_each_cell_unchanged() {
        let src = raster(2, 2, vec![1.0, 2.0, 3.0, 4.0]);
        let out = focal_stats(&src, 0, Neighborhood::Square, FocalStat::Mean);

        assert_eq!(out.data(), src.data());
    }

    #[test]
    fn std_of_a_flat_window_is_zero() {
        let src = raster(3, 3, vec![5.0; 9]);
        let out = focal_stats(&src, 1, Neighborhood::Square, FocalStat::Std);

        assert_eq!(out.get(1, 1).unwrap(), 0.0);
    }
}
