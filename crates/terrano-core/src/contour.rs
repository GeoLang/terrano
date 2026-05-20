use crate::Raster;

/// A contour line segment (two endpoints).
#[derive(Debug, Clone, PartialEq)]
pub struct ContourSegment {
    pub x1: f64,
    pub y1: f64,
    pub x2: f64,
    pub y2: f64,
}

/// A contour line (connected polyline at a single elevation).
#[derive(Debug, Clone)]
pub struct ContourLine {
    pub level: f64,
    pub vertices: Vec<(f64, f64)>,
}

/// Generate contour lines from a DEM raster using marching squares.
///
/// # Arguments
/// * `dem` — input elevation raster
/// * `interval` — elevation interval between contours
/// * `base` — base elevation to start contours from (often 0.0)
///
/// # Returns
/// Vector of contour lines, each with a level and vertex list.
pub fn contours(dem: &Raster, interval: f64, base: f64) -> Vec<ContourLine> {
    if interval <= 0.0 || dem.width() < 2 || dem.height() < 2 {
        return Vec::new();
    }

    // Find elevation range
    let mut min_z = f64::MAX;
    let mut max_z = f64::MIN;
    for row in 0..dem.height() {
        for col in 0..dem.width() {
            let v = dem.get(row, col).unwrap();
            if !dem.is_nodata(v) {
                min_z = min_z.min(v);
                max_z = max_z.max(v);
            }
        }
    }

    if min_z > max_z {
        return Vec::new();
    }

    // Generate contour levels
    let first_level = ((min_z - base) / interval).ceil() * interval + base;
    let mut levels = Vec::new();
    let mut level = first_level;
    while level <= max_z {
        levels.push(level);
        level += interval;
    }

    // Generate segments for each level using marching squares
    let mut result = Vec::new();
    for &level in &levels {
        let segments = march_squares(dem, level);
        let lines = connect_segments(segments, level);
        result.extend(lines);
    }

    result
}

/// Marching squares: generate line segments for a single contour level.
fn march_squares(dem: &Raster, level: f64) -> Vec<ContourSegment> {
    let mut segments = Vec::new();
    let cs = dem.cell_size;

    for row in 0..dem.height() - 1 {
        for col in 0..dem.width() - 1 {
            // Get four corner values (TL, TR, BR, BL)
            let tl = dem.get(row, col).unwrap();
            let tr = dem.get(row, col + 1).unwrap();
            let br = dem.get(row + 1, col + 1).unwrap();
            let bl = dem.get(row + 1, col).unwrap();

            if dem.is_nodata(tl) || dem.is_nodata(tr) || dem.is_nodata(br) || dem.is_nodata(bl) {
                continue;
            }

            // Cell origin (top-left corner in coordinate space)
            let x0 = col as f64 * cs;
            let y0 = (dem.height() - 1 - row) as f64 * cs; // Y increases upward

            // Classify corners: 1 if above or at level, 0 if below
            let case = ((tl >= level) as u8)
                | (((tr >= level) as u8) << 1)
                | (((br >= level) as u8) << 2)
                | (((bl >= level) as u8) << 3);

            if case == 0 || case == 15 {
                continue; // Fully above or below
            }

            // Interpolation positions on edges
            let top = lerp_edge(tl, tr, level) * cs + x0;
            let right = y0 - lerp_edge(tr, br, level) * cs;
            let bottom = lerp_edge(bl, br, level) * cs + x0;
            let left = y0 - lerp_edge(tl, bl, level) * cs;

            // Edge midpoints
            let top_pt = (top, y0);
            let right_pt = (x0 + cs, right);
            let bottom_pt = (bottom, y0 - cs);
            let left_pt = (x0, left);

            match case {
                1 => segments.push(seg(left_pt, top_pt)),
                2 => segments.push(seg(top_pt, right_pt)),
                3 => segments.push(seg(left_pt, right_pt)),
                4 => segments.push(seg(right_pt, bottom_pt)),
                5 => {
                    // Saddle: two segments
                    let center = (tl + tr + br + bl) / 4.0;
                    if center >= level {
                        segments.push(seg(left_pt, top_pt));
                        segments.push(seg(right_pt, bottom_pt));
                    } else {
                        segments.push(seg(left_pt, bottom_pt));
                        segments.push(seg(top_pt, right_pt));
                    }
                }
                6 => segments.push(seg(top_pt, bottom_pt)),
                7 => segments.push(seg(left_pt, bottom_pt)),
                8 => segments.push(seg(bottom_pt, left_pt)),
                9 => segments.push(seg(bottom_pt, top_pt)),
                10 => {
                    // Saddle: two segments
                    let center = (tl + tr + br + bl) / 4.0;
                    if center >= level {
                        segments.push(seg(top_pt, right_pt));
                        segments.push(seg(bottom_pt, left_pt));
                    } else {
                        segments.push(seg(left_pt, top_pt));
                        segments.push(seg(right_pt, bottom_pt));
                    }
                }
                11 => segments.push(seg(bottom_pt, right_pt)),
                12 => segments.push(seg(right_pt, left_pt)),
                13 => segments.push(seg(right_pt, top_pt)),
                14 => segments.push(seg(top_pt, left_pt)),
                _ => {}
            }
        }
    }

    segments
}

fn seg(a: (f64, f64), b: (f64, f64)) -> ContourSegment {
    ContourSegment {
        x1: a.0,
        y1: a.1,
        x2: b.0,
        y2: b.1,
    }
}

/// Linear interpolation along an edge to find where the contour crosses.
fn lerp_edge(v1: f64, v2: f64, level: f64) -> f64 {
    if (v2 - v1).abs() < 1e-12 {
        0.5
    } else {
        (level - v1) / (v2 - v1)
    }
}

/// Connect individual segments into polylines.
fn connect_segments(segments: Vec<ContourSegment>, level: f64) -> Vec<ContourLine> {
    if segments.is_empty() {
        return Vec::new();
    }

    let eps = 1e-8;
    let mut used = vec![false; segments.len()];
    let mut lines = Vec::new();

    for start_idx in 0..segments.len() {
        if used[start_idx] {
            continue;
        }
        used[start_idx] = true;

        let mut vertices = vec![
            (segments[start_idx].x1, segments[start_idx].y1),
            (segments[start_idx].x2, segments[start_idx].y2),
        ];

        // Grow forward
        loop {
            let tail = *vertices.last().unwrap();
            let mut found = false;
            for i in 0..segments.len() {
                if used[i] {
                    continue;
                }
                if close(tail, (segments[i].x1, segments[i].y1), eps) {
                    vertices.push((segments[i].x2, segments[i].y2));
                    used[i] = true;
                    found = true;
                    break;
                } else if close(tail, (segments[i].x2, segments[i].y2), eps) {
                    vertices.push((segments[i].x1, segments[i].y1));
                    used[i] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                break;
            }
        }

        // Grow backward
        loop {
            let head = vertices[0];
            let mut found = false;
            for i in 0..segments.len() {
                if used[i] {
                    continue;
                }
                if close(head, (segments[i].x2, segments[i].y2), eps) {
                    vertices.insert(0, (segments[i].x1, segments[i].y1));
                    used[i] = true;
                    found = true;
                    break;
                } else if close(head, (segments[i].x1, segments[i].y1), eps) {
                    vertices.insert(0, (segments[i].x2, segments[i].y2));
                    used[i] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                break;
            }
        }

        lines.push(ContourLine { level, vertices });
    }

    lines
}

fn close(a: (f64, f64), b: (f64, f64), eps: f64) -> bool {
    (a.0 - b.0).abs() < eps && (a.1 - b.1).abs() < eps
}

/// Fill sinks (depressions) in a DEM using a priority-flood algorithm.
///
/// This is a prerequisite for hydrological analysis — removes local minima
/// so that every cell can drain to the edge of the raster.
///
/// Based on Wang & Liu (2006) "An efficient method for identifying and filling
/// surface depressions in digital elevation models".
pub fn fill_sinks(dem: &Raster) -> Raster {
    use std::cmp::Ordering;
    use std::collections::BinaryHeap;

    let w = dem.width();
    let h = dem.height();
    let mut output = dem.clone();
    let mut in_queue = vec![false; w * h];

    // Min-heap (reverse ordering for BinaryHeap which is max-heap)
    #[derive(PartialEq)]
    struct Cell {
        elevation: f64,
        row: usize,
        col: usize,
    }

    impl Eq for Cell {}

    impl PartialOrd for Cell {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for Cell {
        fn cmp(&self, other: &Self) -> Ordering {
            // Reverse for min-heap behavior
            other
                .elevation
                .partial_cmp(&self.elevation)
                .unwrap_or(Ordering::Equal)
        }
    }

    let mut heap = BinaryHeap::new();

    // Initialize: add all edge cells to the priority queue
    for row in 0..h {
        for col in 0..w {
            let is_edge = row == 0 || row == h - 1 || col == 0 || col == w - 1;
            if is_edge {
                let val = dem.get(row, col).unwrap();
                if !dem.is_nodata(val) {
                    heap.push(Cell {
                        elevation: val,
                        row,
                        col,
                    });
                    in_queue[row * w + col] = true;
                }
            }
        }
    }

    // 8-connected neighbors
    let dr: [i32; 8] = [0, 1, 1, 1, 0, -1, -1, -1];
    let dc: [i32; 8] = [1, 1, 0, -1, -1, -1, 0, 1];

    // Process cells in elevation order
    while let Some(cell) = heap.pop() {
        for d in 0..8 {
            let nr = cell.row as i32 + dr[d];
            let nc = cell.col as i32 + dc[d];

            if nr < 0 || nr >= h as i32 || nc < 0 || nc >= w as i32 {
                continue;
            }

            let nr = nr as usize;
            let nc = nc as usize;
            let idx = nr * w + nc;

            if in_queue[idx] {
                continue;
            }

            in_queue[idx] = true;
            let neighbor_val = dem.get(nr, nc).unwrap();

            if dem.is_nodata(neighbor_val) {
                continue;
            }

            // If neighbor is lower than current cell, raise it (fill the sink)
            let filled_val = neighbor_val.max(cell.elevation);
            output.set(nr, nc, filled_val);

            heap.push(Cell {
                elevation: filled_val,
                row: nr,
                col: nc,
            });
        }
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contours_simple_slope() {
        // 5x5 DEM with linear slope from 0 to 4
        let mut data = Vec::new();
        for row in 0..5 {
            for _col in 0..5 {
                data.push(row as f64);
            }
        }
        let dem = Raster::from_vec(5, 5, data, 1.0, -9999.0).unwrap();
        let result = contours(&dem, 1.0, 0.0);

        // Should have contours at 1, 2, 3 (0 and 4 are exact boundaries)
        let levels: Vec<f64> = result.iter().map(|c| c.level).collect();
        assert!(levels.contains(&1.0));
        assert!(levels.contains(&2.0));
        assert!(levels.contains(&3.0));
    }

    #[test]
    fn test_contours_have_vertices() {
        // 4x4 ramp
        let data: Vec<f64> = (0..16).map(|i| (i / 4) as f64 * 10.0).collect();
        let dem = Raster::from_vec(4, 4, data, 1.0, -9999.0).unwrap();
        let result = contours(&dem, 5.0, 0.0);

        for line in &result {
            assert!(
                line.vertices.len() >= 2,
                "contour at {} has {} vertices",
                line.level,
                line.vertices.len()
            );
        }
    }

    #[test]
    fn test_fill_sinks_flat() {
        // Already flat — no sinks to fill
        let data = vec![5.0; 25];
        let dem = Raster::from_vec(5, 5, data, 1.0, -9999.0).unwrap();
        let filled = fill_sinks(&dem);
        for row in 0..5 {
            for col in 0..5 {
                assert_eq!(filled.get(row, col).unwrap(), 5.0);
            }
        }
    }

    #[test]
    fn test_fill_sinks_single_pit() {
        // 5x5 grid: edges slope from 10 down to 3 at bottom-center,
        // interior at 5, with center pit at 1.
        // The outlet is at bottom edge (value 3), so interior fills to 5 (min surrounding).
        #[rustfmt::skip]
        let data = vec![
            10.0, 10.0, 10.0, 10.0, 10.0,
            10.0,  5.0,  5.0,  5.0, 10.0,
            10.0,  5.0,  1.0,  5.0, 10.0,
            10.0,  5.0,  5.0,  5.0, 10.0,
            10.0, 10.0,  3.0, 10.0, 10.0,
        ];
        let dem = Raster::from_vec(5, 5, data, 1.0, -9999.0).unwrap();
        let filled = fill_sinks(&dem);

        // The center pit should be raised to 5.0 (level of surrounding cells)
        // since water can now drain out through the 3.0 outlet at bottom
        let center = filled.get(2, 2).unwrap();
        assert_eq!(center, 5.0, "pit should be filled to 5.0, got {center}");

        // Edge cells should remain unchanged
        assert_eq!(filled.get(0, 0).unwrap(), 10.0);
        // Outlet should remain unchanged
        assert_eq!(filled.get(4, 2).unwrap(), 3.0);
    }

    #[test]
    fn test_fill_sinks_preserves_drainage() {
        // Slope from top to bottom — no sinks, should be unchanged
        let mut data = Vec::new();
        for row in 0..5 {
            for _col in 0..5 {
                data.push((4 - row) as f64); // decreasing from top to bottom
            }
        }
        let dem = Raster::from_vec(5, 5, data.clone(), 1.0, -9999.0).unwrap();
        let filled = fill_sinks(&dem);
        for (i, &expected) in data.iter().enumerate() {
            let row = i / 5;
            let col = i % 5;
            assert_eq!(filled.get(row, col).unwrap(), expected);
        }
    }
}
