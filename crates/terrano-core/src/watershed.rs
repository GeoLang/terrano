use crate::Raster;

/// Delineate watersheds from a D8 flow direction raster.
///
/// Each cell is assigned a watershed ID based on which pour point it drains to.
/// Pour points are cells at the raster edge or cells with no outflow.
///
/// # Arguments
/// * `flow_dir` — D8 flow direction raster (from `flow_direction()`)
///
/// # Returns
/// A raster where each cell value is a watershed ID (starting from 1).
pub fn watershed(flow_dir: &Raster) -> Raster {
    let w = flow_dir.width();
    let h = flow_dir.height();
    let mut labels = vec![0u32; w * h];
    let mut next_label = 1u32;

    let dr: [i32; 8] = [0, 1, 1, 1, 0, -1, -1, -1];
    let dc: [i32; 8] = [1, 1, 0, -1, -1, -1, 0, 1];
    let codes: [f64; 8] = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0];

    // For each cell, trace downstream until we reach a labeled cell or a pour point
    for row in 0..h {
        for col in 0..w {
            if labels[row * w + col] != 0 {
                continue;
            }

            let dir = flow_dir.get(row, col).unwrap();
            if flow_dir.is_nodata(dir) {
                continue;
            }

            // Trace path downstream
            let mut path = Vec::new();
            let mut cr = row;
            let mut cc = col;
            let mut found_label = 0u32;

            loop {
                let idx = cr * w + cc;
                if labels[idx] != 0 {
                    found_label = labels[idx];
                    break;
                }
                path.push(idx);

                let dir_val = flow_dir.get(cr, cc).unwrap();
                if flow_dir.is_nodata(dir_val) || dir_val == 0.0 {
                    // Pour point (pit or flat)
                    break;
                }

                let mut moved = false;
                for d in 0..8 {
                    if (dir_val - codes[d]).abs() < 0.5 {
                        let nr = cr as i32 + dr[d];
                        let nc = cc as i32 + dc[d];
                        if nr >= 0 && nr < h as i32 && nc >= 0 && nc < w as i32 {
                            cr = nr as usize;
                            cc = nc as usize;
                            moved = true;
                        }
                        break;
                    }
                }

                if !moved {
                    // Edge pour point
                    break;
                }

                // Detect cycles
                if path.len() > w * h {
                    break;
                }
            }

            // Assign label to all cells in the path
            let label = if found_label != 0 {
                found_label
            } else {
                let l = next_label;
                next_label += 1;
                l
            };

            for &idx in &path {
                labels[idx] = label;
            }
        }
    }

    // Convert to raster
    let mut result = Raster::new(w, h, flow_dir.cell_size, flow_dir.nodata);
    for row in 0..h {
        for col in 0..w {
            let label = labels[row * w + col];
            if label > 0 {
                result.set(row, col, label as f64);
            }
        }
    }
    result
}

/// Strahler stream ordering from a flow accumulation raster.
///
/// # Arguments
/// * `flow_accum` — flow accumulation raster
/// * `flow_dir` — D8 flow direction raster
/// * `threshold` — minimum accumulation to be considered a stream cell
///
/// # Returns
/// Raster with Strahler stream order values (1, 2, 3, ...).
pub fn stream_order(flow_accum: &Raster, flow_dir: &Raster, threshold: f64) -> Raster {
    let w = flow_dir.width();
    let h = flow_dir.height();
    let mut order = Raster::new(w, h, flow_dir.cell_size, flow_dir.nodata);

    let dr: [i32; 8] = [0, 1, 1, 1, 0, -1, -1, -1];
    let dc: [i32; 8] = [1, 1, 0, -1, -1, -1, 0, 1];
    let codes: [f64; 8] = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0];

    // Initialize: headwater cells (stream cells with no upstream stream tributaries)
    // get order 1
    let mut stream_cells: Vec<(usize, usize)> = Vec::new();
    for row in 0..h {
        for col in 0..w {
            let accum = flow_accum.get(row, col).unwrap();
            if flow_dir.is_nodata(accum) || accum < threshold {
                continue;
            }
            stream_cells.push((row, col));
        }
    }

    // Find headwater cells (stream cells with no upstream stream cell flowing into them)
    let mut in_degree = vec![0u32; w * h];
    for &(row, col) in &stream_cells {
        let dir_val = flow_dir.get(row, col).unwrap();
        if flow_dir.is_nodata(dir_val) || dir_val == 0.0 {
            continue;
        }
        for d in 0..8 {
            if (dir_val - codes[d]).abs() < 0.5 {
                let nr = row as i32 + dr[d];
                let nc = col as i32 + dc[d];
                if nr >= 0 && nr < h as i32 && nc >= 0 && nc < w as i32 {
                    let target = nr as usize * w + nc as usize;
                    let target_accum = flow_accum.get(nr as usize, nc as usize).unwrap();
                    if target_accum >= threshold {
                        in_degree[target] += 1;
                    }
                }
                break;
            }
        }
    }

    // BFS/topological sort from headwaters
    let mut queue: std::collections::VecDeque<(usize, usize)> = std::collections::VecDeque::new();
    for &(row, col) in &stream_cells {
        if in_degree[row * w + col] == 0 {
            order.set(row, col, 1.0);
            queue.push_back((row, col));
        }
    }

    while let Some((row, col)) = queue.pop_front() {
        let current_order = order.get(row, col).unwrap() as u32;
        let dir_val = flow_dir.get(row, col).unwrap();
        if flow_dir.is_nodata(dir_val) || dir_val == 0.0 {
            continue;
        }

        for d in 0..8 {
            if (dir_val - codes[d]).abs() < 0.5 {
                let nr = row as i32 + dr[d];
                let nc = col as i32 + dc[d];
                if nr >= 0 && nr < h as i32 && nc >= 0 && nc < w as i32 {
                    let nr = nr as usize;
                    let nc = nc as usize;
                    let target_accum = flow_accum.get(nr, nc).unwrap();
                    if target_accum < threshold {
                        break;
                    }

                    let existing = order.get(nr, nc).unwrap_or(0.0) as u32;
                    let new_order = if existing == current_order {
                        // Two tributaries of same order merge → order + 1
                        current_order + 1
                    } else {
                        current_order.max(existing)
                    };
                    order.set(nr, nc, new_order as f64);

                    in_degree[nr * w + nc] -= 1;
                    if in_degree[nr * w + nc] == 0 {
                        queue.push_back((nr, nc));
                    }
                }
                break;
            }
        }
    }

    order
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hydrology::flow_direction;

    #[test]
    fn test_watershed_simple() {
        // 5x5 DEM with two valleys draining to different edges
        #[rustfmt::skip]
        let data = vec![
            10.0,  8.0,  9.0,  8.0, 10.0,
             8.0,  5.0,  9.0,  5.0,  8.0,
             6.0,  3.0,  9.0,  3.0,  6.0,
             4.0,  2.0,  9.0,  2.0,  4.0,
             2.0,  1.0,  9.0,  1.0,  2.0,
        ];
        let dem = Raster::from_vec(5, 5, data, 1.0, -9999.0).unwrap();
        let fd = flow_direction(&dem);
        let ws = watershed(&fd);

        // Left side should have different label than right side
        let left = ws.get(2, 1).unwrap();
        let right = ws.get(2, 3).unwrap();
        assert!(left > 0.0);
        assert!(right > 0.0);
        assert_ne!(
            left, right,
            "different valleys should have different watershed IDs"
        );
    }

    #[test]
    fn test_watershed_single_basin() {
        // Cone-shaped DEM: everything flows to center bottom
        #[rustfmt::skip]
        let data = vec![
            9.0, 8.0, 7.0, 8.0, 9.0,
            8.0, 6.0, 5.0, 6.0, 8.0,
            7.0, 4.0, 3.0, 4.0, 7.0,
            6.0, 3.0, 2.0, 3.0, 6.0,
            5.0, 2.0, 1.0, 2.0, 5.0,
        ];
        let dem = Raster::from_vec(5, 5, data, 1.0, -9999.0).unwrap();
        let fd = flow_direction(&dem);
        let ws = watershed(&fd);

        // All cells should belong to the same watershed (single outlet at center-bottom)
        let first = ws.get(0, 0).unwrap();
        assert!(first > 0.0);
        for row in 0..5 {
            for col in 0..5 {
                let val = ws.get(row, col).unwrap();
                if val > 0.0 {
                    assert_eq!(
                        val, first,
                        "cone DEM should be single watershed at ({},{})",
                        row, col
                    );
                }
            }
        }
    }
}
