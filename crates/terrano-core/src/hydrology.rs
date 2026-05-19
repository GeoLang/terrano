use crate::Raster;

/// D8 flow direction encoding.
/// Encodes direction as power of 2: 1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE
pub fn flow_direction(dem: &Raster) -> Raster {
    let mut result = Raster::new(dem.width(), dem.height(), dem.cell_size, dem.nodata);
    let diag = dem.cell_size * std::f64::consts::SQRT_2;

    // D8 neighbor offsets: E, SE, S, SW, W, NW, N, NE
    let dr: [i32; 8] = [0, 1, 1, 1, 0, -1, -1, -1];
    let dc: [i32; 8] = [1, 1, 0, -1, -1, -1, 0, 1];
    let dist: [f64; 8] = [
        dem.cell_size,
        diag,
        dem.cell_size,
        diag,
        dem.cell_size,
        diag,
        dem.cell_size,
        diag,
    ];
    let codes: [f64; 8] = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0];

    for row in 0..dem.height() {
        for col in 0..dem.width() {
            let center = dem.get(row, col).unwrap();
            if dem.is_nodata(center) {
                continue;
            }

            let mut max_drop = 0.0;
            let mut best_dir = 0.0;

            for d in 0..8 {
                let nr = row as i32 + dr[d];
                let nc = col as i32 + dc[d];
                if nr < 0 || nr >= dem.height() as i32 || nc < 0 || nc >= dem.width() as i32 {
                    continue;
                }
                let neighbor = dem.get(nr as usize, nc as usize).unwrap();
                if dem.is_nodata(neighbor) {
                    continue;
                }
                let drop = (center - neighbor) / dist[d];
                if drop > max_drop {
                    max_drop = drop;
                    best_dir = codes[d];
                }
            }
            result.set(row, col, best_dir);
        }
    }
    result
}

/// Compute flow accumulation from a D8 flow direction raster.
/// Each cell's value represents the number of upstream cells that flow into it.
pub fn flow_accumulation(flow_dir: &Raster) -> Raster {
    let w = flow_dir.width();
    let h = flow_dir.height();
    let mut accum = vec![0u32; w * h];
    let mut visited = vec![false; w * h];

    // Trace flow for each cell
    for row in 0..h {
        for col in 0..w {
            if visited[row * w + col] {
                continue;
            }
            trace_flow(flow_dir, &mut accum, &mut visited, row, col, w, h);
        }
    }

    let mut result = Raster::new(w, h, flow_dir.cell_size, flow_dir.nodata);
    for row in 0..h {
        for col in 0..w {
            let val = flow_dir.get(row, col).unwrap();
            if flow_dir.is_nodata(val) {
                continue;
            }
            result.set(row, col, accum[row * w + col] as f64);
        }
    }
    result
}

fn trace_flow(
    flow_dir: &Raster,
    accum: &mut [u32],
    visited: &mut [bool],
    start_row: usize,
    start_col: usize,
    w: usize,
    h: usize,
) {
    let dr: [i32; 8] = [0, 1, 1, 1, 0, -1, -1, -1];
    let dc: [i32; 8] = [1, 1, 0, -1, -1, -1, 0, 1];
    let codes: [f64; 8] = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0];

    let mut path = Vec::new();
    let mut row = start_row;
    let mut col = start_col;

    // Follow flow downstream, collecting path
    loop {
        let idx = row * w + col;
        if visited[idx] {
            // Already counted, just add our path contribution downstream
            break;
        }
        path.push((row, col));
        visited[idx] = true;

        let dir = flow_dir.get(row, col).unwrap();
        if flow_dir.is_nodata(dir) || dir == 0.0 {
            break; // Pit or edge
        }

        // Find which neighbor this direction points to
        let mut found = false;
        for d in 0..8 {
            if (dir - codes[d]).abs() < 0.5 {
                let nr = row as i32 + dr[d];
                let nc = col as i32 + dc[d];
                if nr >= 0 && nr < h as i32 && nc >= 0 && nc < w as i32 {
                    row = nr as usize;
                    col = nc as usize;
                    found = true;
                }
                break;
            }
        }
        if !found {
            break;
        }
    }

    // Add accumulation from this path
    // Each cell in the path contributes 1 to all downstream cells
    for &(r, c) in &path {
        // Trace downstream from this cell and add 1
        let mut cr = r;
        let mut cc = c;
        // The cell itself gets counted by upstream cells flowing into it
        loop {
            let dir = flow_dir.get(cr, cc).unwrap();
            if flow_dir.is_nodata(dir) || dir == 0.0 {
                break;
            }
            let mut moved = false;
            for d in 0..8 {
                if (dir - codes[d]).abs() < 0.5 {
                    let nr = cr as i32 + dr[d];
                    let nc = cc as i32 + dc[d];
                    if nr >= 0 && nr < h as i32 && nc >= 0 && nc < w as i32 {
                        cr = nr as usize;
                        cc = nc as usize;
                        accum[cr * w + cc] += 1;
                        moved = true;
                    }
                    break;
                }
            }
            if !moved {
                break;
            }
            // Prevent infinite loops
            if cr == r && cc == c {
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flow_direction_simple_slope() {
        // 3x3 DEM sloping from NW to SE
        #[rustfmt::skip]
        let data = vec![
            9.0, 8.0, 7.0,
            6.0, 5.0, 4.0,
            3.0, 2.0, 1.0,
        ];
        let dem = Raster::from_vec(3, 3, data, 1.0, -9999.0).unwrap();
        let fd = flow_direction(&dem);
        // Center cell (5.0): S neighbor (2.0) has drop 3/1=3.0, SE (1.0) has drop 4/√2=2.83
        // Steepest gradient is S (code 4)
        let center_dir = fd.get(1, 1).unwrap();
        assert_eq!(center_dir, 4.0); // S (steepest slope per unit distance)
    }

    #[test]
    fn test_flow_direction_east_slope() {
        // 3x3 DEM sloping east
        #[rustfmt::skip]
        let data = vec![
            3.0, 2.0, 1.0,
            3.0, 2.0, 1.0,
            3.0, 2.0, 1.0,
        ];
        let dem = Raster::from_vec(3, 3, data, 1.0, -9999.0).unwrap();
        let fd = flow_direction(&dem);
        let center_dir = fd.get(1, 1).unwrap();
        assert_eq!(center_dir, 1.0); // E
    }
}
