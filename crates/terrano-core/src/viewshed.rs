use crate::Raster;

/// Cell value where the observer has line of sight to the ground.
pub const VIEWSHED_VISIBLE: f64 = 1.0;

/// Cell value where terrain blocks the observer's line of sight.
pub const VIEWSHED_HIDDEN: f64 = 0.0;

/// Line-of-sight viewshed from one observer cell.
///
/// Casts a ray from the observer's eye to every cell inside `radius` and hides
/// the cell when the terrain the ray crosses stands at a steeper elevation
/// angle than the target itself. Cells past the radius and cells with no
/// elevation stay nodata, as does the whole raster when the observer cell has
/// no elevation or lies outside the grid.
///
/// Each ray is sampled once per row or column it crosses, whichever it spans
/// more of, so the work is the cell count times the ray length.
///
/// # Arguments
/// * `dem` — elevation raster
/// * `observer_row`, `observer_col` — the cell the observer stands on
/// * `observer_height` — eye height above that cell's ground, in elevation units
/// * `radius` — how far the observer looks, in the raster's coordinate units.
///   `f64::INFINITY` covers the whole raster.
pub fn viewshed(
    dem: &Raster,
    observer_row: usize,
    observer_col: usize,
    observer_height: f64,
    radius: f64,
) -> Raster {
    let mut result = Raster::new(dem.width(), dem.height(), dem.cell_size, dem.nodata);
    let Some(ground) = dem.get(observer_row, observer_col) else {
        return result;
    };
    if dem.is_nodata(ground) {
        return result;
    }
    let eye = ground + observer_height;
    let radius_cells = radius / dem.cell_size;
    result.set(observer_row, observer_col, VIEWSHED_VISIBLE);

    for row in 0..dem.height() {
        for col in 0..dem.width() {
            if row == observer_row && col == observer_col {
                continue;
            }
            let target = dem.get(row, col).unwrap();
            if dem.is_nodata(target) {
                continue;
            }
            let row_span = row as f64 - observer_row as f64;
            let col_span = col as f64 - observer_col as f64;
            if row_span.hypot(col_span) > radius_cells {
                continue;
            }
            let clear = has_line_of_sight(
                dem,
                observer_row,
                observer_col,
                eye,
                row_span,
                col_span,
                target,
            );
            let value = if clear {
                VIEWSHED_VISIBLE
            } else {
                VIEWSHED_HIDDEN
            };
            result.set(row, col, value);
        }
    }
    result
}

/// Whether the ray from the eye to the cell `row_span, col_span` away clears
/// the terrain between them.
///
/// Angles are compared as rise over run in cell units, so `cell_size` cancels.
/// A sample exactly on the line of sight grazes it rather than blocking it,
/// which is what keeps a flat surface visible from an eye at ground level. A
/// nodata sample obstructs nothing: there is no elevation there to block with.
fn has_line_of_sight(
    dem: &Raster,
    observer_row: usize,
    observer_col: usize,
    eye: f64,
    row_span: f64,
    col_span: f64,
    target: f64,
) -> bool {
    let length = row_span.hypot(col_span);
    let target_angle = (target - eye) / length;
    let steps = row_span.abs().max(col_span.abs()) as usize;
    for step in 1..steps {
        let along = step as f64 / steps as f64;
        let elevation = sample_bilinear(
            dem,
            observer_row as f64 + row_span * along,
            observer_col as f64 + col_span * along,
        );
        let Some(elevation) = elevation else {
            continue;
        };
        if (elevation - eye) / (length * along) > target_angle {
            return false;
        }
    }
    true
}

/// Bilinear elevation at a fractional cell position inside the grid. `None`
/// when any of the four cells around it is nodata.
fn sample_bilinear(dem: &Raster, row: f64, col: f64) -> Option<f64> {
    let row0 = row.floor() as usize;
    let col0 = col.floor() as usize;
    // a position on the far edge floors onto the last cell, so clamp the second
    // corner rather than reading past it
    let row1 = (row0 + 1).min(dem.height().saturating_sub(1));
    let col1 = (col0 + 1).min(dem.width().saturating_sub(1));
    let down = row - row0 as f64;
    let right = col - col0 as f64;

    let mut elevation = 0.0;
    for (r, c, weight) in [
        (row0, col0, (1.0 - down) * (1.0 - right)),
        (row0, col1, (1.0 - down) * right),
        (row1, col0, down * (1.0 - right)),
        (row1, col1, down * right),
    ] {
        let value = dem.get(r, c)?;
        if dem.is_nodata(value) {
            return None;
        }
        elevation += value * weight;
    }
    Some(elevation)
}

#[cfg(test)]
mod tests {
    use super::*;

    const NODATA: f64 = -9999.0;

    fn flat(width: usize, height: usize) -> Raster {
        Raster::from_vec(width, height, vec![0.0; width * height], 1.0, NODATA).unwrap()
    }

    /// A wall across every row, so nothing beyond it can be seen around it.
    fn walled(width: usize, height: usize, wall_col: usize, wall_height: f64) -> Raster {
        let mut dem = flat(width, height);
        for row in 0..height {
            dem.set(row, wall_col, wall_height);
        }
        dem
    }

    #[test]
    fn flat_ground_is_all_visible() {
        let dem = flat(9, 9);
        let vs = viewshed(&dem, 4, 4, 2.0, f64::INFINITY);
        for row in 0..9 {
            for col in 0..9 {
                assert_eq!(vs.get(row, col).unwrap(), VIEWSHED_VISIBLE, "({row},{col})");
            }
        }
    }

    #[test]
    fn a_wall_hides_the_ground_behind_it() {
        let (width, height) = (11usize, 5usize);
        let dem = walled(width, height, 5, 50.0);
        let vs = viewshed(&dem, 2, 0, 2.0, f64::INFINITY);

        for col in 1..5 {
            assert_eq!(
                vs.get(2, col).unwrap(),
                VIEWSHED_VISIBLE,
                "col {col} stands in front of the wall"
            );
        }
        assert_eq!(
            vs.get(2, 5).unwrap(),
            VIEWSHED_VISIBLE,
            "the wall top is above the eye"
        );
        for row in 0..height {
            for col in 6..width {
                assert_eq!(
                    vs.get(row, col).unwrap(),
                    VIEWSHED_HIDDEN,
                    "({row},{col}) lies behind the wall"
                );
            }
        }
    }

    /// The shadow the wall casts, which is the part a radial sweep of the
    /// terrain profile cannot express: an eye high enough sees over the wall,
    /// but only past the ground the wall top still occludes.
    #[test]
    fn an_eye_above_the_wall_sees_over_it() {
        let dem = walled(11, 5, 5, 50.0);
        let vs = viewshed(&dem, 2, 0, 200.0, f64::INFINITY);
        // eye 200, wall top 50 five cells out: the shadow reaches 6.67 cells
        assert_eq!(
            vs.get(2, 6).unwrap(),
            VIEWSHED_HIDDEN,
            "the shadow just behind the wall"
        );
        for col in 7..11 {
            assert_eq!(
                vs.get(2, col).unwrap(),
                VIEWSHED_VISIBLE,
                "col {col} lies past the shadow"
            );
        }
    }

    #[test]
    fn a_peak_sees_further_than_a_pit() {
        // a cone rising to the middle column, seen from its top and from its foot
        let (width, height) = (11usize, 5usize);
        let mut dem = flat(width, height);
        for row in 0..height {
            for col in 0..width {
                dem.set(row, col, 20.0 - (col as f64 - 5.0).abs() * 4.0);
            }
        }
        let count_visible =
            |vs: &Raster| vs.data().iter().filter(|&&v| v == VIEWSHED_VISIBLE).count();
        let peak = viewshed(&dem, 2, 5, 2.0, f64::INFINITY);
        let foot = viewshed(&dem, 2, 0, 2.0, f64::INFINITY);
        assert!(
            count_visible(&peak) > count_visible(&foot),
            "peak {} vs foot {}",
            count_visible(&peak),
            count_visible(&foot)
        );
    }

    #[test]
    fn cells_past_the_radius_have_no_verdict() {
        let dem = Raster::from_vec(9, 9, vec![0.0; 81], 10.0, NODATA).unwrap();
        let vs = viewshed(&dem, 4, 4, 2.0, 25.0);
        // cells are 10 units wide, so the radius reaches 2.5 of them
        assert_eq!(vs.get(4, 6).unwrap(), VIEWSHED_VISIBLE);
        assert!(vs.is_nodata(vs.get(4, 7).unwrap()));
    }

    #[test]
    fn an_observer_with_no_ground_under_it_sees_nothing() {
        let mut dem = flat(5, 5);
        dem.set(2, 2, NODATA);
        let vs = viewshed(&dem, 2, 2, 2.0, f64::INFINITY);
        assert!(vs.data().iter().all(|&v| vs.is_nodata(v)));

        let outside = viewshed(&flat(5, 5), 9, 9, 2.0, f64::INFINITY);
        assert!(outside.data().iter().all(|&v| outside.is_nodata(v)));
    }

    #[test]
    fn nodata_terrain_blocks_nothing_and_gets_no_verdict() {
        let mut dem = flat(7, 3);
        dem.set(1, 3, NODATA);
        let vs = viewshed(&dem, 1, 0, 2.0, f64::INFINITY);
        assert!(vs.is_nodata(vs.get(1, 3).unwrap()));
        assert_eq!(vs.get(1, 6).unwrap(), VIEWSHED_VISIBLE);
    }
}
