//! Polygons onto a grid, the inverse of `polygonize`. A cell takes a polygon's
//! value when its centre falls inside, which is what makes a set of boundaries
//! usable as the zone raster `zonal_stats` reads.

use crate::{Raster, RegionPolygon};

/// Burn `polygons` onto a `width` x `height` grid spanning `bbox`
/// (xmin, ymin, xmax, ymax), with y running north-down as rasters do.
///
/// A cell belongs to a polygon when its centre is inside, counting crossings
/// over every ring so holes cut out. Overlapping polygons resolve last-wins,
/// and a cell no polygon covers stays nodata.
///
/// Note the y axis: `polygonize` emits rings in cell units counting rows
/// downward, so its output needs flipping to north-up before it comes back
/// through here.
pub fn rasterize(
    polygons: &[RegionPolygon],
    width: usize,
    height: usize,
    bbox: (f64, f64, f64, f64),
    cell_size: f64,
    nodata: f64,
) -> Raster {
    let mut result = Raster::new(width, height, cell_size, nodata);
    if width == 0 || height == 0 {
        return result;
    }
    let (xmin, ymin, xmax, ymax) = bbox;
    let cell_w = (xmax - xmin) / width as f64;
    let cell_h = (ymax - ymin) / height as f64;
    if cell_w <= 0.0 || cell_h <= 0.0 {
        return result;
    }

    let mut crossings: Vec<f64> = Vec::new();
    for polygon in polygons {
        for row in 0..height {
            let y = ymax - (row as f64 + 0.5) * cell_h;
            crossings.clear();
            for ring in &polygon.rings {
                for edge in ring.windows(2) {
                    let ((x1, y1), (x2, y2)) = (edge[0], edge[1]);
                    if (y1 > y) != (y2 > y) {
                        crossings.push(x1 + (y - y1) / (y2 - y1) * (x2 - x1));
                    }
                }
            }
            if crossings.is_empty() {
                continue;
            }
            crossings.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for span in crossings.chunks(2) {
                if span.len() < 2 {
                    continue;
                }
                // the cells whose centre x falls in the span
                let lo = ((span[0] - xmin) / cell_w - 0.5).ceil();
                let hi = ((span[1] - xmin) / cell_w - 0.5).floor();
                if hi < 0.0 || lo > width as f64 - 1.0 {
                    continue;
                }
                let first = lo.max(0.0) as usize;
                let last = hi.min(width as f64 - 1.0) as usize;
                for col in first..=last {
                    result.set(row, col, polygon.value);
                }
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polygonize;

    const NODATA: f64 = -9999.0;

    fn square(value: f64, x0: f64, y0: f64, x1: f64, y1: f64) -> RegionPolygon {
        RegionPolygon {
            value,
            rings: vec![vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]],
        }
    }

    #[test]
    fn a_square_covering_the_grid_fills_every_cell() {
        let out = rasterize(
            &[square(3.0, 0.0, 0.0, 4.0, 4.0)],
            4,
            4,
            (0.0, 0.0, 4.0, 4.0),
            1.0,
            NODATA,
        );

        assert!(out.data().iter().all(|&v| v == 3.0));
    }

    #[test]
    fn a_half_square_fills_the_cells_whose_centres_it_covers() {
        // covers x 0..2 of a 4-wide grid, so columns 0 and 1
        let out = rasterize(
            &[square(1.0, 0.0, 0.0, 2.0, 4.0)],
            4,
            4,
            (0.0, 0.0, 4.0, 4.0),
            1.0,
            NODATA,
        );

        assert_eq!(out.get(0, 0).unwrap(), 1.0);
        assert_eq!(out.get(0, 1).unwrap(), 1.0);
        assert!(out.is_nodata(out.get(0, 2).unwrap()));
        assert!(out.is_nodata(out.get(0, 3).unwrap()));
    }

    #[test]
    fn a_hole_is_cut_out_of_its_polygon() {
        let donut = RegionPolygon {
            value: 5.0,
            rings: vec![
                vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)],
                // the middle 2x2, wound the other way
                vec![(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0), (1.0, 1.0)],
            ],
        };
        let out = rasterize(&[donut], 4, 4, (0.0, 0.0, 4.0, 4.0), 1.0, NODATA);

        assert_eq!(out.get(0, 0).unwrap(), 5.0);
        assert!(out.is_nodata(out.get(1, 1).unwrap()));
        assert!(out.is_nodata(out.get(2, 2).unwrap()));
        assert_eq!(out.get(3, 3).unwrap(), 5.0);
    }

    #[test]
    fn overlapping_polygons_resolve_last_wins() {
        let out = rasterize(
            &[
                square(1.0, 0.0, 0.0, 4.0, 4.0),
                square(2.0, 0.0, 0.0, 4.0, 4.0),
            ],
            2,
            2,
            (0.0, 0.0, 4.0, 4.0),
            1.0,
            NODATA,
        );

        assert!(out.data().iter().all(|&v| v == 2.0));
    }

    #[test]
    fn polygonize_then_rasterize_returns_the_grid_it_started_from() {
        let data = vec![
            1.0, 1.0, 2.0, 2.0, //
            1.0, 1.0, 2.0, 2.0, //
            3.0, 3.0, 3.0, 3.0, //
            3.0, 3.0, 3.0, 3.0,
        ];
        let src = Raster::from_vec(4, 4, data, 1.0, NODATA).unwrap();
        // polygonize counts rows downward, rasterize reads north-up, so the
        // rings flip before they come back
        let polygons: Vec<RegionPolygon> = polygonize(&src)
            .into_iter()
            .map(|p| RegionPolygon {
                value: p.value,
                rings: p
                    .rings
                    .into_iter()
                    .map(|ring| ring.into_iter().map(|(x, y)| (x, 4.0 - y)).collect())
                    .collect(),
            })
            .collect();

        let back = rasterize(&polygons, 4, 4, (0.0, 0.0, 4.0, 4.0), 1.0, NODATA);

        assert_eq!(back.data(), src.data());
    }
}
