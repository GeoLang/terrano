//! Raster to vector: connected runs of equal-valued cells become polygons
//! whose rings follow cell corners. Values compare exactly, so this suits a
//! classified raster (a reclass, a land-cover map) and not a continuous one,
//! where every cell would come back as its own square.

use crate::Raster;
use std::collections::HashMap;

/// One connected run of equal-valued cells.
#[derive(Debug, Clone)]
pub struct RegionPolygon {
    pub value: f64,
    /// Closed rings in coordinate units, exterior first then its holes.
    pub rings: Vec<Vec<(f64, f64)>>,
}

/// Cell-corner grid position, so ring vertices compare and hash exactly.
type Vertex = (usize, usize);

const UNSET: usize = usize::MAX;

/// Polygonize a raster with 4-connected regions of exactly equal value.
/// Nodata cells belong to no region and bound the ones around them.
pub fn polygonize(raster: &Raster) -> Vec<RegionPolygon> {
    let w = raster.width();
    let h = raster.height();
    if w == 0 || h == 0 {
        return Vec::new();
    }

    let (labels, values) = label_regions(raster);
    let mut edges: Vec<Vec<(Vertex, Vertex)>> = vec![Vec::new(); values.len()];
    for i in 0..w * h {
        let id = labels[i];
        if id == UNSET {
            continue;
        }
        let (r, c) = (i / w, i % w);
        let same = |nr: isize, nc: isize| {
            nr >= 0
                && nc >= 0
                && (nr as usize) < h
                && (nc as usize) < w
                && labels[nr as usize * w + nc as usize] == id
        };
        // clockwise in a y-down grid, so an outer ring encloses positive area
        if !same(r as isize - 1, c as isize) {
            edges[id].push(((c, r), (c + 1, r)));
        }
        if !same(r as isize, c as isize + 1) {
            edges[id].push(((c + 1, r), (c + 1, r + 1)));
        }
        if !same(r as isize + 1, c as isize) {
            edges[id].push(((c + 1, r + 1), (c, r + 1)));
        }
        if !same(r as isize, c as isize - 1) {
            edges[id].push(((c, r + 1), (c, r)));
        }
    }

    let cs = raster.cell_size;
    let mut out = Vec::new();
    for (id, value) in values.into_iter().enumerate() {
        for rings in assemble(chain_rings(std::mem::take(&mut edges[id]))) {
            out.push(RegionPolygon {
                value,
                rings: rings
                    .into_iter()
                    .map(|ring| {
                        ring.into_iter()
                            .map(|(x, y)| (x as f64 * cs, y as f64 * cs))
                            .collect()
                    })
                    .collect(),
            });
        }
    }
    out
}

/// Flood fill 4-connected cells of equal value, returning a region id per cell
/// and each region's value.
fn label_regions(raster: &Raster) -> (Vec<usize>, Vec<f64>) {
    let w = raster.width();
    let h = raster.height();
    let mut labels = vec![UNSET; w * h];
    let mut values = Vec::new();

    for start in 0..w * h {
        if labels[start] != UNSET {
            continue;
        }
        let value = raster.get(start / w, start % w).unwrap();
        if raster.is_nodata(value) {
            continue;
        }
        let id = values.len();
        values.push(value);
        labels[start] = id;

        let mut stack = vec![start];
        while let Some(i) = stack.pop() {
            let (r, c) = (i / w, i % w);
            let mut neighbors = [0usize; 4];
            let mut n = 0;
            if r > 0 {
                neighbors[n] = i - w;
                n += 1;
            }
            if r + 1 < h {
                neighbors[n] = i + w;
                n += 1;
            }
            if c > 0 {
                neighbors[n] = i - 1;
                n += 1;
            }
            if c + 1 < w {
                neighbors[n] = i + 1;
                n += 1;
            }
            for &ni in &neighbors[..n] {
                if labels[ni] != UNSET {
                    continue;
                }
                // a nodata neighbour never matches: NaN compares unequal to
                // itself, and a sentinel differs from any region's value
                if raster.get(ni / w, ni % w).unwrap() != value {
                    continue;
                }
                labels[ni] = id;
                stack.push(ni);
            }
        }
    }
    (labels, values)
}

fn direction(a: Vertex, b: Vertex) -> (isize, isize) {
    (b.0 as isize - a.0 as isize, b.1 as isize - a.1 as isize)
}

/// At a pinch, where two diagonal cells of one region meet at a corner, a
/// vertex carries two outgoing edges. Taking the sharpest right turn closes
/// the nearer lobe first, so each lobe becomes its own ring instead of one
/// self-crossing ring.
fn pick_next(cur: Vertex, incoming: Option<(isize, isize)>, nexts: &[Vertex]) -> usize {
    if nexts.len() == 1 {
        return 0;
    }
    let Some((dx, dy)) = incoming else { return 0 };
    for want in [(-dy, dx), (dx, dy), (dy, -dx), (-dx, -dy)] {
        if let Some(i) = nexts.iter().position(|&n| direction(cur, n) == want) {
            return i;
        }
    }
    0
}

/// Walk the directed boundary edges into closed rings.
fn chain_rings(edges: Vec<(Vertex, Vertex)>) -> Vec<Vec<Vertex>> {
    let mut outgoing: HashMap<Vertex, Vec<Vertex>> = HashMap::new();
    for (a, b) in edges {
        outgoing.entry(a).or_default().push(b);
    }

    let mut rings = Vec::new();
    let starts: Vec<Vertex> = outgoing.keys().copied().collect();
    for start in starts {
        while outgoing.get(&start).is_some_and(|v| !v.is_empty()) {
            let mut ring = vec![start];
            let mut cur = start;
            let mut incoming = None;
            let closed = loop {
                let Some(nexts) = outgoing.get_mut(&cur).filter(|v| !v.is_empty()) else {
                    // every boundary vertex has equal in- and out-degree, so
                    // this only trips on a bug; drop the ring rather than
                    // abort the module
                    break false;
                };
                let next = nexts.remove(pick_next(cur, incoming, nexts));
                incoming = Some(direction(cur, next));
                ring.push(next);
                cur = next;
                if cur == start {
                    break true;
                }
            };
            if closed {
                rings.push(drop_collinear(ring));
            } else {
                break;
            }
        }
    }
    rings
}

/// The walk emits one vertex per cell edge, so a straight run along a region's
/// side arrives as a string of collinear points. Keeping only the corners is
/// the difference between five vertices for a rectangle and one per cell.
fn drop_collinear(ring: Vec<Vertex>) -> Vec<Vertex> {
    let n = ring.len() - 1; // the last vertex repeats the first
    if n < 3 {
        return ring;
    }
    let mut kept: Vec<Vertex> = (0..n)
        .filter(|&i| {
            let cur = ring[i];
            direction(ring[(i + n - 1) % n], cur) != direction(cur, ring[(i + 1) % n])
        })
        .map(|i| ring[i])
        .collect();
    if kept.len() < 3 {
        return ring;
    }
    kept.push(kept[0]);
    kept
}

/// Twice the signed area, positive for a clockwise ring in the y-down grid.
fn signed_area2(ring: &[Vertex]) -> f64 {
    let mut sum = 0.0;
    for pair in ring.windows(2) {
        let (x1, y1) = (pair[0].0 as f64, pair[0].1 as f64);
        let (x2, y2) = (pair[1].0 as f64, pair[1].1 as f64);
        sum += x1 * y2 - x2 * y1;
    }
    sum
}

fn contains(ring: &[Vertex], point: Vertex) -> bool {
    let (px, py) = (point.0 as f64, point.1 as f64);
    let mut inside = false;
    for pair in ring.windows(2) {
        let (x1, y1) = (pair[0].0 as f64, pair[0].1 as f64);
        let (x2, y2) = (pair[1].0 as f64, pair[1].1 as f64);
        if (y1 > py) != (y2 > py) && px < (x2 - x1) * (py - y1) / (y2 - y1) + x1 {
            inside = !inside;
        }
    }
    inside
}

/// Split a region's rings into one polygon per outer ring, each carrying the
/// holes it contains. A region normally has a single outer ring; a pinched one
/// has several, and then each lobe is its own polygon.
fn assemble(rings: Vec<Vec<Vertex>>) -> Vec<Vec<Vec<Vertex>>> {
    let (outers, holes): (Vec<_>, Vec<_>) = rings.into_iter().partition(|r| signed_area2(r) > 0.0);
    let mut polygons: Vec<Vec<Vec<Vertex>>> = outers.into_iter().map(|r| vec![r]).collect();

    for hole in holes {
        let probe = hole[0];
        let owner = polygons
            .iter()
            .enumerate()
            .filter(|(_, p)| contains(&p[0], probe))
            .min_by(|(_, a), (_, b)| {
                signed_area2(&a[0])
                    .partial_cmp(&signed_area2(&b[0]))
                    .unwrap()
            })
            .map(|(i, _)| i);
        if let Some(i) = owner {
            polygons[i].push(hole);
        }
    }
    polygons
}

#[cfg(test)]
mod tests {
    use super::*;

    const NODATA: f64 = -9999.0;

    fn raster(width: usize, height: usize, data: Vec<f64>) -> Raster {
        Raster::from_vec(width, height, data, 1.0, NODATA).unwrap()
    }

    fn area(poly: &RegionPolygon) -> f64 {
        poly.rings
            .iter()
            .map(|ring| {
                let mut sum = 0.0;
                for pair in ring.windows(2) {
                    sum += pair[0].0 * pair[1].1 - pair[1].0 * pair[0].1;
                }
                sum / 2.0
            })
            .sum()
    }

    #[test]
    fn one_block_becomes_one_closed_ring() {
        // 2x2 block of 5s inside a 4x4 of 1s
        let mut data = vec![1.0; 16];
        for (r, c) in [(1, 1), (1, 2), (2, 1), (2, 2)] {
            data[r * 4 + c] = 5.0;
        }
        let polys = polygonize(&raster(4, 4, data));

        let block: Vec<_> = polys.iter().filter(|p| p.value == 5.0).collect();
        assert_eq!(block.len(), 1);
        assert_eq!(block[0].rings.len(), 1);
        let ring = &block[0].rings[0];
        assert_eq!(ring.first(), ring.last());
        assert_eq!(ring.len(), 5);
        assert_eq!(area(block[0]), 4.0);
    }

    #[test]
    fn the_surrounding_region_carries_the_block_as_a_hole() {
        let mut data = vec![1.0; 16];
        for (r, c) in [(1, 1), (1, 2), (2, 1), (2, 2)] {
            data[r * 4 + c] = 5.0;
        }
        let polys = polygonize(&raster(4, 4, data));

        let ring1: Vec<_> = polys.iter().filter(|p| p.value == 1.0).collect();
        assert_eq!(ring1.len(), 1);
        assert_eq!(ring1[0].rings.len(), 2);
        // 16 cells less the 4 the hole takes out
        assert_eq!(area(ring1[0]), 12.0);
    }

    #[test]
    fn equal_values_that_do_not_touch_stay_separate() {
        // 7s in opposite corners of a 3x3
        let mut data = vec![1.0; 9];
        data[0] = 7.0;
        data[8] = 7.0;
        let polys = polygonize(&raster(3, 3, data));

        assert_eq!(polys.iter().filter(|p| p.value == 7.0).count(), 2);
    }

    #[test]
    fn nodata_bounds_a_region_instead_of_joining_it() {
        let data = vec![
            NODATA, NODATA, NODATA, //
            NODATA, 3.0, NODATA, //
            NODATA, NODATA, NODATA,
        ];
        let polys = polygonize(&raster(3, 3, data));

        assert_eq!(polys.len(), 1);
        assert_eq!(polys[0].value, 3.0);
        assert_eq!(area(&polys[0]), 1.0);
        // the single cell spans corner (1,1) to (2,2)
        assert!(polys[0].rings[0].contains(&(1.0, 1.0)));
        assert!(polys[0].rings[0].contains(&(2.0, 2.0)));
    }

    #[test]
    fn nan_cells_are_nodata_whatever_the_sentinel_says() {
        let data = vec![f64::NAN, 2.0, 2.0, f64::NAN];
        let polys = polygonize(&raster(4, 1, data));

        assert_eq!(polys.len(), 1);
        assert_eq!(polys[0].value, 2.0);
        assert_eq!(area(&polys[0]), 2.0);
    }

    #[test]
    fn diagonal_neighbours_are_separate_under_four_connectivity() {
        let polys = polygonize(&raster(2, 2, vec![2.0, f64::NAN, f64::NAN, 2.0]));

        assert_eq!(polys.len(), 2);
        assert!(polys.iter().all(|p| p.value == 2.0));
    }

    #[test]
    fn a_straight_run_keeps_only_its_corners() {
        // 4x1 strip: five vertices, not one per cell edge
        let polys = polygonize(&raster(4, 1, vec![9.0; 4]));

        assert_eq!(polys.len(), 1);
        assert_eq!(polys[0].rings[0].len(), 5);
        assert_eq!(area(&polys[0]), 4.0);
    }

    #[test]
    fn ring_coordinates_scale_with_cell_size() {
        let src = Raster::from_vec(2, 1, vec![4.0, 4.0], 30.0, NODATA).unwrap();
        let polys = polygonize(&src);

        assert_eq!(polys.len(), 1);
        assert_eq!(area(&polys[0]), 2.0 * 30.0 * 30.0);
        assert!(polys[0].rings[0].contains(&(60.0, 30.0)));
    }

    #[test]
    fn an_empty_raster_yields_nothing() {
        assert!(polygonize(&raster(2, 2, vec![NODATA; 4])).is_empty());
    }
}
