//! wasm surface over terrano-core: free functions over flat f64 buffers, so a
//! browser hands typed arrays across the boundary without a Raster type on the
//! JS side. Every function takes (data, width, height, cell_size, nodata) plus
//! its own parameters and returns a new buffer of the same shape. Contours
//! come back flat-encoded, [level, vertex_count, x0, y0, x1, y1, ...] per
//! line, because nested structs would otherwise be serialized per vertex.
//! Polygons nest one level deeper: [value, ring_count, (vertex_count, x0, y0,
//! ...) per ring], exterior ring first.

use terrano_core::{BinaryOp, Raster, UnaryOp};
use wasm_bindgen::prelude::*;

fn raster(
    data: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
) -> Result<Raster, JsError> {
    Raster::from_vec(width, height, data.to_vec(), cell_size, nodata)
        .map_err(|e| JsError::new(&e.to_string()))
}

#[wasm_bindgen]
pub fn hillshade(
    data: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
    azimuth: f64,
    altitude: f64,
) -> Result<Vec<f64>, JsError> {
    let dem = raster(data, width, height, cell_size, nodata)?;
    Ok(terrano_core::hillshade(&dem, azimuth, altitude).into_data())
}

#[wasm_bindgen]
pub fn slope(
    data: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
) -> Result<Vec<f64>, JsError> {
    let dem = raster(data, width, height, cell_size, nodata)?;
    Ok(terrano_core::slope(&dem).into_data())
}

#[wasm_bindgen]
pub fn aspect(
    data: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
) -> Result<Vec<f64>, JsError> {
    let dem = raster(data, width, height, cell_size, nodata)?;
    Ok(terrano_core::aspect(&dem).into_data())
}

#[wasm_bindgen(js_name = fillSinks)]
pub fn fill_sinks(
    data: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
) -> Result<Vec<f64>, JsError> {
    let dem = raster(data, width, height, cell_size, nodata)?;
    Ok(terrano_core::fill_sinks(&dem).into_data())
}

/// `classes` is flat (min_inclusive, max_exclusive, new_value) triples.
#[wasm_bindgen]
pub fn reclassify(
    data: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
    classes: &[f64],
) -> Result<Vec<f64>, JsError> {
    if classes.len() % 3 != 0 {
        return Err(JsError::new(
            "classes must be flat (min, max, value) triples",
        ));
    }
    let triples: Vec<(f64, f64, f64)> = classes.chunks(3).map(|c| (c[0], c[1], c[2])).collect();
    let src = raster(data, width, height, cell_size, nodata)?;
    Ok(terrano_core::reclassify(&src, &triples).into_data())
}

fn unary_op(op: &str, operand: f64) -> Result<UnaryOp, JsError> {
    Ok(match op {
        "add" => UnaryOp::Add(operand),
        "multiply" => UnaryOp::Multiply(operand),
        "sqrt" => UnaryOp::Sqrt,
        "abs" => UnaryOp::Abs,
        "log" => UnaryOp::Log,
        _ => return Err(JsError::new(&format!("unknown unary op {op:?}"))),
    })
}

#[wasm_bindgen(js_name = applyUnary)]
pub fn apply_unary(
    data: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
    op: &str,
    operand: f64,
) -> Result<Vec<f64>, JsError> {
    let src = raster(data, width, height, cell_size, nodata)?;
    Ok(src.apply_unary(&unary_op(op, operand)?).into_data())
}

fn binary_op(op: &str) -> Result<BinaryOp, JsError> {
    Ok(match op {
        "add" => BinaryOp::Add,
        "subtract" => BinaryOp::Subtract,
        "multiply" => BinaryOp::Multiply,
        "divide" => BinaryOp::Divide,
        "min" => BinaryOp::Min,
        "max" => BinaryOp::Max,
        _ => return Err(JsError::new(&format!("unknown binary op {op:?}"))),
    })
}

#[wasm_bindgen(js_name = applyBinary)]
pub fn apply_binary(
    a: &[f64],
    b: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
    op: &str,
) -> Result<Vec<f64>, JsError> {
    let left = raster(a, width, height, cell_size, nodata)?;
    let right = raster(b, width, height, cell_size, nodata)?;
    Ok(left
        .apply_binary(&right, &binary_op(op)?)
        .map_err(|e| JsError::new(&e.to_string()))?
        .into_data())
}

/// (a - b) / (a + b), the NDVI/NDWI family, composed from the same binary ops
/// a calculator would use so nodata and zero-sum cells behave identically.
#[wasm_bindgen(js_name = normalizedDifference)]
pub fn normalized_difference(
    a: &[f64],
    b: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
) -> Result<Vec<f64>, JsError> {
    let left = raster(a, width, height, cell_size, nodata)?;
    let right = raster(b, width, height, cell_size, nodata)?;
    let err = |e: terrano_core::Error| JsError::new(&e.to_string());
    let num = left
        .apply_binary(&right, &BinaryOp::Subtract)
        .map_err(err)?;
    let den = left.apply_binary(&right, &BinaryOp::Add).map_err(err)?;
    Ok(num
        .apply_binary(&den, &BinaryOp::Divide)
        .map_err(err)?
        .into_data())
}

#[wasm_bindgen]
pub fn contours(
    data: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
    interval: f64,
    base: f64,
) -> Result<Vec<f64>, JsError> {
    let dem = raster(data, width, height, cell_size, nodata)?;
    let lines = terrano_core::contours(&dem, interval, base);
    let mut flat = Vec::with_capacity(lines.iter().map(|l| 2 + l.vertices.len() * 2).sum());
    for line in lines {
        flat.push(line.level);
        flat.push(line.vertices.len() as f64);
        for (x, y) in line.vertices {
            flat.push(x);
            flat.push(y);
        }
    }
    Ok(flat)
}

#[wasm_bindgen]
pub fn polygonize(
    data: &[f64],
    width: usize,
    height: usize,
    cell_size: f64,
    nodata: f64,
) -> Result<Vec<f64>, JsError> {
    let src = raster(data, width, height, cell_size, nodata)?;
    let polygons = terrano_core::polygonize(&src);
    let mut flat = Vec::new();
    for polygon in polygons {
        flat.push(polygon.value);
        flat.push(polygon.rings.len() as f64);
        for ring in polygon.rings {
            flat.push(ring.len() as f64);
            for (x, y) in ring {
                flat.push(x);
                flat.push(y);
            }
        }
    }
    Ok(flat)
}

// happy paths only: constructing a JsError aborts off wasm, so the error arms
// are exercised by the browser tests in the consuming viewer instead
#[cfg(test)]
mod tests {
    use super::*;

    const NODATA: f64 = -9999.0;

    fn ramp(width: usize, height: usize) -> Vec<f64> {
        (0..width * height)
            .map(|i| (i / width) as f64 * 10.0)
            .collect()
    }

    #[test]
    fn hillshade_matches_core_and_keeps_shape() {
        let dem = ramp(5, 5);
        let out = hillshade(&dem, 5, 5, 30.0, NODATA, 315.0, 45.0).unwrap();
        assert_eq!(out.len(), 25);
        let core = terrano_core::hillshade(
            &Raster::from_vec(5, 5, dem, 30.0, NODATA).unwrap(),
            315.0,
            45.0,
        );
        assert_eq!(out, core.into_data());
        // interior of a north-south ramp under a north-west sun is lit
        assert!(out[12] > 0.0 && out[12] <= 255.0);
    }

    #[test]
    fn slope_of_a_flat_raster_is_zero_inside_the_border() {
        let out = slope(&[7.0; 16], 4, 4, 10.0, NODATA).unwrap();
        assert_eq!(out[5], 0.0);
        assert_eq!(out[0], NODATA, "border cells have no gradient");
    }

    #[test]
    fn normalized_difference_is_the_ndvi_formula() {
        let nir = [0.8, 0.5, NODATA, 0.0];
        let red = [0.2, 0.5, 0.1, 0.0];
        let out = normalized_difference(&nir, &red, 2, 2, 1.0, NODATA).unwrap();
        assert!((out[0] - 0.6).abs() < 1e-12);
        assert_eq!(out[1], 0.0);
        assert_eq!(out[2], NODATA, "nodata propagates");
        assert_eq!(out[3], NODATA, "zero denominator is nodata, not inf");
    }

    #[test]
    fn reclassify_takes_flat_triples() {
        let out = reclassify(
            &[1.0, 5.0, 9.0, NODATA],
            2,
            2,
            1.0,
            NODATA,
            &[
                0.0, 4.0, 100.0, //
                4.0, 8.0, 200.0,
            ],
        )
        .unwrap();
        assert_eq!(&out[..2], &[100.0, 200.0]);
        assert_eq!(out[2], NODATA, "a value outside every class is nodata");
        assert_eq!(out[3], NODATA);
    }

    #[test]
    fn contours_flat_encode_level_count_then_vertices() {
        // two flat halves meeting in the middle: levels at 5 and at the top
        let mut dem = vec![0.0; 20];
        dem[10..].fill(10.0);
        let flat = contours(&dem, 5, 4, 1.0, NODATA, 5.0, 0.0).unwrap();
        let mut i = 0;
        let mut levels = Vec::new();
        while i < flat.len() {
            levels.push(flat[i]);
            let n = flat[i + 1] as usize;
            assert!(n >= 2, "a line has at least two vertices");
            i += 2 + n * 2;
        }
        assert_eq!(i, flat.len(), "flat encoding parses exactly");
        assert!(levels.contains(&5.0), "{levels:?}");
    }

    #[test]
    fn polygonize_flat_encodes_value_rings_then_vertices() {
        // a 5 block inside a 1 field, so one polygon carries a hole
        let mut data = vec![1.0; 16];
        for i in [5, 6, 9, 10] {
            data[i] = 5.0;
        }
        let flat = polygonize(&data, 4, 4, 1.0, NODATA).unwrap();

        let mut i = 0;
        let mut shapes = Vec::new();
        while i < flat.len() {
            let value = flat[i];
            let rings = flat[i + 1] as usize;
            i += 2;
            for _ in 0..rings {
                let n = flat[i] as usize;
                assert!(n >= 4, "a ring closes on at least four vertices");
                i += 1 + n * 2;
            }
            shapes.push((value, rings));
        }
        assert_eq!(i, flat.len(), "flat encoding parses exactly");
        assert!(shapes.contains(&(5.0, 1)), "{shapes:?}");
        assert!(shapes.contains(&(1.0, 2)), "{shapes:?}");
    }

    #[test]
    fn binary_and_unary_ops_parse_their_names() {
        let a = [4.0, 9.0];
        let b = [2.0, 3.0];
        let div = apply_binary(&a, &b, 2, 1, 1.0, NODATA, "divide").unwrap();
        assert_eq!(div, vec![2.0, 3.0]);
        let sqrt = apply_unary(&a, 2, 1, 1.0, NODATA, "sqrt", 0.0).unwrap();
        assert_eq!(sqrt, vec![2.0, 3.0]);
    }
}
