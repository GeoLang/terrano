//! Terrano — Raster algebra and terrain analysis engine.
//!
//! Grid-based raster operations including map algebra, hillshade, slope, aspect,
//! and hydrological analysis.

mod algebra;
mod error;
mod hydrology;
mod raster;
mod terrain;

pub use algebra::{BinaryOp, UnaryOp, reclassify};
pub use error::Error;
pub use hydrology::{flow_accumulation, flow_direction};
pub use raster::Raster;
pub use terrain::{aspect, hillshade, slope};
