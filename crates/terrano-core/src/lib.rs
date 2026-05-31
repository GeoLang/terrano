//! Terrano — Raster algebra and terrain analysis engine.
//!
//! Grid-based raster operations including map algebra, hillshade, slope, aspect,
//! and hydrological analysis.

mod algebra;
mod cog;
mod contour;
mod error;
mod geotiff;
pub mod grib;
mod hydrology;
pub mod netcdf;
mod raster;
mod terrain;
mod timeseries;
mod watershed;

pub use algebra::{BinaryOp, UnaryOp, reclassify};
pub use cog::{CogParams, Overview, extract_tile, generate_overviews, write_cog};
pub use contour::{ContourLine, ContourSegment, contours, fill_sinks};
pub use error::Error;
pub use geotiff::{GeoTiffMetadata, read_geotiff, write_geotiff};
pub use hydrology::{flow_accumulation, flow_direction};
pub use raster::Raster;
pub use terrain::{aspect, hillshade, slope};
pub use timeseries::{ChangeResult, CompositeMethod, PhenologyMetrics, RasterStack, TrendResult};
pub use watershed::{stream_order, watershed};
