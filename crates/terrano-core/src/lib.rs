//! Terrano — Raster algebra and terrain analysis engine.
//!
//! Grid-based raster operations including map algebra, hillshade, slope, aspect,
//! and hydrological analysis.

mod algebra;
mod banded;
mod cog;
mod contour;
mod error;
mod focal;
mod geotiff;
pub mod grib;
mod hydrology;
pub mod netcdf;
mod polygonize;
mod raster;
mod rasterize;
mod terrain;
mod timeseries;
mod watershed;
mod zonal;

pub use algebra::{BinaryOp, UnaryOp, reclassify};
pub use banded::BandedRaster;
pub use cog::{
    CogLevel, CogParams, CogReader, Overview, RangeRead, extract_tile, generate_overviews,
    write_cog, write_cog_bands,
};
pub use contour::{ContourLine, ContourSegment, contours, fill_sinks};
pub use error::Error;
pub use focal::{FocalStat, Neighborhood, focal_stats};
pub use geotiff::{
    GeoTiffMetadata, SampleFormat, read_geotiff, read_geotiff_bands, write_geotiff,
    write_geotiff_bands,
};
pub use hydrology::{flow_accumulation, flow_direction};
pub use polygonize::{RegionPolygon, polygonize};
pub use raster::Raster;
pub use rasterize::rasterize;
pub use terrain::{aspect, hillshade, slope};
pub use timeseries::{ChangeResult, CompositeMethod, PhenologyMetrics, RasterStack, TrendResult};
pub use watershed::{stream_order, watershed};
pub use zonal::{ZoneStats, zonal_stats};
