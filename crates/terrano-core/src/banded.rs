//! Multi-band rasters, e.g. an RGB or RGBA image.

use crate::{Error, Raster};

/// An ordered set of bands sharing one grid: width, height, cell size and georeferencing.
///
/// Each band is a plain [`Raster`], so every single-band algorithm works per band.
#[derive(Debug, Clone)]
pub struct BandedRaster {
    bands: Vec<Raster>,
    names: Vec<Option<String>>,
}

impl BandedRaster {
    /// Create from unnamed bands. All bands must have the same dimensions.
    pub fn new(bands: Vec<Raster>) -> Result<Self, Error> {
        let names = vec![None; bands.len()];
        Self::assemble(bands, names)
    }

    /// Create from named bands, e.g. `["red", "green", "blue"]`.
    pub fn with_names(bands: Vec<Raster>, names: Vec<String>) -> Result<Self, Error> {
        if names.len() != bands.len() {
            return Err(Error::InvalidInput(format!(
                "{} band names for {} bands",
                names.len(),
                bands.len()
            )));
        }
        Self::assemble(bands, names.into_iter().map(Some).collect())
    }

    fn assemble(bands: Vec<Raster>, names: Vec<Option<String>>) -> Result<Self, Error> {
        let first = bands
            .first()
            .ok_or_else(|| Error::InvalidInput("BandedRaster requires at least one band".into()))?;
        let (w, h) = (first.width(), first.height());
        for (i, band) in bands.iter().enumerate().skip(1) {
            if band.width() != w || band.height() != h {
                return Err(Error::InvalidInput(format!(
                    "band {i} is {}x{}, first band is {w}x{h}",
                    band.width(),
                    band.height()
                )));
            }
        }
        Ok(Self { bands, names })
    }

    pub fn band_count(&self) -> usize {
        self.bands.len()
    }

    pub fn width(&self) -> usize {
        self.bands[0].width()
    }

    pub fn height(&self) -> usize {
        self.bands[0].height()
    }

    pub fn cell_size(&self) -> f64 {
        self.bands[0].cell_size
    }

    pub fn band(&self, index: usize) -> Option<&Raster> {
        self.bands.get(index)
    }

    pub fn band_mut(&mut self, index: usize) -> Option<&mut Raster> {
        self.bands.get_mut(index)
    }

    pub fn bands(&self) -> &[Raster] {
        &self.bands
    }

    pub fn band_name(&self, index: usize) -> Option<&str> {
        self.names.get(index)?.as_deref()
    }

    pub fn into_bands(self) -> Vec<Raster> {
        self.bands
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn band(width: usize, height: usize, fill: f64) -> Raster {
        Raster::from_vec(width, height, vec![fill; width * height], 1.0, -9999.0).unwrap()
    }

    #[test]
    fn test_banded_new() {
        let b = BandedRaster::new(vec![band(3, 2, 1.0), band(3, 2, 2.0)]).unwrap();
        assert_eq!(b.band_count(), 2);
        assert_eq!(b.width(), 3);
        assert_eq!(b.height(), 2);
        assert_eq!(b.cell_size(), 1.0);
        assert_eq!(b.band(1).unwrap().get(1, 2), Some(2.0));
        assert!(b.band(2).is_none());
        assert_eq!(b.band_name(0), None);
    }

    #[test]
    fn test_banded_mismatched_dimensions_is_error() {
        let result = BandedRaster::new(vec![band(3, 2, 1.0), band(2, 2, 1.0)]);
        assert!(result.is_err());
    }

    #[test]
    fn test_banded_empty_is_error() {
        assert!(BandedRaster::new(vec![]).is_err());
    }

    #[test]
    fn test_banded_with_names() {
        let b = BandedRaster::with_names(
            vec![band(2, 2, 1.0), band(2, 2, 2.0), band(2, 2, 3.0)],
            vec!["red".into(), "green".into(), "blue".into()],
        )
        .unwrap();
        assert_eq!(b.band_name(0), Some("red"));
        assert_eq!(b.band_name(2), Some("blue"));
        assert_eq!(b.band_name(3), None);
    }

    #[test]
    fn test_banded_name_count_mismatch_is_error() {
        let result = BandedRaster::with_names(vec![band(2, 2, 1.0)], vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn test_banded_band_mut() {
        let mut b = BandedRaster::new(vec![band(2, 2, 0.0)]).unwrap();
        b.band_mut(0).unwrap().set(0, 1, 7.0);
        assert_eq!(b.band(0).unwrap().get(0, 1), Some(7.0));
    }
}
