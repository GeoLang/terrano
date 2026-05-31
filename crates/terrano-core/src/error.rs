use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("dimension mismatch: expected {expected} cells, got {got}")]
    DimensionMismatch { expected: usize, got: usize },

    #[error("format error: {0}")]
    Format(String),

    #[error("invalid input: {0}")]
    InvalidInput(String),

    #[error("rasters have incompatible dimensions")]
    IncompatibleRasters,

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}
