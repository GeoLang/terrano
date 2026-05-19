use clap::{Parser, Subcommand};
use terrano_core::{Raster, aspect, hillshade, slope};

#[derive(Parser)]
#[command(
    name = "terrano",
    version,
    about = "Raster algebra and terrain analysis CLI"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compute summary statistics of a synthetic DEM
    Stats {
        /// Width of synthetic DEM
        #[arg(long, default_value_t = 10)]
        width: usize,
        /// Height of synthetic DEM
        #[arg(long, default_value_t = 10)]
        height: usize,
    },
    /// Compute hillshade of a synthetic DEM
    Hillshade {
        /// Sun azimuth in degrees
        #[arg(long, default_value_t = 315.0)]
        azimuth: f64,
        /// Sun altitude in degrees
        #[arg(long, default_value_t = 45.0)]
        altitude: f64,
    },
}

fn synthetic_dem(width: usize, height: usize) -> Raster {
    let mut data = vec![0.0; width * height];
    for row in 0..height {
        for col in 0..width {
            data[row * width + col] = (row as f64) * 10.0 + (col as f64) * 5.0;
        }
    }
    Raster::from_vec(width, height, data, 1.0, -9999.0).unwrap()
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::Stats { width, height } => {
            let dem = synthetic_dem(width, height);
            let s = slope(&dem);
            let a = aspect(&dem);
            // Report interior cell stats
            let mid_r = height / 2;
            let mid_c = width / 2;
            println!(
                "DEM center elevation: {:.2}",
                dem.get(mid_r, mid_c).unwrap()
            );
            println!("Slope at center: {:.2}°", s.get(mid_r, mid_c).unwrap());
            println!("Aspect at center: {:.2}°", a.get(mid_r, mid_c).unwrap());
        }
        Commands::Hillshade { azimuth, altitude } => {
            let dem = synthetic_dem(10, 10);
            let hs = hillshade(&dem, azimuth, altitude);
            println!("Hillshade at center: {:.2}", hs.get(5, 5).unwrap());
        }
    }
}
