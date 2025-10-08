pub mod air_pls;
pub use air_pls::air_pls;

pub mod calculate_eic;
pub use calculate_eic::{Eic, EicOptions, calculate_eic_from_bin1, calculate_eic_from_mzml};

pub mod find_noise_level;
pub use find_noise_level::find_noise_level;

pub mod find_peaks;

pub mod functions;

pub mod calculate_baseline;
pub use calculate_baseline::calculate_baseline;

pub mod get_boundaries;

pub mod get_peak;

pub mod get_peaks_from_chrom;

pub mod get_peaks_from_eic;

pub mod kmeans;
pub use kmeans::{Point, kmeans};

pub mod lm;
pub use lm::lm;

pub mod noise_san_plot;
pub use noise_san_plot::noise_san_plot;

pub mod parse;

pub mod scan_for_peaks;

pub mod sgg;

pub mod structs;

pub mod utilities;
pub use utilities::{
    closest_index, mean_step, min_positive_step, min_sep, odd_in_range, quad_peak,
};
