use crate::utilities::air_pls::{AirPlsOptions, air_pls};

#[derive(Clone, Copy, Debug)]
pub struct BaselineOptions {
    pub baseline_window: Option<f64>,
    pub baseline_window_factor: Option<usize>,
    pub level: Option<u8>,
}

impl Default for BaselineOptions {
    fn default() -> Self {
        Self {
            baseline_window: Some(10.0),
            baseline_window_factor: Some(100),
            level: Some(1),
        }
    }
}

pub fn calculate_baseline(y: &[f64], options: BaselineOptions) -> Vec<f64> {
    let n = y.len();
    if n == 0 {
        return Vec::new();
    }

    let x: Vec<f64> = (0..n).map(|i| i as f64).collect();

    let defaults = BaselineOptions::default();
    let window = options
        .baseline_window
        .unwrap_or(defaults.baseline_window.unwrap_or(10.0));
    let factor = options
        .baseline_window_factor
        .unwrap_or(defaults.baseline_window_factor.unwrap_or(100));

    let result = air_pls(
        &x,
        y,
        AirPlsOptions {
            lambda: window,
            max_iterations: factor,
            ..Default::default()
        },
    );
    let baseline = result.baseline;

    baseline
}
