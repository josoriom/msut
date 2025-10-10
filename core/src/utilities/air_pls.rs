//! airPLS (Adaptive Iteratively Reweighted Penalized Least Squares) — baseline correction
//!
//! Rust implementation of JavaScript version from:
//! https://github.com/mljs/airpls
//!
//! # References
//! * Z. M. Zhang, S. Chen, Y. Z. Liang (2010).
//! “Baseline correction using adaptive iteratively reweighted penalized least squares.”
//!   https://doi.org/10.1039/B922045C
//! * MLJS `airpls` repository.

use crate::utilities::noise_san_plot::{NoiseSanPlotOptions, noise_san_plot};

#[derive(Clone, Copy, Debug)]
pub struct Zone {
    pub from: f64,
    pub to: f64,
}

#[derive(Clone, Debug)]
pub struct AirPlsOptions {
    pub max_iterations: usize,
    pub lambda: f64,
    pub tolerance: f64,
    pub control_points: Option<Vec<u8>>,
    pub zones: Option<Vec<Zone>>,
    pub weights: Option<Vec<f64>>,
}

impl Default for AirPlsOptions {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            lambda: 10.0,
            tolerance: 0.01,
            control_points: None,
            zones: None,
            weights: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct AirPlsResult {
    pub corrected: Vec<f64>,
    pub baseline: Vec<f64>,
    pub iteration: usize,
    pub error: f64,
}

pub fn air_pls(x: &[f64], y: &[f64], opts: AirPlsOptions) -> AirPlsResult {
    let n = y.len();
    assert_eq!(x.len(), n);

    let (mut weights, control_points) = get_control_points(x, y, &opts);
    let mut baseline = vec![0.0; n];
    let mut corrected = y.to_vec();

    let mut iteration = 0usize;
    let mut sum_neg_differences = f64::MAX;
    let mut prev_neg_sum = f64::MAX;
    let stop_criterion = get_stop_criterion(y, opts.tolerance);

    let mut positive_threshold = 1.0;

    while iteration < opts.max_iterations && sum_neg_differences.abs() > stop_criterion {
        let rhs: Vec<f64> = y.iter().zip(weights.iter()).map(|(yy, w)| yy * w).collect();
        baseline = solve_tridiagonal_system(&weights, opts.lambda, &rhs);

        sum_neg_differences = apply_correction(y, &baseline, &mut corrected);

        if iteration == 1 {
            let noise_result = noise_san_plot(&corrected, NoiseSanPlotOptions::default());
            positive_threshold = noise_result.positive;
        } else {
            let abs_change = (prev_neg_sum / sum_neg_differences).abs();
            if abs_change < 1.01 && abs_change > 0.99 {
                break;
            }
        }

        prev_neg_sum = sum_neg_differences;

        let denom = sum_neg_differences.abs().max(1e-16);
        let last = n - 1;

        for i in 1..last {
            let diff = corrected[i];
            if control_points[i] < 1 && diff.abs() > positive_threshold {
                weights[i] = 0.0;
            } else {
                let factor = if diff > 0.0 { -1.0 } else { 1.0 };
                weights[i] = ((factor * (iteration as f64 * diff)) / denom).exp();
            }
        }
        weights[0] = 1.0;
        weights[last] = 1.0;

        iteration += 1;
    }

    AirPlsResult {
        corrected,
        baseline,
        iteration,
        error: sum_neg_differences,
    }
}

fn get_control_points(x: &[f64], _y: &[f64], opts: &AirPlsOptions) -> (Vec<f64>, Vec<u8>) {
    let n = x.len();
    let mut control_points = opts.control_points.clone().unwrap_or_else(|| vec![0u8; n]);

    let mut weights = opts.weights.clone().unwrap_or_else(|| vec![1.0; n]);

    if let Some(zones) = &opts.zones {
        for z in zones {
            let mut index_from = get_close_index(x, z.from);
            let mut index_to = get_close_index(x, z.to);
            if index_from > index_to {
                std::mem::swap(&mut index_from, &mut index_to);
            }
            for i in index_from..index_to.min(n) {
                control_points[i] = 1;
            }
        }
    }

    if opts.control_points.is_some() || opts.zones.is_some() {
        for i in 0..n {
            weights[i] *= control_points[i] as f64;
        }
    }

    (weights, control_points)
}

fn get_stop_criterion(y: &[f64], tolerance: f64) -> f64 {
    y.iter().map(|v| v.abs()).sum::<f64>() * tolerance
}

fn get_close_index(xs: &[f64], goal: f64) -> usize {
    let mut best = 0usize;
    let mut best_d = (xs[0] - goal).abs();
    for i in 1..xs.len() {
        let d = (xs[i] - goal).abs();
        if d < best_d {
            best = i;
            best_d = d;
        }
    }
    best
}

fn apply_correction(y: &[f64], baseline: &[f64], corrected: &mut [f64]) -> f64 {
    let mut sum_neg_differences = 0.0;
    for i in 0..y.len() {
        let diff = y[i] - baseline[i];
        if diff < 0.0 {
            sum_neg_differences += diff;
        }
        corrected[i] = diff;
    }
    sum_neg_differences
}

fn solve_tridiagonal_system(weights: &[f64], lambda: f64, rhs: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    let mut a = vec![0.0; n];
    let mut b = vec![0.0; n];
    let mut c = vec![0.0; n];

    if n == 1 {
        b[0] = lambda + weights[0];
        return vec![rhs[0] / b[0]];
    }

    b[0] = lambda + weights[0];
    c[0] = -lambda;

    for i in 1..n - 1 {
        a[i] = -lambda;
        b[i] = 2.0 * lambda + weights[i];
        c[i] = -lambda;
    }

    a[n - 1] = -lambda;
    b[n - 1] = lambda + weights[n - 1];

    let mut cp = vec![0.0; n - 1];
    let mut dp = vec![0.0; n];
    let mut denom = b[0];
    cp[0] = c[0] / denom;
    dp[0] = rhs[0] / denom;

    for i in 1..n {
        denom = b[i] - a[i] * cp[i - 1];
        if i < n - 1 {
            cp[i] = c[i] / denom;
        }
        dp[i] = (rhs[i] - a[i] * dp[i - 1]) / denom;
    }

    let mut x = vec![0.0; n];
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    x
}
