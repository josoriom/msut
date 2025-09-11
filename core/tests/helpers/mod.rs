// tests/helpers.rs
use msut::utilities::structs::DataXY;

/// Get one Gaussian value at x (bell curve). `base` lifts the whole curve.
#[inline]
pub fn gaussian_value(x: f64, mu: f64, sigma: f64, amp: f64, base: f64) -> f64 {
    base + amp * (-0.5 * ((x - mu) / sigma).powi(2)).exp()
}

/// Build a signal made of one or more Gaussians, plus a constant `base` and optional small fake `noise`. Returns f32 values for plotting/algorithms.
pub fn gaussian_mixture_f32(
    xs: &[f64],
    peaks: &[(f64, f64, f64)],
    base: f64,
    noise: f64,
) -> Vec<f32> {
    xs.iter()
        .map(|&x| {
            let mut y = base;
            for &(mu, sigma, amp) in peaks {
                y += gaussian_value(x, mu, sigma, amp, 0.0);
            }
            if noise > 0.0 {
                let z = ((x * 137.13).sin() + (x * 73.7).cos()) * 0.5;
                y += z * noise;
            }
            y as f32
        })
        .collect()
}

/// Make an even grid from start to end with n points (inclusive).
pub fn make_grid(start: f64, end: f64, n: usize) -> Vec<f64> {
    if n <= 1 {
        return vec![start];
    }
    (0..n)
        .map(|i| start + (end - start) * (i as f64) / ((n - 1) as f64))
        .collect()
}

/// Same as `make_grid`, named like NumPy for familiarity.
#[allow(dead_code)]
pub fn linspace(from: f64, to: f64, n: usize) -> Vec<f64> {
    make_grid(from, to, n)
}

/// Small deterministic wiggle in [-0.5, 0.5] you can use as repeatable “noise”.
#[allow(dead_code)]
pub fn jitter(i: u32) -> f64 {
    let mut x = i.wrapping_mul(1664525).wrapping_add(1013904223);
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    (x as f64 / (u32::MAX as f64)) - 0.5
}

/// Check two numbers are close within a simple absolute tolerance.
#[allow(dead_code)]
pub fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

/// Build a `DataXY` from x and y vectors.
#[allow(dead_code)]
pub fn data_xy(xs: Vec<f64>, ys: Vec<f32>) -> DataXY {
    DataXY { x: xs, y: ys }
}

// Deterministic uniform values in [lo, hi]
#[allow(dead_code)]
pub fn uniform_vec_f32(n: usize, lo: f32, hi: f32, seed: u64) -> Vec<f32> {
    assert!(hi > lo);
    let mut out = Vec::with_capacity(n);
    let mut s = seed | 1; // odd
    for _ in 0..n {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let u = ((s >> 11) as f64) * (1.0 / (1u64 << 53) as f64);
        let v = lo as f64 + (hi - lo) as f64 * u;
        out.push(v as f32);
    }
    out
}

// Deterministic shuffle (Fisher–Yates)
#[allow(dead_code)]
pub fn shuffle_with_seed<T>(xs: &mut [T], seed: u64) {
    let mut s = seed | 1;
    let n = xs.len();
    if n <= 1 {
        return;
    }
    for i in (1..n).rev() {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let u = ((s >> 11) as f64) * (1.0 / (1u64 << 53) as f64);
        let j = (u * ((i + 1) as f64)).floor() as usize;
        xs.swap(i, j);
    }
}
