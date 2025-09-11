use msut::utilities::sgg::{SggOptions, sgg};
use std::f64::consts::PI;

fn assert_close(actual: f64, expected: f64, decimals: i32) {
    let tol = 10f64.powi(-decimals.max(0));
    assert!(
        (actual - expected).abs() <= tol,
        "expected {actual} ≈ {expected} (±{tol})"
    );
}

fn assert_close_tol(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol,
        "expected {actual} ≈ {expected} (±{tol})"
    );
}

fn jitter(i: u32) -> f64 {
    let mut x = i.wrapping_mul(1664525).wrapping_add(1013904223);
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    (x as f64 / (u32::MAX as f64)) - 0.5
}

#[test]
fn simple_triangle_check_symmetry_with_x_array() {
    let mut xs = Vec::with_capacity(101);
    for i in 0..=100 {
        xs.push(i as f64);
    }

    let mut ys = vec![0.0f32; 101];
    let last = ys.len() - 1;
    for i in 0..=50 {
        let v = i as f32;
        ys[i] = v;
        ys[last - i] = v;
    }

    let out = sgg(
        &ys,
        &xs,
        SggOptions {
            window_size: 9,
            derivative: 0,
            polynomial: 3,
        },
    );

    for i in 0..50 {
        let a = out[i] as f64;
        let b = out[last - i] as f64;
        assert_close(a, b, 2);
    }
}

#[test]
fn smoothing_test() {
    let opts = SggOptions {
        window_size: 15,
        derivative: 0,
        polynomial: 3,
    };
    let n = 200usize;
    let noise = 0.01;

    let dx = 2.0 * PI / n as f64;
    let mut xs = Vec::with_capacity(n);
    for i in 0..n {
        xs.push(i as f64 * dx);
    }

    let mut data = vec![0.0f32; n];
    for i in 0..n {
        let angle = i as f64 * dx;
        let v = angle.sin() + jitter(i as u32) * noise;
        data[i] = v as f32;
    }

    let out = sgg(&data, &xs, opts);
    let half = opts.window_size / 2;

    let decimals = ((-noise.log10()) - 1.0).round() as i32;
    let tol = 10f64.powi(-decimals.max(0));

    for j in half..(n - half) {
        assert_close_tol(out[j] as f64, data[j] as f64, tol);
    }
}

#[test]
fn first_derivative_test() {
    let opts = SggOptions {
        window_size: 47,
        derivative: 1,
        polynomial: 3,
    };
    let noise = 0.1;
    let n = 200usize;
    let dx = 2.0 * PI / n as f64;

    let mut xs = Vec::with_capacity(n);
    for i in 0..n {
        xs.push(i as f64 * dx);
    }

    let mut data = vec![0.0f32; n];
    for i in 0..n {
        let angle = i as f64 * dx;
        let v = angle.sin() + jitter(i as u32) * noise;
        data[i] = v as f32;
    }

    let out = sgg(&data, &xs, opts);
    let half = opts.window_size / 2;

    let decimals = ((-noise.log10()) - 1.0).round() as i32;
    let tol = 10f64.powi(-decimals.max(0));

    for j in half..(n - half) {
        let expected = (j as f64 * dx).cos();
        assert_close_tol(out[j] as f64, expected, tol);
    }
}

#[test]
fn first_derivative_x_as_vector_equivalence() {
    let opts = SggOptions {
        window_size: 47,
        derivative: 1,
        polynomial: 3,
    };
    let n = 200usize;
    let dx = 2.0 * PI / n as f64;

    let mut xs = Vec::with_capacity(n);
    for i in 0..n {
        xs.push(i as f64 * dx);
    }

    let mut data = Vec::with_capacity(n);
    for i in 0..n {
        data.push((i as f64 * dx).sin() as f32);
    }

    let a = sgg(&data, &xs, opts);
    let b = sgg(&data, &xs, opts);

    let half = opts.window_size / 2;
    for j in half..(n - half) {
        let diff = (a[j] as f64) - (b[j] as f64);
        assert!(diff.abs() <= 1e-6);
    }
}

#[test]
fn border_test() {
    let opts = SggOptions {
        window_size: 9,
        derivative: 1,
        polynomial: 3,
    };
    let n = 20usize;

    let mut xs = Vec::with_capacity(n);
    for i in 0..n {
        xs.push(i as f64);
    }

    let mut data = Vec::with_capacity(n);
    for i in 0..n {
        let x = i as f64;
        data.push((x * x * x - 4.0 * x * x + 5.0 * x) as f32);
    }

    let out = sgg(&data, &xs, opts);

    for j in 0..n {
        let x = j as f64;
        let expected = (3.0 * x - 8.0) * x + 5.0;
        assert_close(out[j] as f64, expected, 4);
    }
}
