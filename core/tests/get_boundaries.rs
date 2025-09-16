mod helpers;
use helpers::{data_xy, gaussian_mixture_f32, jitter, linspace};

use msut::utilities::get_boundaries::{BoundariesOptions, get_boundaries};
use msut::utilities::structs::DataXY;

#[test]
fn empty_mismatch_single_return_none() {
    let d = DataXY {
        x: vec![],
        y: vec![],
    };
    let b = get_boundaries(&d, 0.0, None);
    assert!(b.from.index.is_none() && b.to.index.is_none());

    let d = DataXY {
        x: vec![0.0, 1.0],
        y: vec![1.0],
    };
    let b = get_boundaries(&d, 0.5, None);
    assert!(b.from.index.is_none() && b.to.index.is_none());

    let d = DataXY {
        x: vec![0.0],
        y: vec![1.0],
    };
    let b = get_boundaries(&d, 0.0, None);
    assert!(b.from.index.is_none() && b.to.index.is_none());
}

#[test]
fn small_noise_with_upward_baseline_triggers_stop() {
    let n = 601;
    let xs = linspace(0.0, 12.0, n);
    let mu = 6.0;
    let mut ys = gaussian_mixture_f32(&xs, &[(mu, 0.6, 120.0)], 0.0, 0.0);
    for (i, y) in ys.iter_mut().enumerate() {
        let x = xs[i];
        *y = (*y as f64 + 0.05 * x + 0.2 * jitter(i as u32)) as f32;
    }
    let opts = Some(BoundariesOptions {
        epsilon: 1e-5,
        window_size: 11,
    });

    let d = data_xy(xs.clone(), ys);
    let b = get_boundaries(&d, mu, opts);
    assert!(b.to.index.is_some());
    assert!(xs[b.to.index.unwrap()] > mu);
}

#[test]
fn right_boundary_of_left_peak_in_valley() {
    let n = 601;
    let xs = linspace(0.0, 12.0, n);
    let mu1 = 5.0;
    let mu2 = 7.0;
    let mut y1 = gaussian_mixture_f32(&xs, &[(mu1, 0.35, 100.0)], 0.0, 0.0);
    let y2 = gaussian_mixture_f32(&xs, &[(mu2, 0.35, 70.0)], 0.0, 0.0);
    for i in 0..n {
        y1[i] += y2[i];
    }
    let d = data_xy(xs.clone(), y1);

    let b = get_boundaries(&d, mu1, None);
    assert!(b.to.index.is_some());
    let x_to = xs[b.to.index.unwrap()];
    assert!(x_to > mu1 && x_to < mu2);
}

#[test]
fn quadratic_x_spacing_no_panic() {
    let n = 401;
    let xs: Vec<f64> = (0..n).map(|i| (i as f64).powi(2) / (n as f64)).collect();
    let mu = xs[n / 2];
    let ys = gaussian_mixture_f32(
        &xs,
        &[(mu, (xs.last().unwrap() - xs[0]) / 20.0, 100.0)],
        1.0,
        0.0,
    );
    let d = data_xy(xs, ys);
    let _ = get_boundaries(&d, mu, None);
}

#[test]
fn constant_series_skips_peak_pair() {
    let n = 21;
    let xs = linspace(0.0, 10.0, n);
    let ys = vec![7.0f32; n];
    let d = data_xy(xs.clone(), ys);
    let peak_idx = n / 2;
    let peak_x = xs[peak_idx];

    let b = get_boundaries(&d, peak_x, None);
    assert_eq!(b.from.index, Some(peak_idx - 1));
    assert_eq!(b.to.index, Some(peak_idx + 1));
}
