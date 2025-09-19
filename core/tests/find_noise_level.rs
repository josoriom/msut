use msut::utilities::find_noise_level::find_noise_level;

mod helpers;
use helpers::{approx_eq, gaussian_mixture_f32, make_grid, shuffle_with_seed, uniform_vec_f32};

// --- Basic validity & degeneracy ---

// #[test]
// fn empty_input_returns_infinite() {
//     let noise = find_noise_level(&[]);
//     assert!(noise.is_infinite());
// }

#[test]
fn all_zero_is_noise_noise() {
    let data = vec![0.0, 0.0, 0.0, 0.0];
    let noise = find_noise_level(&data);
    assert_eq!(noise, 0.0);
}

// --- Noise-only / unimodal ---

#[test]
fn noise_only_unimodal_has_finite_positive_noise() {
    let data = uniform_vec_f32(50_000, 5.0, 8.0, 12345);
    let noise = find_noise_level(&data);
    assert!(noise.is_finite() && noise > 0.0);
}

// --- Clear bimodal separation → confident split between clusters ---

#[test]
fn clear_bimodal_yields_noise_between_clusters() {
    let mut lows = uniform_vec_f32(40_000, 3.0, 4.0, 7);
    let highs = uniform_vec_f32(2_000, 120.0, 200.0, 13);
    lows.extend(highs);

    let noise = find_noise_level(&lows);
    assert!(noise < 200.0, "noise={} too large", noise);
}

// --- Pre-noised (zeros + only highs) → finite noise well above baseline ---

#[test]
fn prenoised_sparse_highs_has_fallback_with_high_noise() {
    let mut data = vec![0.0f32; 80_000];
    data.extend(uniform_vec_f32(2_000, 900.0, 1200.0, 29));

    let noise = find_noise_level(&data);
    assert!(noise.is_finite());
    assert!(noise > 500.0);
}

// --- Robustness to bad values ---

#[test]
fn ignores_nan_inf_and_negatives() {
    let mut data = vec![f32::NAN, f32::INFINITY, f32::NEG_INFINITY, -5.0, -0.1, 0.0];
    data.extend(uniform_vec_f32(20_000, 2.5, 3.5, 31));

    let noise = find_noise_level(&data);
    // println!("-=--::>{}", noise);
    assert!(noise > 2.5 && noise < 120.0, "noise={}", noise);
}

// --- Invariance to order ---

#[test]
#[ignore = "redo"]
fn noise_invariant_to_permutation() {
    let mut v1 = uniform_vec_f32(50_000, 3.0, 4.0, 41);
    v1.extend(uniform_vec_f32(2_000, 160.0, 240.0, 43));

    let mut v2 = v1.clone();
    shuffle_with_seed(&mut v2, 999);

    let a = find_noise_level(&v1);
    let b = find_noise_level(&v2);

    assert!(approx_eq(a as f64, b as f64, 1e-3_f64), "a={} b={}", a, b);
}

// --- Chromatogram-like (Gaussian over baseline) ---

#[test]
fn gaussian_peak_over_baseline_produces_finite_noise() {
    let xs = make_grid(0.0, 10.0, 10_001);
    let ys = gaussian_mixture_f32(&xs, &[(5.0, 0.08, 1.0)], 0.02, 0.0);
    let noise = find_noise_level(&ys);
    assert!(noise.is_finite());
}

// --- Large-N smoke test ---

#[test]
fn large_input_smoke_test_finishes_and_looks_reasonable() {
    let mut v = uniform_vec_f32(1_000_000, 2.0, 3.0, 71);
    v.extend(uniform_vec_f32(50_000, 120.0, 200.0, 73));
    let noise = find_noise_level(&v);
    assert!(noise.is_finite());
    assert!(noise > 2.0 && noise < 200.0, "noise={}", noise);
}
