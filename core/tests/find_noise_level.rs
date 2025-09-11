use msut::utilities::find_noise_level::{find_noise_level, find_noise_level_with_options};
use msut::utilities::scan_for_peaks::ScanPeaksOptions;

mod helpers;
use helpers::{approx_eq, gaussian_mixture_f32, make_grid, shuffle_with_seed, uniform_vec_f32};

// --- Basic validity & degeneracy ---

#[test]
fn empty_input_returns_infinite() {
    let noise = find_noise_level(&[]);
    assert!(noise.is_infinite());
}

#[test]
fn all_zero_or_nonpositive_is_infinite_noise() {
    let data = vec![0.0, -1.0, -5.0, 0.0];
    let noise = find_noise_level(&data);
    assert!(noise.is_infinite());
}

#[test]
fn identical_positive_values_are_infinite_noise() {
    let data = vec![123.0f32; 10_000];
    let noise = find_noise_level(&data);
    assert!(noise.is_infinite());
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

// --- Impulse-like noise: should peg near noise ceiling ---

#[test]
fn impulse_like_noise_is_high_noise() {
    let expected_noise = 150.0;
    let mut base = uniform_vec_f32(9_700, 0.0, 2.0, 880);
    let spikes = uniform_vec_f32(300, 80.0, expected_noise, 881);
    base.extend(spikes);
    shuffle_with_seed(&mut base, 13579);
    let noise = find_noise_level(&base);
    assert!(noise.is_finite(), "noise must be finite");
    assert!(noise < expected_noise, "noise={} too large", noise);
}

// --- Overlapped / shallow valley → noise near noise band ---

#[test]
fn shallow_overlap_leads_to_noise_near_noise_band() {
    let noise_band = 25.0;
    let mut noise = uniform_vec_f32(40_000, 10.0, 25.0, 19);
    let analyte = uniform_vec_f32(4_000, 22.0, noise_band, 23);
    noise.extend(analyte);

    let n = find_noise_level(&noise);
    assert!(n.is_finite());
    assert!(n < noise_band, "noise={} too large", n);
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
fn noise_invariant_to_permutation() {
    let mut v1 = uniform_vec_f32(50_000, 3.0, 4.0, 41);
    v1.extend(uniform_vec_f32(2_000, 160.0, 240.0, 43));

    let mut v2 = v1.clone();
    shuffle_with_seed(&mut v2, 999);

    let a = find_noise_level(&v1);
    let b = find_noise_level(&v2);

    assert!(approx_eq(a as f64, b as f64, 1e-3_f64), "a={} b={}", a, b);
}

// --- SGG window clamping (via FindPeaksOptions) ---

#[test]
fn sgg_window_is_clamped_to_bins_and_no_panic() {
    let mut v = uniform_vec_f32(8_000, 3.0, 4.0, 55);
    v.extend(uniform_vec_f32(500, 300.0, 600.0, 56));

    let opts = ScanPeaksOptions {
        epsilon: 1e-5,
        window_size: 100_000,
    };
    let noise = find_noise_level_with_options(&v, Some(opts));
    assert!(noise.is_finite());
}

// --- Multimodal: choose valley between noise & dominant analyte mode ---

#[test]
fn multimodal_three_clusters_has_reasonable_noise() {
    let mut v = uniform_vec_f32(30_000, 3.0, 4.0, 61);
    v.extend(uniform_vec_f32(3_000, 200.0, 300.0, 63));
    v.extend(uniform_vec_f32(1_00, 1000.0, 1200.0, 62));

    let noise = find_noise_level(&v);
    assert!(noise > 4.0 && noise < 500.0, "noise={}", noise);
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
    assert!(noise > 2.0 && noise < 120.0, "noise={}", noise);
}

#[test]
fn chromatogram_with_tall_spikes_keeps_noise_near_baseline() {
    let mut v = uniform_vec_f32(60_000, 150.0, 200.0, 7001);
    v.extend(uniform_vec_f32(300, 1400.0, 2200.0, 7002));
    shuffle_with_seed(&mut v, 7003);

    let noise = find_noise_level(&v);
    assert!(noise.is_finite(), "noise must be finite");
    assert!(noise > 150.0 && noise < 250.0, "noise={}", noise);
}

#[test]
fn chromatogram_baseline_around_200_with_central_spikes_caps_near_baseline() {
    let mut baseline = uniform_vec_f32(60_000, 150.0, 200.0, 9001);
    let hump = uniform_vec_f32(5_000, 260.0, 340.0, 9002);
    baseline.extend(hump);
    baseline.extend(uniform_vec_f32(200, 900.0, 1500.0, 9003));
    shuffle_with_seed(&mut baseline, 9004);

    let noise = find_noise_level(&baseline);
    println!("--::>>{}", noise);
    assert!(noise.is_finite(), "noise must be finite");
    assert!(
        noise > 150.0 && noise < 260.0,
        "noise={} not near baseline",
        noise
    );
}
