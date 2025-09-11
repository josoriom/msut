use msut::utilities::scan_for_peaks::{ScanPeaksOptions, scan_for_peaks};
use msut::utilities::structs::DataXY;

mod helpers;
use helpers::{approx_eq, gaussian_mixture_f32, make_grid};

// Test: One wide LC-like Gaussian is detected at the right location.
// Why: Checks default settings for common LC peaks.
#[test]
fn single_gaussian_peak_lc_style_detected() {
    let xs = make_grid(0.0, 10.0, 2001);
    let ys = gaussian_mixture_f32(&xs, &[(5.0, 0.45, 100.0)], 0.0, 0.0);
    let data = DataXY { x: xs, y: ys };

    let peaks = scan_for_peaks(&data, None);
    assert_eq!(peaks.len(), 1);
    assert!(approx_eq(peaks[0], 5.0, 0.01));
}

// Test: One narrow Gaussian is detected precisely.
// Why: Ensures narrow peaks survive smoothing and diff steps.
#[test]
fn single_gaussian_peak_gc_style_detected() {
    let xs = make_grid(0.0, 10.0, 4001);
    let ys = gaussian_mixture_f32(&xs, &[(3.2, 0.08, 80.0)], 0.0, 0.0);
    let data = DataXY { x: xs, y: ys };

    let peaks = scan_for_peaks(&data, None);
    assert_eq!(peaks.len(), 1);
    assert!(approx_eq(peaks[0], 3.2, 0.005));
}

// Test: Two well-separated peaks are both found.
// Why: Basic multi-peak detection must work reliably.
#[test]
fn two_well_separated_peaks_detected() {
    let xs = make_grid(0.0, 10.0, 3001);
    let ys = gaussian_mixture_f32(&xs, &[(3.0, 0.25, 60.0), (7.0, 0.35, 90.0)], 0.0, 0.0);
    let data = DataXY { x: xs, y: ys };

    let peaks = scan_for_peaks(&data, None);
    assert_eq!(peaks.len(), 2);
    assert!(approx_eq(peaks[0], 3.0, 0.01));
    assert!(approx_eq(peaks[1], 7.0, 0.01));
}

// Test: Overlapping peaks are still split correctly.
// Why: Ensures the algorithm honors valleys between peaks.
#[test]
fn overlapping_peaks_still_split() {
    let xs = make_grid(0.0, 10.0, 4001);
    let ys = gaussian_mixture_f32(&xs, &[(4.4, 0.1, 80.0), (4.8, 0.1, 70.0)], 0.0, 0.0);
    let data = DataXY { x: xs, y: ys };

    let peaks = scan_for_peaks(
        &data,
        Some(ScanPeaksOptions {
            window_size: 5,
            ..Default::default()
        }),
    );
    assert_eq!(peaks.len(), 2);
    assert!(approx_eq(peaks[0], 4.4, 0.01));
    assert!(approx_eq(peaks[1], 4.8, 0.01));
}

// Test: Works when sampling is slightly non-uniform.
// Why: Real data often has small spacing variations.
#[test]
fn non_uniform_sampling_works() {
    let mut xs = make_grid(0.0, 10.0, 2501);
    for (i, x) in xs.iter_mut().enumerate() {
        *x += 0.002 * ((i as f64) * 0.003).sin();
    }
    let ys = gaussian_mixture_f32(&xs, &[(2.0, 0.18, 40.0), (8.0, 0.25, 70.0)], 0.0, 0.0);
    let data = DataXY { x: xs, y: ys };

    let peaks = scan_for_peaks(&data, None);
    assert_eq!(peaks.len(), 2);
    assert!(approx_eq(peaks[0], 2.0, 0.02));
    assert!(approx_eq(peaks[1], 8.0, 0.02));
}

// Test: Flat signal yields no peaks.
// Why: Baseline-only data should not produce detections.
#[test]
fn no_false_peaks_on_flat_signal() {
    let xs = make_grid(0.0, 10.0, 1001);
    let ys = vec![0.0; xs.len()];
    let data = DataXY { x: xs, y: ys };

    let peaks = scan_for_peaks(&data, None);
    assert!(peaks.is_empty());
}

// Test: Tuned SG window handles very narrow GC peaks.
// Why: Shows configuration can preserve tiny features.
#[test]
fn tuning_window_size_for_gc_narrow_peaks() {
    let xs = make_grid(0.0, 3.0, 6001);
    let ys = gaussian_mixture_f32(&xs, &[(1.5, 0.01, 50.0)], 0.0, 0.0);
    let data = DataXY { x: xs, y: ys };

    let opts = ScanPeaksOptions {
        window_size: 9,
        ..Default::default()
    };
    let peaks = scan_for_peaks(&data, Some(opts));
    assert_eq!(peaks.len(), 1);
    assert!(approx_eq(peaks[0], 1.5, 0.002));
}

#[test]
fn flat_baseline_with_threshold_nonzero() {
    let xs = make_grid(0.0, 10.0, 1001);
    let ys = vec![0.0f32; xs.len()];
    let data = DataXY { x: xs, y: ys };
    let opts = ScanPeaksOptions {
        ..Default::default()
    };
    let peaks = scan_for_peaks(&data, Some(opts));
    println!("---::>>{}", serde_json::to_string(&peaks).unwrap());
    assert!(peaks.is_empty());
}

#[test]
#[ignore]
fn single_peak_with_small_added_noise_and_high_threshold() {
    let xs = make_grid(0.0, 10.0, 3001);
    let ys = gaussian_mixture_f32(&xs, &[(5.0, 0.3, 20.0)], 0.0, 0.3);
    let data = DataXY {
        x: xs,
        y: ys.clone(),
    };
    let opts = ScanPeaksOptions {
        ..Default::default()
    };
    let peaks = scan_for_peaks(&data, Some(opts));
    assert_eq!(peaks.len(), 1);
}

#[test]
fn very_close_peaks_should_merge() {
    let xs = make_grid(0.0, 4.0, 4001);
    let ys = gaussian_mixture_f32(&xs, &[(2.0, 0.20, 60.0), (2.12, 0.20, 58.0)], 0.0, 0.0);
    let data = DataXY { x: xs, y: ys };
    let peaks = scan_for_peaks(&data, None);
    assert_eq!(peaks.len(), 1);
}

#[test]
fn endpoint_peak_detected() {
    let xs = make_grid(0.0, 10.0, 2001);
    let ys = gaussian_mixture_f32(&xs, &[(0.3, 0.12, 50.0)], 0.0, 0.0);
    let data = DataXY { x: xs, y: ys };
    let peaks = scan_for_peaks(&data, None);
    assert_eq!(peaks.len(), 1);
    assert!(approx_eq(peaks[0], 0.3, 0.02));
}
