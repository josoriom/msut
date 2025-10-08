use msut::utilities::calculate_baseline::BaselineOptions;
use msut::utilities::find_peaks::{FilterPeaksOptions, FindPeaksOptions, find_peaks};
use msut::utilities::scan_for_peaks::ScanPeaksOptions;
use msut::utilities::structs::DataXY;

mod helpers;
use helpers::{approx_eq, data_xy, gaussian_mixture_f32, make_grid};

fn nonuniform_grid(start: f64, end: f64, n: usize, a: f64, b: f64) -> Vec<f64> {
    let mut xs = Vec::with_capacity(n.max(2));
    let mut x = start;
    let mut i = 0usize;
    while x <= end && xs.len() < n {
        xs.push(x);
        let step = if i % 2 == 0 { a } else { b };
        x += step;
        i += 1;
    }
    if *xs.last().unwrap() < end {
        xs.push(end);
    }
    xs
}

// empty input should return no peaks
#[test]
fn empty_input_returns_empty() {
    let data = DataXY {
        x: vec![],
        y: vec![],
    };
    let res = find_peaks(&data, None);
    assert!(res.is_empty());
}

// single point cannot define a peak
#[test]
fn single_point_returns_empty() {
    let data = DataXY {
        x: vec![0.0],
        y: vec![0.0],
    };
    let res = find_peaks(&data, None);
    assert!(res.is_empty());
}

// perfectly flat series has no peaks
#[test]
fn constant_series_no_peaks() {
    let xs = make_grid(0.0, 10.0, 101);
    let ys = vec![0.0f32; xs.len()];
    let data = data_xy(xs, ys);
    let res = find_peaks(&data, None);
    assert!(res.is_empty());
}

// one well-formed gaussian should be detected near its center
#[test]
fn single_gaussian_peak_detected() {
    let xs = make_grid(0.0, 10.0, 1001);
    let ys = gaussian_mixture_f32(&xs, &[(4.0, 0.1, 1.0)], 0.0, 0.0);
    let data = data_xy(xs.clone(), ys);
    let res = find_peaks(&data, None);
    assert_eq!(res.len(), 1);
    let p = &res[0];
    assert!(p.from < p.to);
    assert!(p.intensity > 0.5);
    assert!(approx_eq(p.rt, 4.0, 0.06));
}

// two far-apart peaks should both be found and ordered by time
#[test]
fn two_separated_peaks_detected() {
    let xs = make_grid(0.0, 10.0, 2001);
    let ys = gaussian_mixture_f32(&xs, &[(3.0, 0.09, 1.0), (7.0, 0.09, 0.8)], 0.0, 0.0);
    let data = data_xy(xs, ys);
    let res = find_peaks(&data, None);
    assert_eq!(res.len(), 2);
    assert!(res[0].rt < res[1].rt);
    assert!(approx_eq(res[0].rt, 3.0, 0.06));
    assert!(approx_eq(res[1].rt, 7.0, 0.06));
}

// very close peaks should merge into a single detection due to min separation
#[test]
fn close_peaks_merge_into_one() {
    let xs = make_grid(0.0, 2.0, 801);
    let ys = gaussian_mixture_f32(
        &xs,
        &[(1.0 - 0.01, 0.03, 1.0), (1.0 + 0.01, 0.03, 0.9)],
        0.0,
        0.0,
    );
    let data = data_xy(xs, ys);
    let res = find_peaks(&data, None);
    assert_eq!(res.len(), 1);
}

// non-uniform sampling should still detect the correct number of peaks
#[test]
fn non_uniform_sampling_works() {
    let xs = nonuniform_grid(0.0, 10.0, 900, 0.011, 0.019);
    let ys = gaussian_mixture_f32(&xs, &[(2.5, 0.12, 1.0), (6.5, 0.15, 0.8)], 0.0, 0.0);
    let data = data_xy(xs, ys);
    let res = find_peaks(&data, None);
    assert_eq!(res.len(), 2);
}

// thresholded mode should detect a peak only where the smoothed signal exceeds threshold
#[test]
fn threshold_mode_detects_peak() {
    let xs = make_grid(0.0, 10.0, 1201);
    let ys = gaussian_mixture_f32(&xs, &[(5.0, 0.12, 1.0)], 0.2, 0.0);
    let data = data_xy(xs, ys);
    let find_opts = ScanPeaksOptions {
        epsilon: 1e-5,
        window_size: 15,
        ..Default::default()
    };
    let opts = FindPeaksOptions {
        get_boundaries_options: None,
        filter_peaks_options: Some(FilterPeaksOptions::default()),
        scan_peaks_options: Some(find_opts),
        baseline_options: Some(BaselineOptions {
            ..Default::default()
        }),
    };
    let res = find_peaks(&data, Some(opts));
    assert_eq!(res.len(), 1);
    assert!(approx_eq(res[0].rt, 5.0, 0.07));
}

// baseline offset should not prevent detection and intensity should reflect base+amp
#[test]
fn baseline_offset_ok() {
    let xs = make_grid(0.0, 8.0, 1601);
    let ys = gaussian_mixture_f32(&xs, &[(4.0, 0.12, 0.8)], 0.5, 0.0);
    let data = data_xy(xs, ys);
    let res = find_peaks(&data, None);
    assert_eq!(res.len(), 1);
    assert!(res[0].intensity > 0.8);
}

// overlapping but distinct peaks should still be split when separation exceeds minimum
#[test]
fn overlapping_peaks_still_split() {
    let xs = make_grid(0.0, 2.0, 1601);
    let ys = gaussian_mixture_f32(&xs, &[(0.9, 0.08, 1.0), (1.1, 0.08, 0.95)], 0.0, 0.0);
    let data = data_xy(xs, ys);
    let res = find_peaks(&data, None);
    assert_eq!(res.len(), 2);
    assert!(res[0].rt < res[1].rt);
}

#[test]
fn chromatogram_detects_tall_narrow_peak_on_200_baseline() {
    use msut::utilities::structs::DataXY;

    let xs = make_grid(0.0, 5.0, 6000);
    let ys = gaussian_mixture_f32(&xs, &[(3.45, 0.007, 480.0)], 200.0, 6.0);
    let data = DataXY { x: xs, y: ys };

    let peaks = find_peaks(&data, Some(FindPeaksOptions::default()));
    assert!(
        peaks.iter().any(|p| (p.rt - 3.45).abs() < 0.02),
        "missing tall central peak near rt=3.45"
    );
}

#[test]
fn noise_defined_filters_below_threshold() {
    let xs = make_grid(0.0, 10.0, 2001);
    let ys = gaussian_mixture_f32(&xs, &[(5.0, 0.08, 0.05)], 0.0, 0.0);
    let data = data_xy(xs, ys);
    let opts = FindPeaksOptions {
        get_boundaries_options: None,
        filter_peaks_options: Some(FilterPeaksOptions {
            integral_threshold: None,
            width_threshold: None,
            noise: Some(0.06),
            auto_noise: Some(false),
            allow_overlap: Some(true),
            ..Default::default()
        }),
        scan_peaks_options: None,
        baseline_options: Some(BaselineOptions {
            ..Default::default()
        }),
    };
    let res = find_peaks(&data, Some(opts));
    assert!(res.is_empty());
}

#[test]
fn noise_defined_keeps_just_above_threshold() {
    let xs = make_grid(0.0, 10.0, 2001);
    let ys = gaussian_mixture_f32(&xs, &[(5.0, 0.08, 0.061)], 0.0, 0.0);
    let data = data_xy(xs.clone(), ys);
    let opts = FindPeaksOptions {
        get_boundaries_options: None,
        filter_peaks_options: Some(FilterPeaksOptions {
            integral_threshold: None,
            width_threshold: None,
            noise: Some(0.04),
            auto_noise: Some(false),
            allow_overlap: Some(true),
            ..Default::default()
        }),
        scan_peaks_options: None,
        baseline_options: Some(BaselineOptions {
            ..Default::default()
        }),
    };
    let res = find_peaks(&data, Some(opts));
    assert_eq!(res.len(), 1);
    assert!(approx_eq(res[0].rt, 5.0, 0.08));
}

#[test]
fn negative_noise_is_clamped_to_zero() {
    let xs = make_grid(0.0, 6.0, 1201);
    let ys = gaussian_mixture_f32(&xs, &[(3.0, 0.06, 0.2)], 0.0, 0.0);
    let data = data_xy(xs.clone(), ys);
    let opts = FindPeaksOptions {
        get_boundaries_options: None,
        filter_peaks_options: Some(FilterPeaksOptions {
            integral_threshold: None,
            width_threshold: None,
            noise: Some(-5.0),
            auto_noise: Some(false),
            allow_overlap: Some(true),
            ..Default::default()
        }),
        scan_peaks_options: None,
        baseline_options: Some(BaselineOptions {
            ..Default::default()
        }),
    };
    let res = find_peaks(&data, Some(opts));
    assert_eq!(res.len(), 1);
    assert!(approx_eq(res[0].rt, 3.0, 0.07));
}

#[test]
#[should_panic]
fn auto_noise_with_defined_noise_panics() {
    let xs = make_grid(0.0, 4.0, 801);
    let ys = gaussian_mixture_f32(&xs, &[(2.0, 0.05, 1.0)], 0.0, 0.0);
    let data = data_xy(xs, ys);
    let opts = FindPeaksOptions {
        get_boundaries_options: None,
        filter_peaks_options: Some(FilterPeaksOptions {
            integral_threshold: None,
            width_threshold: None,
            noise: Some(0.2),
            auto_noise: Some(true),
            allow_overlap: Some(true),
            ..Default::default()
        }),
        scan_peaks_options: None,
        baseline_options: Some(BaselineOptions {
            ..Default::default()
        }),
    };
    let _ = find_peaks(&data, Some(opts));
}

#[test]
fn extremely_high_noise_removes_all_peaks() {
    let xs = make_grid(0.0, 6.0, 1201);
    let ys = gaussian_mixture_f32(&xs, &[(3.0, 0.08, 1000.0)], 0.0, 0.0);
    let data = data_xy(xs, ys);
    let opts = FindPeaksOptions {
        get_boundaries_options: None,
        filter_peaks_options: Some(FilterPeaksOptions {
            integral_threshold: None,
            width_threshold: None,
            noise: Some(1e9),
            auto_noise: Some(false),
            allow_overlap: Some(true),
            ..Default::default()
        }),
        scan_peaks_options: None,
        baseline_options: Some(BaselineOptions {
            ..Default::default()
        }),
    };
    let res = find_peaks(&data, Some(opts));
    assert!(res.is_empty());
}

#[test]
fn edge_peak_removed_by_noise() {
    let xs = make_grid(0.0, 2.0, 801);
    let ys = gaussian_mixture_f32(&xs, &[(0.2, 0.05, 0.5)], 0.0, 0.0);
    let data = data_xy(xs, ys);
    let opts = FindPeaksOptions {
        get_boundaries_options: None,
        filter_peaks_options: Some(FilterPeaksOptions {
            integral_threshold: None,
            width_threshold: None,
            noise: Some(0.6),
            auto_noise: Some(false),
            allow_overlap: Some(true),
            ..Default::default()
        }),
        scan_peaks_options: None,
        baseline_options: Some(BaselineOptions {
            ..Default::default()
        }),
    };
    let res = find_peaks(&data, Some(opts));
    assert!(res.is_empty());
}

#[test]
fn non_uniform_sampling_with_noise_filters_low_peak_keeps_high_peak() {
    let xs = nonuniform_grid(0.0, 10.0, 900, 0.011, 0.019);
    let ys = gaussian_mixture_f32(&xs, &[(2.5, 0.12, 0.20), (6.5, 0.15, 1.00)], 0.0, 0.0);
    let data = data_xy(xs.clone(), ys);
    let opts = FindPeaksOptions {
        get_boundaries_options: None,
        filter_peaks_options: Some(FilterPeaksOptions {
            integral_threshold: None,
            width_threshold: None,
            noise: Some(0.25),
            auto_noise: Some(false),
            allow_overlap: Some(true),
            ..Default::default()
        }),
        scan_peaks_options: None,
        baseline_options: Some(BaselineOptions {
            ..Default::default()
        }),
    };
    let res = find_peaks(&data, Some(opts));
    assert_eq!(res.len(), 1);
    assert!(approx_eq(res[0].rt, 6.5, 0.08));
}

#[test]
fn default_noise_zero_when_disabled_and_undefined_detects_peak() {
    let xs = make_grid(0.0, 6.0, 1201);
    let ys = gaussian_mixture_f32(&xs, &[(3.0, 0.07, 0.03)], 0.0, 0.0);
    let data = data_xy(xs.clone(), ys);
    let opts = FindPeaksOptions {
        get_boundaries_options: None,
        filter_peaks_options: Some(FilterPeaksOptions {
            integral_threshold: Some(0.004),
            width_threshold: None,
            noise: None,
            auto_noise: Some(false),
            allow_overlap: Some(true),
            ..Default::default()
        }),
        scan_peaks_options: None,
        baseline_options: Some(BaselineOptions {
            ..Default::default()
        }),
    };
    let res = find_peaks(&data, Some(opts));
    assert_eq!(res.len(), 1);
    assert!(approx_eq(res[0].rt, 3.0, 0.08));
}
