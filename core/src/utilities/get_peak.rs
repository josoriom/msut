use std::panic::{AssertUnwindSafe, catch_unwind};

use crate::utilities::find_peaks::FilterPeaksOptions;
use crate::utilities::find_peaks::{FindPeaksOptions, find_peaks};
use crate::utilities::scan_for_peaks::ScanPeaksOptions;
use crate::utilities::structs::{DataXY, Peak, Roi};

const DEFAULT_WINDOW_SIZES: &[usize] = &[5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33];

fn series_len(data: &DataXY) -> usize {
    // adjust if your DataXY field names differ
    data.x.len().min(data.y.len())
}

fn odd_at_most(ws: usize, n: usize) -> Option<usize> {
    if n < 3 {
        return None;
    }
    if ws > n {
        return None;
    }
    let mut w = ws;
    if w % 2 == 0 {
        w = w.saturating_sub(1);
    }
    if w >= 3 { Some(w) } else { None }
}

pub fn get_peak_across_windows(
    data: &DataXY,
    roi: Roi,
    options: Option<FindPeaksOptions>,
    window_sizes: Option<&[usize]>,
) -> Option<Peak> {
    let base = options.unwrap_or_default();
    let wss = window_sizes.unwrap_or(DEFAULT_WINDOW_SIZES);
    let n = series_len(data);
    if n < 3 {
        return None;
    }

    let mut candidates: Vec<Peak> = Vec::new();

    for &ws in wss {
        let Some(ws_eff) = odd_at_most(ws, n) else {
            continue;
        };

        let base_scan = base
            .scan_peaks_options
            .as_ref()
            .map(|s| ScanPeaksOptions {
                epsilon: s.epsilon,
                window_size: s.window_size,
            })
            .unwrap_or_default();

        let base_filter = base
            .filter_peaks_options
            .as_ref()
            .map(|f| FilterPeaksOptions {
                integral_threshold: f.integral_threshold,
                width_threshold: f.width_threshold,
                intensity_threshold: f.intensity_threshold,
                noise: f.noise,
                auto_noise: f.auto_noise,
                allow_overlap: f.allow_overlap,
                sn_ratio: f.sn_ratio,
            });

        let opts = FindPeaksOptions {
            get_boundaries_options: base.get_boundaries_options,
            filter_peaks_options: base_filter,
            scan_peaks_options: Some(ScanPeaksOptions {
                epsilon: base_scan.epsilon,
                window_size: ws_eff,
            }),
        };

        let peaks = match catch_unwind(AssertUnwindSafe(|| find_peaks(data, Some(opts)))) {
            Ok(v) => v,
            Err(_) => continue,
        };
        if peaks.is_empty() {
            continue;
        }

        if let Some(best_for_ws) = closest_to_rt(&peaks, roi.rt) {
            candidates.push(best_for_ws);
        }
    }

    if candidates.is_empty() {
        return None;
    }
    best_closest_to_rt(&candidates, roi.rt)
}
pub fn get_peak(data: &DataXY, roi: Roi, options: Option<FindPeaksOptions>) -> Option<Peak> {
    get_peak_across_windows(data, roi, options, None)
}

fn closest_to_rt(list: &[Peak], rt: f64) -> Option<Peak> {
    if list.is_empty() || !rt.is_finite() {
        return None;
    }
    list.iter()
        .filter(|p| p.rt.is_finite())
        .min_by(|a, b| (a.rt - rt).abs().total_cmp(&(b.rt - rt).abs()))
        .cloned()
}

fn best_closest_to_rt(list: &[Peak], rt: f64) -> Option<Peak> {
    if list.is_empty() || !rt.is_finite() {
        return None;
    }

    let in_window: Vec<&Peak> = list
        .iter()
        .filter(|p| {
            let a = p.from.min(p.to);
            let b = p.from.max(p.to);
            p.rt.is_finite()
                && p.from.is_finite()
                && p.to.is_finite()
                && (b - a) > 0.0
                && p.rt >= a
                && p.rt <= b
        })
        .collect();

    if in_window.is_empty() {
        return None;
    }

    let min_delta = in_window
        .iter()
        .map(|p| (p.rt - rt).abs())
        .fold(f64::INFINITY, f64::min);

    const EPS: f64 = 0.1;
    let tied: Vec<&Peak> = in_window
        .into_iter()
        .filter(|p| ((p.rt - rt).abs() - min_delta).abs() <= EPS)
        .collect();

    if tied.len() == 1 {
        return tied.into_iter().next().cloned();
    }

    let rts: Vec<f64> = tied.iter().map(|p| p.rt).collect();
    let rt_median = median(&rts);
    let rt_mad = mad(&rts, rt_median);

    let rt_band: Vec<&Peak> = if rt_mad > 0.0 && rt_mad.is_finite() {
        tied.into_iter()
            .filter(|p| (p.rt - rt_median).abs() <= rt_mad)
            .collect()
    } else {
        tied
    };

    if rt_band.len() == 1 {
        return rt_band.into_iter().next().cloned();
    }

    let pool = if !rt_band.is_empty() { rt_band } else { vec![] };
    let base_pool = if pool.is_empty() {
        list.iter().collect::<Vec<_>>()
    } else {
        pool
    };

    let widths: Vec<f64> = base_pool.iter().map(|p| (p.to - p.from).abs()).collect();
    let w_mean = mean(&widths);
    let w_sd = stddev(&widths);

    let width_band: Vec<&Peak> = if w_sd > 0.0 && w_sd.is_finite() {
        base_pool
            .into_iter()
            .filter(|p| ((p.to - p.from).abs() - w_mean).abs() <= w_sd)
            .collect()
    } else {
        Vec::new()
    };

    let final_pool = if !width_band.is_empty() {
        width_band
    } else {
        list.iter().collect()
    };

    final_pool
        .into_iter()
        .max_by(|a, b| ((a.to - a.from).abs()).total_cmp(&((b.to - b.from).abs())))
        .cloned()
}

fn mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let sum: f64 = v.iter().copied().sum();
    sum / (v.len() as f64)
}

fn stddev(v: &[f64]) -> f64 {
    let n = v.len();
    if n < 2 {
        return 0.0;
    }
    let m = mean(v);
    let var = v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (n as f64);
    var.sqrt()
}

fn median(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mut s: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
    if s.is_empty() {
        return 0.0;
    }
    s.sort_by(|a, b| a.total_cmp(b));
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        (s[n / 2 - 1] + s[n / 2]) / 2.0
    }
}

fn mad(v: &[f64], med: f64) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let devs: Vec<f64> = v
        .iter()
        .copied()
        .filter(|x| x.is_finite())
        .map(|x| (x - med).abs())
        .collect();
    median(&devs)
}
