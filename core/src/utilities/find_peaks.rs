use std::cmp::Ordering;

use crate::utilities::calculate_baseline::{BaselineOptions, calculate_baseline};
use crate::utilities::closest_index;
use crate::utilities::find_noise_level::find_noise_level;
use crate::utilities::get_boundaries::{Boundaries, BoundariesOptions, get_boundaries};

use crate::utilities::scan_for_peaks::{ScanPeaksOptions, scan_for_peaks_across_windows};
use crate::utilities::structs::{DataXY, Peak};
use crate::utilities::utilities::xy_integration;

const DEFAULT_WINDOW_SIZES: &[usize] = &[5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33];

#[derive(Clone, Copy, Debug)]
pub struct FilterPeaksOptions {
    pub integral_threshold: Option<f64>,
    pub width_threshold: Option<usize>,
    pub intensity_threshold: Option<f64>,
    pub noise: Option<f64>,
    pub auto_noise: Option<bool>,
    pub auto_baseline: Option<bool>,
    pub allow_overlap: Option<bool>,
    pub sn_ratio: Option<f64>,
}

impl Default for FilterPeaksOptions {
    fn default() -> Self {
        Self {
            integral_threshold: None,
            width_threshold: Some(5),
            intensity_threshold: None,
            noise: None,
            auto_noise: Some(false),
            auto_baseline: Some(false),
            allow_overlap: Some(false),
            sn_ratio: Some(1.0),
        }
    }
}

#[derive(Clone, Debug)]
struct PeakCandidate {
    from: f64,
    to: f64,
    rt: f64,
    integral: f64,
    intensity: f64,
    number_of_points: usize,
    ratio: f64,
    noise: f64,
}

impl From<PeakCandidate> for Peak {
    fn from(c: PeakCandidate) -> Self {
        Peak {
            from: c.from,
            to: c.to,
            rt: c.rt,
            integral: c.integral,
            intensity: c.intensity,
            ratio: c.ratio,
            np: c.number_of_points,
            noise: c.noise,
        }
    }
}

#[derive(Clone, Debug)]
pub struct FindPeaksOptions {
    pub get_boundaries_options: Option<BoundariesOptions>,
    pub filter_peaks_options: Option<FilterPeaksOptions>,
    pub scan_peaks_options: Option<ScanPeaksOptions>,
    pub baseline_options: Option<BaselineOptions>,
}
impl Default for FindPeaksOptions {
    fn default() -> Self {
        Self {
            get_boundaries_options: Some(BoundariesOptions::default()),
            filter_peaks_options: Some(FilterPeaksOptions::default()),
            scan_peaks_options: Some(ScanPeaksOptions::default()),
            baseline_options: Some(BaselineOptions::default()),
        }
    }
}

pub fn find_peaks(data: &DataXY, options: Option<FindPeaksOptions>) -> Vec<Peak> {
    let o = options.unwrap_or_default();
    let filter_opts = o.filter_peaks_options.unwrap_or_default();
    let base_opts = o.baseline_options.unwrap_or_default();

    let y64: Vec<f64> = data.y.clone();
    let floor = if filter_opts.auto_baseline.unwrap_or(false) {
        let mut b = base_opts.clone();
        b.level = Some(0);
        calculate_baseline(&y64, b)
    } else {
        vec![0.0; y64.len()]
    };

    let y_center: Vec<f64> = y64
        .iter()
        .zip(&floor)
        .map(|(a, m)| (a - m).max(0.0))
        .collect();

    let auto_noise = filter_opts.auto_noise.unwrap_or(false);
    if auto_noise && filter_opts.noise.is_some() {
        panic!("auto_noise=true cannot be used with noise");
    }

    let noise = if auto_noise {
        find_noise_level(&y_center)
    } else {
        filter_opts.noise.unwrap_or(0.0).max(0.0)
    };

    let normalized_data = DataXY {
        x: data.x.clone(),
        y: y_center,
    };

    let positions = scan_for_peaks_across_windows(
        &normalized_data,
        o.scan_peaks_options,
        Some(DEFAULT_WINDOW_SIZES),
    );
    if positions.is_empty() {
        return Vec::new();
    }

    let mut bopt = o.get_boundaries_options.unwrap_or_default();
    bopt.noise = noise;
    let mut candidates: Vec<PeakCandidate> = Vec::with_capacity(positions.len());
    for seed_rt in positions {
        let b = get_boundaries(&normalized_data, seed_rt, Some(bopt));
        let seed_idx = closest_index(&normalized_data.x, seed_rt);
        let (rt, apex_y) = apex_in_window(&normalized_data, &b)
            .unwrap_or((normalized_data.x[seed_idx], normalized_data.y[seed_idx]));

        if apex_y <= noise {
            continue;
        }

        match (b.from.index, b.from.value, b.to.index, b.to.value) {
            (Some(fi), Some(fx), Some(ti), Some(tx)) if fi < ti => {
                let (integral, intensity) = xy_integration(&data.x[fi..=ti], &data.y[fi..=ti]);
                let cand = PeakCandidate {
                    from: fx,
                    to: tx,
                    rt,
                    integral,
                    intensity,
                    number_of_points: ti - fi + 1,
                    ratio: 0.0,
                    noise,
                };
                candidates.push(cand);
            }
            _ => {}
        }
    }
    if candidates.is_empty() {
        return Vec::new();
    }

    let mut max_intensity = 0.0f64;
    for c in &candidates {
        if c.intensity > max_intensity {
            max_intensity = c.intensity;
        }
    }
    let mut peaks = filter_peak_candidates(candidates, filter_opts);

    peaks.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));
    peaks = dedupe_near_identical(peaks);

    if !peaks.is_empty() {
        let mut cutoff = 0.0_f64;
        if noise > 0.0 {
            let sn_mult = filter_opts.sn_ratio.unwrap_or(1.0) as f64;
            cutoff = sn_mult * noise;
        }
        if let Some(user_int) = filter_opts.intensity_threshold {
            cutoff = cutoff.max(user_int);
        }
        if cutoff > 0.0 {
            peaks.retain(|p| p.intensity > cutoff);
        }
    }

    if peaks.len() > 1 {
        peaks = suppress_contained_peaks(data, peaks);
    }
    peaks
}

pub fn apex_in_window(data: &DataXY, b: &Boundaries) -> Option<(f64, f64)> {
    let l = b.from.index?;
    let r = b.to.index?;
    if l >= r {
        return None;
    }
    let (off, &ymax) = data.y[l..=r]
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(Ordering::Equal))?;
    let i = l + off;
    Some((data.x[i], ymax))
}

fn filter_peak_candidates(peaks: Vec<PeakCandidate>, opt: FilterPeaksOptions) -> Vec<Peak> {
    let mut out: Vec<Peak> = Vec::with_capacity(peaks.len());

    for p in peaks {
        let mut pass = true;

        if let Some(min_i) = opt.intensity_threshold {
            if p.intensity < min_i {
                pass = false;
            }
        }
        if pass {
            if let Some(wth) = opt.width_threshold {
                if p.number_of_points <= wth {
                    pass = false;
                }
            }
        }
        if pass {
            out.push(Peak::from(p));
        }
    }
    out
}

fn dedupe_near_identical(peaks: Vec<Peak>) -> Vec<Peak> {
    if peaks.len() <= 1 {
        return peaks;
    }
    let mut out: Vec<Peak> = Vec::with_capacity(peaks.len());
    let eps_rt = 1e-6;
    let eps_w = 1e-6;
    let mut i = 0usize;
    while i < peaks.len() {
        let p = peaks[i].clone();
        let mut j = i + 1;
        let mut keep_p = p.clone();
        while j < peaks.len() {
            let q = &peaks[j];
            let same = (p.from - q.from).abs() <= eps_w
                && (p.to - q.to).abs() <= eps_w
                && (p.rt - q.rt).abs() <= eps_rt;
            if same {
                if q.intensity > keep_p.intensity {
                    keep_p = q.clone();
                }
                j += 1;
            } else {
                break;
            }
        }
        out.push(keep_p);
        i = j;
    }
    out
}

fn suppress_contained_peaks(data: &DataXY, peaks: Vec<Peak>) -> Vec<Peak> {
    if peaks.len() <= 1 {
        return peaks;
    }
    let mut order: Vec<usize> = (0..peaks.len()).collect();
    order.sort_by(|&i, &j| {
        peaks[j]
            .intensity
            .partial_cmp(&peaks[i].intensity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut dx_min = f64::INFINITY;
    for w in data.x.windows(2) {
        let d = w[1] - w[0];
        if d > 0.0 && d < dx_min {
            dx_min = d;
        }
    }
    let eps = if dx_min.is_finite() {
        2.0 * dx_min
    } else {
        0.01
    };
    let mut keep = vec![true; peaks.len()];
    for (a_idx, &ia) in order.iter().enumerate() {
        if !keep[ia] {
            continue;
        }
        let la = peaks[ia].from;
        let ra = peaks[ia].to;
        for &ib in order.iter().skip(a_idx + 1) {
            if !keep[ib] {
                continue;
            }
            let lb = peaks[ib].from;
            let rb = peaks[ib].to;
            let wb = (rb - lb).abs();
            if wb <= eps {
                keep[ib] = false;
                continue;
            }
            let l = la.max(lb);
            let r = ra.min(rb);
            let overlap = (r - l).max(0.0);
            if overlap >= 0.90 * wb {
                keep[ib] = false;
            }
        }
    }
    let mut out = Vec::<Peak>::new();
    for i in 0..peaks.len() {
        if keep[i] {
            out.push(peaks[i].clone());
        }
    }
    out.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));
    out
}
