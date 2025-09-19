use crate::utilities::closest_index;
use crate::utilities::find_noise_level::find_noise_level;
use crate::utilities::get_boundaries::{BoundariesOptions, get_boundaries};
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
    pub allow_overlap: Option<bool>,
    pub sn_ratio: Option<usize>,
}
impl Default for FilterPeaksOptions {
    fn default() -> Self {
        Self {
            integral_threshold: None,
            width_threshold: Some(5),
            intensity_threshold: None,
            noise: None,
            auto_noise: Some(false),
            allow_overlap: Some(false),
            sn_ratio: Some(1),
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
            np: c.number_of_points as i32,
            noise: c.noise,
        }
    }
}

#[derive(Clone, Debug)]
pub struct FindPeaksOptions {
    pub get_boundaries_options: Option<BoundariesOptions>,
    pub filter_peaks_options: Option<FilterPeaksOptions>,
    pub scan_peaks_options: Option<ScanPeaksOptions>,
}
impl Default for FindPeaksOptions {
    fn default() -> Self {
        Self {
            get_boundaries_options: Some(BoundariesOptions::default()),
            filter_peaks_options: Some(FilterPeaksOptions::default()),
            scan_peaks_options: None,
        }
    }
}

pub fn find_peaks(data: &DataXY, options: Option<FindPeaksOptions>) -> Vec<Peak> {
    let o = options.unwrap_or_default();
    let filter_opts = o.filter_peaks_options.unwrap_or_default();

    let auto_noise = filter_opts.auto_noise.unwrap_or(false);
    let manual_noise_opt = filter_opts.noise;
    if auto_noise && manual_noise_opt.is_some() {
        panic!("Filter peaks: `auto_noise=true` is incompatible with a defined `noise`.");
    }
    let resolved_noise: f64 = if auto_noise {
        let n = find_noise_level(&data.y) as f64;
        if n.is_finite() && n > 0.0 { n } else { 0.0 }
    } else if let Some(n_raw) = manual_noise_opt {
        let n = if n_raw.is_finite() { n_raw } else { 0.0 };
        if n >= 0.0 { n } else { 0.0 }
    } else {
        0.0
    };

    let positions =
        scan_for_peaks_across_windows(data, o.scan_peaks_options, Some(DEFAULT_WINDOW_SIZES));
    if positions.is_empty() {
        return Vec::new();
    }

    let mut bopt = o.get_boundaries_options.unwrap_or_default();
    bopt.noise = resolved_noise;

    let mut candidates: Vec<PeakCandidate> = Vec::with_capacity(positions.len());
    let mut sum_integrals = 0.0f64;
    let use_ratio = filter_opts.integral_threshold.is_some();

    for rt in positions {
        let ai = closest_index(&data.x, rt);
        let apex_h = data.y[ai] as f64;
        if apex_h < resolved_noise {
            continue;
        }

        let b = get_boundaries(data, rt, Some(bopt));
        match (b.from.index, b.from.value, b.to.index, b.to.value) {
            (Some(fi), Some(fx), Some(ti), Some(tx)) if fi < ti => {
                let (integral, intensity) = xy_integration(&data.x[fi..=ti], &data.y[fi..=ti]);
                if use_ratio {
                    sum_integrals += integral;
                }
                candidates.push(PeakCandidate {
                    from: fx,
                    to: tx,
                    rt,
                    integral,
                    intensity,
                    number_of_points: ti - fi + 1,
                    ratio: 0.0,
                    noise: resolved_noise,
                });
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

    let mut peaks = filter_peak_candidates(
        candidates,
        filter_opts,
        if use_ratio { sum_integrals } else { 0.0 },
        max_intensity,
    );

    peaks.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));
    peaks = dedupe_near_identical(peaks);

    // if !filter_opts.allow_overlap.unwrap_or(false) && peaks.len() > 1 {
    //     peaks = prune_overlaps_by_valley(data, peaks, resolved_noise);
    // }

    if !peaks.is_empty() {
        let mut cutoff = 0.0_f64;
        if resolved_noise > 0.0 {
            let sn_mult = filter_opts.sn_ratio.unwrap_or(1) as f64;
            cutoff = sn_mult * resolved_noise;
        }
        if let Some(user_int) = filter_opts.intensity_threshold {
            cutoff = cutoff.max(user_int);
        }
        if cutoff > 0.0 {
            peaks.retain(|p| p.intensity > cutoff);
        }
    }

    if peaks.len() > 1 {
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
        let mut keep = vec![true; peaks.len()];
        for (a_idx, &ia) in order.iter().enumerate() {
            if !keep[ia] {
                continue;
            }
            let a = &peaks[ia];
            let width = (a.to - a.from).abs().max(if dx_min.is_finite() {
                2.0 * dx_min
            } else {
                0.01
            });
            for &ib in order.iter().skip(a_idx + 1) {
                if !keep[ib] {
                    continue;
                }
                let b = &peaks[ib];
                let inside_apex = b.rt >= a.from && b.rt <= a.to;
                if !inside_apex {
                    continue;
                }
                let close = (a.rt - b.rt).abs() <= 0.75 * width;
                if !close {
                    continue;
                }
                let small_rel = if a.intensity > 0.0 {
                    b.intensity / a.intensity < 0.25
                } else {
                    true
                };
                if !small_rel {
                    continue;
                }
                let iax = closest_index(&data.x, a.rt);
                let ibx = closest_index(&data.x, b.rt);
                let l = iax.min(ibx);
                let r = iax.max(ibx);
                let mut valley = f32::INFINITY;
                let mut j = l;
                while j <= r {
                    let v = data.y[j];
                    if v < valley {
                        valley = v;
                    }
                    j += 1;
                }
                let valley_low = (valley as f64) <= (1.5 * resolved_noise);
                if valley_low {
                    keep[ib] = false;
                }
            }
        }
        let mut tmp = Vec::<Peak>::new();
        for i in 0..peaks.len() {
            if keep[i] {
                tmp.push(peaks[i].clone());
            }
        }
        tmp.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));
        peaks = tmp;
    }

    if peaks.len() > 1 {
        peaks = suppress_contained_peaks(data, peaks);
    }

    peaks
}

fn filter_peak_candidates(
    peaks: Vec<PeakCandidate>,
    opt: FilterPeaksOptions,
    sum: f64,
    max_intensity: f64,
) -> Vec<Peak> {
    let mut out: Vec<Peak> = Vec::with_capacity(peaks.len());
    let tall_gate = 0.60 * max_intensity;

    for mut p in peaks {
        let mut pass = true;

        if let Some(min_i) = opt.intensity_threshold {
            if p.intensity < min_i {
                pass = false;
            }
        }
        if pass {
            if let Some(th) = opt.integral_threshold {
                let ratio = if sum > 0.0 { p.integral / sum } else { 0.0 };
                if ratio < th {
                    pass = false;
                } else {
                    p.ratio = ratio;
                }
            }
        }
        if pass {
            if let Some(wth) = opt.width_threshold {
                if p.number_of_points <= wth && p.intensity < tall_gate {
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
