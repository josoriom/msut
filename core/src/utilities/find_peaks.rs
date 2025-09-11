use crate::utilities::find_noise_level::find_noise_level;
use crate::utilities::get_boundaries::{BoundariesOptions, get_boundaries};
use crate::utilities::scan_for_peaks::{ScanPeaksOptions, scan_for_peaks};
use crate::utilities::structs::{DataXY, Peak};

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
            integral_threshold: Some(0.004),
            width_threshold: Some(10),
            noise: None,
            auto_noise: Some(false),
            allow_overlap: Some(false),
            sn_ratio: Some(2),
            intensity_threshold: None,
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
    if auto_noise && filter_opts.noise.is_some() {
        panic!("Filter peaks: `auto_noise=true` is incompatible with a defined `noise`.");
    }
    let resolved_noise: f64 = if let Some(n) = filter_opts.noise {
        if auto_noise {
            panic!("Filter peaks: `auto_noise=true` is incompatible with a defined `noise`.");
        }
        n.max(0.0)
    } else if auto_noise {
        let n = find_noise_level(&data.y) as f64;
        if n.is_finite() && n > 0.0 { n } else { 0.0 }
    } else {
        0.0
    };

    let positions = scan_for_peaks(data, o.scan_peaks_options);

    let mut candidates = Vec::with_capacity(positions.len());
    let mut sum_integrals = 0.0f64;
    let use_ratio = filter_opts.integral_threshold.is_some();

    for &rt in positions.iter() {
        let b = get_boundaries(data, rt, o.get_boundaries_options);
        let (fi, fv, ti, tv) = match (b.from.index, b.from.value, b.to.index, b.to.value) {
            (Some(a), Some(av), Some(bi), Some(bv)) if a < bi => (a, av, bi, bv),
            _ => {
                let idx = closest_index(&data.x, rt);
                let (lf, rti) = fallback_bounds(&data.y, idx);
                if lf >= rti {
                    continue;
                }
                (lf, data.x[lf], rti, data.x[rti])
            }
        };

        let (integral, intensity) = xy_integration(&data.x[fi..=ti], &data.y[fi..=ti]);
        if use_ratio {
            sum_integrals += integral;
        }
        candidates.push(PeakCandidate {
            from: fv,
            to: tv,
            rt,
            integral,
            intensity,
            number_of_points: ti - fi + 1,
            ratio: integral / sum_integrals,
        });
    }

    let mut filtered = filter_peak_candidates(
        candidates,
        filter_opts,
        if use_ratio { sum_integrals } else { 0.0 },
        resolved_noise,
    );

    if !filter_opts.allow_overlap.unwrap_or(false) {
        filtered = remove_nested_peaks(filtered);
    }
    filtered
}

fn filter_peak_candidates(
    peaks: Vec<PeakCandidate>,
    opt: FilterPeaksOptions,
    sum: f64,
    noise_threshold: f64,
) -> Vec<Peak> {
    let mut out: Vec<Peak> = Vec::with_capacity(peaks.len());
    let eff_noise = if noise_threshold.is_finite() && noise_threshold > 0.0 {
        noise_threshold
    } else {
        0.0
    };
    for mut p in peaks {
        let mut pass = true;

        if eff_noise > 0.0 && p.intensity < 2.0 * eff_noise {
            pass = false;
        }

        if let Some(th) = opt.integral_threshold {
            let ratio = if sum > 0.0 { p.integral / sum } else { 0.0 };
            if ratio < th {
                pass = false;
            } else {
                p.ratio = ratio;
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

fn remove_nested_peaks(peaks: Vec<Peak>) -> Vec<Peak> {
    if peaks.len() <= 1 {
        return peaks;
    }
    let mut idxs: Vec<usize> = (0..peaks.len()).collect();
    idxs.sort_by(|&i, &j| {
        let wi = peaks[i].to - peaks[i].from;
        let wj = peaks[j].to - peaks[j].from;
        wj.partial_cmp(&wi)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                peaks[j]
                    .intensity
                    .partial_cmp(&peaks[i].intensity)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                peaks[j]
                    .integral
                    .partial_cmp(&peaks[i].integral)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    let mut keep = vec![false; peaks.len()];
    let eps = 1e-12;
    let mut kept: Vec<usize> = Vec::new();
    for i in idxs {
        let p = &peaks[i];
        let mut contained = false;
        for &k in &kept {
            let q = &peaks[k];
            let inside_start = p.from >= q.from - eps;
            let inside_end = p.to <= q.to + eps;
            let strictly_inside = (p.from > q.from + eps) || (p.to < q.to - eps);
            if inside_start && inside_end && strictly_inside {
                contained = true;
                break;
            }
        }
        if !contained {
            keep[i] = true;
            kept.push(i);
        }
    }
    let mut out: Vec<Peak> = Vec::with_capacity(kept.len());
    for (i, k) in keep.iter().enumerate() {
        if *k {
            out.push(peaks[i].clone());
        }
    }
    out.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));
    out
}

pub fn xy_integration(x: &[f64], y: &[f32]) -> (f64, f64) {
    let n = x.len();
    if n == 0 || n != y.len() {
        return (0.0, f64::NEG_INFINITY);
    }
    if n == 1 {
        return (0.0, y[0] as f64);
    }
    let mut s = 0.0f64;
    let mut m = y[0];
    for i in 0..(n - 1) {
        let dx = x[i + 1] - x[i];
        let yi = y[i] as f64;
        let yj = y[i + 1] as f64;
        s += dx * (yi + yj) * 0.5;
        if y[i + 1] > m {
            m = y[i + 1];
        }
    }
    (s, m as f64)
}

fn closest_index(xs: &[f64], v: f64) -> usize {
    if xs.is_empty() {
        return 0;
    }
    let mut lo = 0usize;
    let mut hi = xs.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if xs[mid] < v {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo == 0 {
        0
    } else if lo >= xs.len() {
        xs.len() - 1
    } else if (v - xs[lo - 1]).abs() <= (xs[lo] - v).abs() {
        lo - 1
    } else {
        lo
    }
}

fn fallback_bounds(y: &[f32], idx: usize) -> (usize, usize) {
    let n = y.len();
    if n == 0 {
        return (0, 0);
    }
    let base = y[idx].abs().max(1.0_f32);
    let tol = ((0.02_f32 * base).max(3.0_f32)).min(250.0_f32);
    let mut l = idx;
    while l > 0 && y[l - 1] + tol <= y[l] {
        l -= 1;
    }
    let mut r = idx;
    while r + 1 < n && y[r + 1] + tol <= y[r] {
        r += 1;
    }
    if l >= r {
        if idx > 0 && idx + 1 < n {
            (idx - 1, idx + 1)
        } else if idx + 1 < n {
            (idx, idx + 1)
        } else if idx > 0 {
            (idx - 1, idx)
        } else {
            (0, 0)
        }
    } else {
        (l, r)
    }
}
