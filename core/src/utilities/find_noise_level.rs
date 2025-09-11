use crate::utilities::scan_for_peaks::{ScanPeaksOptions, scan_for_peaks};
use crate::utilities::sgg::{SggOptions, sgg};
use crate::utilities::structs::DataXY;

#[inline]
fn is_finite_positive(x: f32) -> bool {
    x.is_finite() && x > 0.0
}

fn choose_bins(n: usize) -> usize {
    let mut b = 1usize;
    let mut t = (n as f64).sqrt().ceil() as usize;
    if t < 256 {
        t = 256;
    }
    while b < t {
        b <<= 1;
    }
    if b > 16384 { 16384 } else { b }
}

#[inline]
fn clamp01(x: f32) -> f32 {
    if x < 0.0 {
        0.0
    } else if x > 1.0 {
        1.0
    } else {
        x
    }
}

#[inline]
fn bin_center_log(min_log: f32, bin_w: f32, i: usize) -> f32 {
    min_log + (i as f32 + 0.5) * bin_w
}

#[inline]
fn x_to_nearest_bin(x: f64, min_log: f32, bin_w: f32, bins: usize) -> usize {
    let mut i = (((x as f32 - min_log) / bin_w) - 0.5).round() as isize;
    if i < 0 {
        i = 0;
    }
    let mut u = i as usize;
    if u >= bins {
        u = bins - 1;
    }
    u
}

fn valley_between(h: &[f32], l: usize, r: usize) -> Option<usize> {
    if r <= l + 1 {
        return None;
    }
    let mut idx = l + 1;
    let mut v = h[idx];
    for i in (l + 1)..r {
        if h[i] < v {
            v = h[i];
            idx = i;
        }
    }
    let start = idx.saturating_sub(2);
    let end = (idx + 2).min(h.len() - 1);
    let mut best = idx;
    let mut bv = h[idx];
    for i in start..=end {
        if h[i] < bv {
            bv = h[i];
            best = i;
        }
    }
    Some(best)
}

#[inline]
fn quantile_hist_range(
    hist: &[u32],
    min_log: f32,
    bin_w: f32,
    lo: usize,
    hi: usize,
    q: f32,
) -> f32 {
    if lo >= hist.len() || hi == 0 || lo > hi {
        return f32::INFINITY;
    }
    let mut total: u64 = 0;
    for i in lo..=hi {
        total += hist[i] as u64;
    }
    if total == 0 {
        return f32::INFINITY;
    }
    let target = (total as f64 * q as f64).ceil() as u64;
    let mut acc: u64 = 0;
    let mut idx = hi;
    for i in lo..=hi {
        acc += hist[i] as u64;
        if acc >= target {
            idx = i;
            break;
        }
    }
    10f32.powf(bin_center_log(min_log, bin_w, idx))
}

fn otsu_threshold_log(hist: &[u32], min_log: f32, bin_w: f32) -> f32 {
    let mut total: u64 = 0;
    let mut sum1: u64 = 0;
    for (i, &h) in hist.iter().enumerate() {
        total += h as u64;
        sum1 += (i as u64) * (h as u64);
    }
    if total == 0 {
        return f32::INFINITY;
    }
    let mut w_b: u64 = 0;
    let mut sum_b: u64 = 0;
    let mut max_var: f64 = -1.0;
    let mut idx = 0usize;
    for i in 0..hist.len() {
        w_b += hist[i] as u64;
        if w_b == 0 {
            continue;
        }
        let w_f = total - w_b;
        if w_f == 0 {
            break;
        }
        sum_b += (i as u64) * (hist[i] as u64);
        let m_b = sum_b as f64 / w_b as f64;
        let m_f = (sum1 - sum_b) as f64 / w_f as f64;
        let var_b = (w_b as f64) * (w_f as f64) * (m_b - m_f) * (m_b - m_f);
        if var_b > max_var {
            max_var = var_b;
            idx = i;
        }
    }
    10f32.powf(bin_center_log(min_log, bin_w, idx))
}

pub fn find_noise_level(intensities: &[f32]) -> f32 {
    find_noise_level_with_options(intensities, None)
}

pub fn find_noise_level_with_options(
    intensities: &[f32],
    peaks_opts: Option<ScanPeaksOptions>,
) -> f32 {
    const MIN_SEP_LOG: f32 = 0.20;
    const MIN_VALLEY_DEPTH: f32 = 0.30;
    const MIN_SCORE: f32 = 0.55;
    const MIN_REL_HEIGHT: f32 = 0.01;
    const LOW_BASELINE_MAX: f32 = 2.2;

    let mut n_pos = 0usize;
    let mut min_log = f32::INFINITY;
    let mut max_log = f32::NEG_INFINITY;
    for &a in intensities {
        if is_finite_positive(a) {
            n_pos += 1;
            let lg = a.log10();
            if lg < min_log {
                min_log = lg;
            }
            if lg > max_log {
                max_log = lg;
            }
        }
    }
    if n_pos == 0 || !min_log.is_finite() || !max_log.is_finite() || max_log <= min_log {
        return f32::INFINITY;
    }

    let bins = choose_bins(n_pos);
    let range = max_log - min_log;
    let inv_bin_w = (bins as f32) / range;
    let mut hist = vec![0u32; bins];
    for &a in intensities {
        if is_finite_positive(a) {
            let lg = a.log10();
            let mut idx = ((lg - min_log) * inv_bin_w) as isize;
            if idx < 0 {
                idx = 0;
            }
            let mut idx = idx as usize;
            if idx >= bins {
                idx = bins - 1;
            }
            hist[idx] = hist[idx].saturating_add(1);
        }
    }

    let bin_w = range / bins as f32;
    let mut xs = Vec::with_capacity(bins);
    let mut ys = Vec::with_capacity(bins);
    for i in 0..bins {
        xs.push(bin_center_log(min_log, bin_w, i) as f64);
        ys.push(hist[i] as f32);
    }

    let base_opts = peaks_opts.unwrap_or_default();
    let mut ws = base_opts.window_size.max(5);
    if ws > bins {
        ws = if bins % 2 == 1 { bins } else { bins - 1 };
    }
    if ws % 2 == 0 {
        ws -= 1;
    }
    let fp = ScanPeaksOptions {
        epsilon: base_opts.epsilon,
        window_size: ws,
    };

    let data = DataXY { x: xs, y: ys };
    let peak_xs = scan_for_peaks(&data, Some(fp));
    if peak_xs.len() < 2 {
        return quantile_hist_range(&hist, min_log, bin_w, 0, bins - 1, 0.99);
    }

    let s0 = SggOptions {
        window_size: ws,
        derivative: 0,
        polynomial: 3,
    };
    let ys_sm = if bins >= 5 {
        sgg(&data.y, &data.x, s0)
    } else {
        data.y.clone()
    };

    let mut peak_idx: Vec<usize> = peak_xs
        .iter()
        .map(|&x| x_to_nearest_bin(x, min_log, bin_w, bins))
        .collect();
    peak_idx.sort_unstable();
    peak_idx.dedup();
    if peak_idx.len() < 2 {
        return quantile_hist_range(&hist, min_log, bin_w, 0, bins - 1, 0.99);
    }

    let mut noise_pk = peak_idx[0];
    let mut noise_h = ys_sm[noise_pk];
    for &p in &peak_idx {
        let h = ys_sm[p];
        if h > noise_h {
            noise_pk = p;
            noise_h = h;
        }
    }

    let noise_center = 10f32.powf(bin_center_log(min_log, bin_w, noise_pk));

    let mut rights: Vec<(usize, f32, f32)> = Vec::new();
    for &p in &peak_idx {
        if p > noise_pk {
            let sep = (p - noise_pk) as f32 * bin_w;
            let h = ys_sm[p];
            if sep >= MIN_SEP_LOG && h >= noise_h * MIN_REL_HEIGHT {
                rights.push((p, h, sep));
            }
        }
    }

    if rights.len() == 1 && noise_center < LOW_BASELINE_MAX {
        return quantile_hist_range(&hist, min_log, bin_w, 0, bins - 1, 0.999);
    }
    if rights.is_empty() {
        return quantile_hist_range(&hist, min_log, bin_w, 0, bins - 1, 0.99);
    }

    rights.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
    let (an_pk, an_h, _) = rights[0];

    let Some(v_idx) = valley_between(&ys_sm, noise_pk, an_pk) else {
        return quantile_hist_range(&hist, min_log, bin_w, 0, bins - 1, 0.99);
    };

    let ref_peak_max = noise_h.max(an_h).max(1.0);
    let valley_depth_ratio = clamp01((ref_peak_max - ys_sm[v_idx]) / ref_peak_max);
    let sep_log = (an_pk - noise_pk) as f32 * bin_w;
    let sep_score = clamp01((sep_log - MIN_SEP_LOG) / 0.70);
    let reliable = 0.6 * valley_depth_ratio + 0.4 * sep_score >= MIN_SCORE
        && valley_depth_ratio >= MIN_VALLEY_DEPTH;

    let noise_log = bin_center_log(min_log, bin_w, noise_pk);
    let analyte_log = bin_center_log(min_log, bin_w, an_pk);
    let valley_thr = 10f32.powf(bin_center_log(min_log, bin_w, v_idx));
    let gm_thr = 10f32.powf(0.5 * (noise_log + analyte_log));

    let base = if reliable {
        valley_thr.max(gm_thr)
    } else {
        quantile_hist_range(&hist, min_log, bin_w, 0, v_idx, 0.995).max(valley_thr.min(gm_thr))
    };

    let otsu = otsu_threshold_log(&hist, min_log, bin_w);
    if base > otsu * 1.8 { otsu } else { base }
}
