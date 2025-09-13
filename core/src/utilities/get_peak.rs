use std::panic::{AssertUnwindSafe, catch_unwind};

use crate::utilities::find_noise_level::find_noise_level;
use crate::utilities::find_peaks::FilterPeaksOptions;
use crate::utilities::find_peaks::{FindPeaksOptions, find_peaks};
use crate::utilities::scan_for_peaks::ScanPeaksOptions;
use crate::utilities::sgg::{SggOptions, sgg};
use crate::utilities::structs::{DataXY, Peak, Roi};

const DEFAULT_WINDOW_SIZES: &[usize] = &[5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33];

fn series_len(data: &DataXY) -> usize {
    data.x.len().min(data.y.len())
}

fn odd_at_most(ws: usize, n: usize) -> Option<usize> {
    if n < 5 || ws < 5 || ws > n {
        return None;
    }
    let w = if ws % 2 == 0 { ws - 1 } else { ws };
    if w >= 5 { Some(w) } else { None }
}

fn min_spacing(xs: &[f64]) -> f64 {
    if xs.len() < 2 {
        return 0.01;
    }
    let mut dx = f64::INFINITY;
    for w in xs.windows(2) {
        let d = w[1] - w[0];
        if d > 0.0 && d < dx {
            dx = d;
        }
    }
    if dx.is_finite() { dx } else { 0.01 }
}

fn closest_index(xs: &[f64], v: f64) -> usize {
    if xs.is_empty() {
        return 0;
    }
    let mut lo = 0usize;
    let mut hi = xs.len();
    while lo < hi {
        let m = (lo + hi) / 2;
        if xs[m] < v {
            lo = m + 1;
        } else {
            hi = m;
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

    let w = if roi.window.is_finite() && roi.window > 0.0 {
        roi.window
    } else {
        0.0
    };
    let lo = roi.rt - w;
    let hi = roi.rt + w;

    let mut sizes = Vec::new();
    for &ws in wss {
        if let Some(e) = odd_at_most(ws, n) {
            sizes.push(e);
        }
    }
    if sizes.is_empty() {
        return None;
    }

    let raw_noise = find_noise_level(&data.y).max(0.0) as f64;
    let gate_noise = 3.0 * raw_noise;

    let mut peaks_all: Vec<Peak> = Vec::new();
    let mut win_id: Vec<usize> = Vec::new();

    for (wi, &ws_eff) in sizes.iter().enumerate() {
        let base_scan = base
            .scan_peaks_options
            .as_ref()
            .map(|s| ScanPeaksOptions {
                epsilon: s.epsilon,
                window_size: s.window_size,
            })
            .unwrap_or_default();

        let tuned_filter = Some(FilterPeaksOptions {
            integral_threshold: None,
            width_threshold: Some(3),
            intensity_threshold: None,
            noise: Some(gate_noise),
            auto_noise: Some(false),
            allow_overlap: Some(true),
            sn_ratio: Some(0),
        });

        let opts = FindPeaksOptions {
            get_boundaries_options: base.get_boundaries_options,
            filter_peaks_options: tuned_filter,
            scan_peaks_options: Some(ScanPeaksOptions {
                epsilon: base_scan.epsilon,
                window_size: ws_eff,
            }),
        };

        let peaks = match catch_unwind(AssertUnwindSafe(|| find_peaks(data, Some(opts)))) {
            Ok(v) => v,
            Err(_) => Vec::new(),
        };

        let in_roi: Vec<Peak> = if w > 0.0 {
            peaks
                .into_iter()
                .filter(|p| p.rt.is_finite() && p.rt >= lo && p.rt <= hi)
                .collect()
        } else {
            peaks
        };

        for p in in_roi {
            peaks_all.push(p);
            win_id.push(wi);
        }
    }
    if peaks_all.is_empty() {
        return None;
    }

    let ys_sm: Vec<f32> = sgg(
        &data.y,
        &data.x,
        SggOptions {
            window_size: 11.min(n.max(1)),
            derivative: 0,
            polynomial: 3,
        },
    );
    let noise = find_noise_level(&data.y).max(0.0) as f32;
    let valley_floor = 6.0 * noise;
    const VALLEY_FRAC: f32 = 0.7;

    let dx = min_spacing(&data.x);
    let rt_gate = (1.2 * dx).max(0.006);

    let mut order: Vec<usize> = (0..peaks_all.len()).collect();
    order.sort_by(|&i, &j| peaks_all[i].rt.total_cmp(&peaks_all[j].rt));

    let mut centers: Vec<f64> = Vec::new();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut group_maxy: Vec<f64> = Vec::new();
    let mut group_ws: Vec<Vec<usize>> = Vec::new();

    for idx in order {
        let p = &peaks_all[idx];
        let mut placed = false;
        if let Some(g) = centers.len().checked_sub(1) {
            let dt = (p.rt - centers[g]).abs();
            if dt <= rt_gate {
                let i0 = closest_index(&data.x, p.rt);
                let i1 = closest_index(&data.x, centers[g]);
                let l = i0.min(i1);
                let r = i0.max(i1);
                let mut valley = f32::INFINITY;
                for j in l..=r {
                    let v = ys_sm[j];
                    if v < valley {
                        valley = v;
                    }
                }
                let thr = (VALLEY_FRAC * (p.intensity.min(group_maxy[g]) as f32)).max(valley_floor);
                if valley >= thr {
                    groups[g].push(idx);
                    centers[g] = (centers[g] * (groups[g].len() as f64 - 1.0) + p.rt)
                        / (groups[g].len() as f64);
                    if p.intensity > group_maxy[g] {
                        group_maxy[g] = p.intensity;
                    }
                    let wtag = win_id[idx];
                    if !group_ws[g].contains(&wtag) {
                        group_ws[g].push(wtag);
                    }
                    placed = true;
                }
            }
        }
        if !placed {
            centers.push(p.rt);
            groups.push(vec![idx]);
            group_maxy.push(p.intensity);
            group_ws.push(vec![win_id[idx]]);
        }
    }

    let r_search = if w > 0.0 {
        (0.25 * w).max(3.0 * dx)
    } else {
        3.0 * dx
    };
    let i_lo = closest_index(&data.x, roi.rt - r_search);
    let i_hi = closest_index(&data.x, roi.rt + r_search);
    let mut im = i_lo;
    let mut ym = f32::NEG_INFINITY;
    let mut i = i_lo;
    while i <= i_hi && i < ys_sm.len() {
        if ys_sm[i] > ym {
            ym = ys_sm[i];
            im = i;
        }
        i += 1;
    }
    let ort_rt = data.x[im];
    let mut ort_cluster: Option<usize> = None;
    for gi in 0..groups.len() {
        for &ix in &groups[gi] {
            if (peaks_all[ix].rt - ort_rt).abs() <= rt_gate {
                ort_cluster = Some(gi);
                break;
            }
        }
        if ort_cluster.is_some() {
            break;
        }
    }
    if ort_cluster.is_none() {
        let mut bd = f64::INFINITY;
        let mut bi = None;
        for gi in 0..groups.len() {
            let d = (centers[gi] - ort_rt).abs();
            if d < bd {
                bd = d;
                bi = Some(gi);
            }
        }
        ort_cluster = bi;
    }

    let m = sizes.len();
    let need_support = if m <= 3 {
        1
    } else {
        ((m as f64) * 0.25).ceil() as usize
    };

    let mut inside: Vec<usize> = Vec::new();
    for i in 0..groups.len() {
        let c = centers[i];
        if w > 0.0 && c >= lo && c <= hi {
            inside.push(i);
        }
    }
    let base_cands: Vec<usize> = if !inside.is_empty() {
        inside
    } else {
        (0..groups.len()).collect()
    };

    let mut roi_max_h = f64::NEG_INFINITY;
    for &i in &base_cands {
        if group_maxy[i] > roi_max_h {
            roi_max_h = group_maxy[i];
        }
    }
    let min_h_abs = 6.0 * noise as f64;
    let min_h_rel = 0.10 * roi_max_h;
    let min_h = if min_h_abs > min_h_rel {
        min_h_abs
    } else {
        min_h_rel
    };

    let mut filtered: Vec<usize> = Vec::new();
    for &i in &base_cands {
        if group_ws[i].len() >= need_support && group_maxy[i] >= min_h {
            filtered.push(i);
        }
    }
    let mut pool: Vec<usize> = if !filtered.is_empty() {
        filtered
    } else {
        base_cands
    };

    if let Some(gi) = ort_cluster {
        if !pool.contains(&gi) {
            pool.push(gi);
        }
    }

    let mut best_i: Option<usize> = None;
    let mut best_d = f64::INFINITY;
    let mut best_sup = 0usize;
    let mut best_h = f64::NEG_INFINITY;
    let dist_eps = 0.25 * dx;

    for i in pool {
        let d = (centers[i] - roi.rt).abs();
        let sup = group_ws[i].len();
        let h = group_maxy[i];
        let better = if best_i.is_none() {
            true
        } else if d + dist_eps < best_d {
            true
        } else if (d - best_d).abs() <= dist_eps {
            if sup != best_sup {
                sup > best_sup
            } else {
                h.total_cmp(&best_h).is_gt()
            }
        } else {
            false
        };
        if better {
            best_i = Some(i);
            best_d = d;
            best_sup = sup;
            best_h = h;
        }
    }

    let g = match best_i {
        Some(v) => v,
        None => return None,
    };

    let mut best: Option<Peak> = None;
    let mut h = f64::NEG_INFINITY;
    for &ix in &groups[g] {
        let y = peaks_all[ix].intensity;
        if y > h {
            h = y;
            best = Some(peaks_all[ix].clone());
        }
    }
    best
}

pub fn get_peak(data: &DataXY, roi: Roi, options: Option<FindPeaksOptions>) -> Option<Peak> {
    get_peak_across_windows(data, roi, options, None)
}
