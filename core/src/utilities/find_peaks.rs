use crate::utilities::closest_index;
use crate::utilities::find_noise_level::find_noise_level;
use crate::utilities::get_boundaries::{BoundariesOptions, get_boundaries};
use crate::utilities::scan_for_peaks::{ScanPeaksOptions, scan_for_peaks_across_windows};
use crate::utilities::sgg::{SggOptions, sgg};
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

    let mut candidates: Vec<PeakCandidate> = Vec::with_capacity(positions.len());
    let mut sum_integrals = 0.0f64;
    let use_ratio = filter_opts.integral_threshold.is_some();

    for rt in positions {
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
            ratio: 0.0,
            noise: resolved_noise,
        });
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

    if peaks.len() > 1 {
        peaks = merge_apex_wiggles(data, peaks, resolved_noise);
    }
    if peaks.len() > 1 {
        peaks = merge_tail_groups(data, peaks, resolved_noise);
    }

    if !filter_opts.allow_overlap.unwrap_or(false) && peaks.len() > 1 {
        peaks = prune_overlaps_by_valley(data, peaks, resolved_noise);
    }

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
            peaks.retain(|p| p.intensity >= cutoff);
        }
    }

    if !peaks.is_empty() {
        peaks = extend_tails(data, peaks, resolved_noise);
        peaks = coalesce_tiny_tail_bumps(data, peaks, resolved_noise);
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

fn merge_apex_wiggles(data: &DataXY, peaks: Vec<Peak>, noise: f64) -> Vec<Peak> {
    if peaks.len() <= 1 {
        return peaks;
    }

    let n = data.y.len();
    let mut ws = 33usize;
    if ws > n || ws < 5 {
        ws = if n >= 5 {
            if n % 2 == 1 { n } else { n - 1 }
        } else {
            21
        };
    }
    if ws % 2 == 0 {
        ws -= 1;
    }
    if ws < 5 {
        ws = 21;
    }

    let ys = if ws <= n {
        sgg(
            &data.y,
            &data.x,
            SggOptions {
                window_size: ws,
                derivative: 0,
                polynomial: 3,
            },
        )
    } else {
        data.y.clone()
    };

    let mut out: Vec<Peak> = Vec::with_capacity(peaks.len());
    let mut cur_opt: Option<Peak> = None;
    let min_rel_drop: f32 = 0.08;

    for mut nxt in peaks.into_iter() {
        if let Some(mut cur) = cur_opt.take() {
            let i0 = closest_index(&data.x, cur.rt);
            let i1 = closest_index(&data.x, nxt.rt);
            let l = i0.min(i1);
            let r = i0.max(i1);

            let mut valley = f32::INFINITY;
            let mut j = l;
            while j <= r {
                let v = ys[j];
                if v < valley {
                    valley = v;
                }
                j += 1;
            }

            let h_max = (cur.intensity.max(nxt.intensity)) as f32;
            let rel_drop = if h_max > 0.0 {
                (h_max - valley) / h_max
            } else {
                0.0
            };
            let valley_gate = (noise as f32).max(0.92 * h_max);
            let merge = valley > valley_gate || rel_drop < min_rel_drop;

            if merge {
                if nxt.intensity > cur.intensity {
                    nxt.from = nxt.from.min(cur.from);
                    nxt.to = nxt.to.max(cur.to);
                    nxt.integral += cur.integral;
                    nxt.np += cur.np;
                    cur_opt = Some(nxt);
                } else {
                    cur.from = cur.from.min(nxt.from);
                    cur.to = cur.to.max(nxt.to);
                    cur.integral += nxt.integral;
                    cur.np += nxt.np;
                    cur_opt = Some(cur);
                }
            } else {
                out.push(cur);
                cur_opt = Some(nxt);
            }
        } else {
            cur_opt = Some(nxt);
        }
    }

    if let Some(last) = cur_opt {
        out.push(last);
    }
    out
}

fn merge_tail_groups(data: &DataXY, peaks: Vec<Peak>, noise: f64) -> Vec<Peak> {
    if peaks.len() <= 1 {
        return peaks;
    }

    let n = data.y.len();
    let mut ws = 33usize;
    if ws > n || ws < 5 {
        ws = if n >= 5 {
            if n % 2 == 1 { n } else { n - 1 }
        } else {
            21
        };
    }
    if ws % 2 == 0 {
        ws -= 1;
    }
    if ws < 5 {
        ws = 21;
    }

    let ys = if ws <= n {
        sgg(
            &data.y,
            &data.x,
            SggOptions {
                window_size: ws,
                derivative: 0,
                polynomial: 3,
            },
        )
    } else {
        data.y.clone()
    };
    let dy1 = if ws <= n {
        sgg(
            &data.y,
            &data.x,
            SggOptions {
                window_size: ws,
                derivative: 1,
                polynomial: 3,
            },
        )
    } else {
        vec![0.0f32; n]
    };
    let d2 = if ws <= n {
        sgg(
            &data.y,
            &data.x,
            SggOptions {
                window_size: ws,
                derivative: 2,
                polynomial: 3,
            },
        )
    } else {
        vec![0.0f32; n]
    };

    let dx_avg = if n > 1 {
        ((data.x[n - 1] - data.x[0]).abs()) / ((n as f64) - 1.0)
    } else {
        1.0
    };
    let near_gap = 2.0 * dx_avg;

    let mut ps = peaks.clone();
    ps.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));

    let mut out: Vec<Peak> = Vec::with_capacity(ps.len());
    let mut cur = ps[0].clone();

    let slope_base = (3.0 * noise / dx_avg.max(f64::EPSILON)) as f32;
    let curv_base = (1.5 * noise / (dx_avg * dx_avg).max(f64::EPSILON)) as f32;

    for k in 1..ps.len() {
        let nxt = &ps[k];

        let i0 = closest_index(&data.x, cur.rt);
        let i1 = closest_index(&data.x, nxt.rt);
        let l = i0.min(i1);
        let r = i0.max(i1);

        let mut tvp = 0.0f32;
        let mut tvn = 0.0f32;
        let mut minv = f32::INFINITY;
        let mut j = l;
        while j + 1 <= r {
            let d = ys[j + 1] - ys[j];
            if d > 0.0 {
                tvp += d;
            } else {
                tvn += -d;
            }
            if ys[j] < minv {
                minv = ys[j];
            }
            j += 1;
        }
        if ys[r] < minv {
            minv = ys[r];
        }
        let frac_pos = tvp / (tvp + tvn + f32::EPSILON);

        let kd = if cur.intensity > 0.0 {
            ((cur.intensity as f32) - minv).max(0.0) / (cur.intensity as f32)
        } else {
            0.0
        };
        let pd = if nxt.intensity > 0.0 {
            ((nxt.intensity as f32) - minv).max(0.0) / (nxt.intensity as f32)
        } else {
            0.0
        };
        let deep_valley = kd >= 0.35 && pd >= 0.35;

        let base_thr = (2.0 * noise) as f32;
        let mut base_run = 0usize;
        let mut base_max = 0usize;
        let mut t = l;
        while t <= r {
            if ys[t] <= base_thr {
                base_run += 1;
            } else {
                base_max = base_max.max(base_run);
                base_run = 0;
            }
            t += 1;
        }
        base_max = base_max.max(base_run);
        let long_baseline = base_max as f64 >= 10.0;

        let cur_width = (cur.to - cur.from).abs().max(dx_avg);
        let nxt_width = (nxt.to - nxt.from).abs().max(dx_avg);
        let span_rt = nxt.rt - cur.rt;
        let within_span =
            span_rt <= 4.0 * cur_width && (nxt.from - cur.to).abs() <= 2.5 * cur_width;

        let local_ratio = if cur.intensity > 0.0 {
            nxt.intensity / cur.intensity
        } else {
            1.0
        };
        let sn_next = if noise > 0.0 {
            nxt.intensity / noise
        } else {
            f64::INFINITY
        };
        let weak_local = local_ratio < 0.45 && sn_next < 10.0 && nxt_width < 0.55 * cur_width;

        let contiguous = (nxt.from - cur.to).abs() <= near_gap;
        let valley_close = if nxt.intensity > 0.0 {
            (((nxt.intensity as f32) - minv).max(0.0) / (nxt.intensity as f32)) < 0.30
        } else {
            false
        };
        let tail_like = frac_pos < 0.22;

        let tiny_rel = local_ratio < 0.02;
        let tiny_abs = sn_next < 6.0;
        let narrow = nxt_width < 0.70 * cur_width;
        let tail_bump_override = within_span && contiguous && narrow && (tiny_rel || tiny_abs);

        let mut rise_run = 0usize;
        let mut rise_veto = false;
        if !tail_bump_override {
            let need_rise = 6usize;
            let rise_amp = ((0.25 * cur.intensity) as f32).max((12.0 * noise) as f32);
            let mut u = l;
            while u <= r {
                let v = ys[u];
                let rising = dy1[u] > slope_base && d2[u] > curv_base && (v - minv) > rise_amp;
                if rising {
                    rise_run += 1;
                } else {
                    rise_run = 0;
                }
                if rise_run >= need_rise {
                    rise_veto = true;
                    break;
                }
                u += 1;
            }
        }

        let base_merge = !rise_veto
            && !deep_valley
            && !long_baseline
            && within_span
            && weak_local
            && tail_like
            && (contiguous || valley_close);
        let do_merge = tail_bump_override || base_merge;

        if do_merge {
            let mut merged = cur.clone();
            if nxt.intensity > merged.intensity {
                merged.intensity = nxt.intensity;
                merged.rt = nxt.rt;
            }
            merged.to = merged.to.max(nxt.to);
            merged.integral += nxt.integral;
            merged.np += nxt.np;
            cur = merged;
        } else {
            out.push(cur);
            cur = nxt.clone();
        }
    }
    out.push(cur);
    out
}

fn prune_overlaps_by_valley(data: &DataXY, peaks: Vec<Peak>, noise: f64) -> Vec<Peak> {
    if peaks.len() <= 1 {
        return peaks;
    }

    let n = data.y.len();
    let ws0 = 21usize;
    let ws = if ws0 <= n {
        if ws0 % 2 == 1 { ws0 } else { ws0 - 1 }
    } else {
        if n % 2 == 1 { n } else { n.saturating_sub(1) }
    };
    let ys = if ws >= 5 && ws <= n {
        sgg(
            &data.y,
            &data.x,
            SggOptions {
                window_size: ws,
                derivative: 0,
                polynomial: 3,
            },
        )
    } else {
        data.y.clone()
    };

    let mut dx_min = f64::INFINITY;
    for w in data.x.windows(2) {
        let d = w[1] - w[0];
        if d > 0.0 && d < dx_min {
            dx_min = d;
        }
    }
    let min_sep = if dx_min.is_finite() {
        1.05 * dx_min
    } else {
        0.02
    };
    let mut ps = peaks.clone();
    ps.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));

    let mut out: Vec<Peak> = Vec::with_capacity(ps.len());
    let mut cur: Option<Peak> = None;

    for mut p in ps {
        if let Some(mut k) = cur.take() {
            let overlaps = (p.from < k.to) && (p.to > k.from);
            let very_close = (p.rt - k.rt).abs() < min_sep;

            if overlaps || very_close {
                let i0 = closest_index(&data.x, k.rt);
                let i1 = closest_index(&data.x, p.rt);
                let l = i0.min(i1);
                let r = i0.max(i1);

                let mut valley = f32::INFINITY;
                let mut vpos = l;
                let mut j = l;
                while j <= r {
                    let v = ys[j];
                    if v < valley {
                        valley = v;
                        vpos = j;
                    }
                    j += 1;
                }

                let kd = if k.intensity > 0.0 {
                    ((k.intensity as f32) - valley) / (k.intensity as f32)
                } else {
                    0.0
                };
                let pd = if p.intensity > 0.0 {
                    ((p.intensity as f32) - valley) / (p.intensity as f32)
                } else {
                    0.0
                };

                let ratio = if k.intensity > p.intensity {
                    p.intensity / k.intensity
                } else {
                    k.intensity / p.intensity
                };
                let req_drop = if ratio < 0.40 { 0.25 } else { 0.08 };

                let keep_both = ((kd >= req_drop) && (pd >= req_drop)) || (valley as f64) <= noise;

                if keep_both {
                    let vx = data.x[vpos];

                    if vx > k.from && vx < k.to {
                        k.to = vx;
                        let k_li = closest_index(&data.x, k.from);
                        let k_ri = closest_index(&data.x, k.to);
                        if k_ri >= k_li {
                            k.np = (k_ri - k_li + 1) as i32;
                        }
                    }

                    if vx > p.from && vx < p.to {
                        p.from = vx;
                        let p_li = closest_index(&data.x, p.from);
                        let p_ri = closest_index(&data.x, p.to);
                        if p_ri >= p_li {
                            p.np = (p_ri - p_li + 1) as i32;
                        }
                    }

                    out.push(k);
                    cur = Some(p);
                } else {
                    let choose_p = if p.intensity > k.intensity {
                        true
                    } else if (p.intensity - k.intensity).abs() < f64::EPSILON {
                        if p.integral > k.integral {
                            true
                        } else if (p.integral - k.integral).abs() < f64::EPSILON {
                            let p_w = (p.to - p.from).abs();
                            let k_w = (k.to - k.from).abs();
                            p_w < k_w
                        } else {
                            false
                        }
                    } else {
                        false
                    };
                    cur = Some(if choose_p { p } else { k });
                }
            } else {
                out.push(k);
                cur = Some(p);
            }
        } else {
            cur = Some(p);
        }
    }

    if let Some(k) = cur {
        out.push(k);
    }
    out
}

fn extend_tails(data: &DataXY, mut peaks: Vec<Peak>, noise: f64) -> Vec<Peak> {
    let n = data.y.len();
    if n == 0 || peaks.is_empty() {
        return peaks;
    }

    let mut ws = 51usize;
    if ws > n || ws < 5 {
        ws = if n >= 5 {
            if n % 2 == 1 { n } else { n - 1 }
        } else {
            21
        };
    }
    if ws % 2 == 0 {
        ws -= 1;
    }
    if ws < 5 {
        ws = 21;
    }

    let ys0 = if ws <= n {
        sgg(
            &data.y,
            &data.x,
            SggOptions {
                window_size: ws,
                derivative: 0,
                polynomial: 3,
            },
        )
    } else {
        data.y.clone()
    };
    let dy1 = if ws <= n {
        sgg(
            &data.y,
            &data.x,
            SggOptions {
                window_size: ws,
                derivative: 1,
                polynomial: 3,
            },
        )
    } else {
        vec![0.0f32; n]
    };
    let d2 = if ws <= n {
        sgg(
            &data.y,
            &data.x,
            SggOptions {
                window_size: ws,
                derivative: 2,
                polynomial: 3,
            },
        )
    } else {
        vec![0.0f32; n]
    };

    let dx_avg = if n > 1 {
        ((data.x[n - 1] - data.x[0]).abs()) / ((n as f64) - 1.0)
    } else {
        1.0
    };

    peaks.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));

    for i in 0..peaks.len() {
        let mut p = peaks[i].clone();

        let floor = ((1.5 * noise) as f32).max((0.003 * p.intensity) as f32);
        let rise_amp = ((8.0 * noise) as f32).max((0.10 * p.intensity) as f32);
        let slope_gate = (3.0 * noise / dx_avg.max(f64::EPSILON)) as f32;
        let curv_gate = (1.5 * noise / (dx_avg * dx_avg).max(f64::EPSILON)) as f32;

        let need_below = 14usize;
        let need_rise = 7usize;

        let ri0 = closest_index(&data.x, p.to);
        let li0 = closest_index(&data.x, p.from);

        let guard_left = if i > 0 {
            closest_index(&data.x, peaks[i - 1].to)
        } else {
            0
        };
        let guard_right = if i + 1 < peaks.len() {
            closest_index(&data.x, peaks[i + 1].from)
        } else {
            n - 1
        };

        if li0 > guard_left {
            let mut minv = ys0[li0];
            let mut minpos = li0;
            let mut below_run = 0usize;
            let mut rise_run = 0usize;
            let mut stop_at = li0;
            let mut j = li0;
            while j > guard_left {
                let v = ys0[j];
                if v < minv {
                    minv = v;
                    minpos = j;
                }
                if v <= floor {
                    below_run += 1;
                } else {
                    below_run = 0;
                }
                if below_run >= need_below {
                    stop_at = j;
                    break;
                }
                let s_left = -dy1[j];
                let c2 = d2[j];
                let rising_left = s_left > slope_gate && c2 > curv_gate && (v - minv) > rise_amp;
                if rising_left {
                    rise_run += 1;
                } else {
                    rise_run = 0;
                }
                if rise_run >= need_rise {
                    stop_at = minpos;
                    break;
                }
                stop_at = j.saturating_sub(1);
                j = j.saturating_sub(1);
            }
            let new_from_i = stop_at.max(guard_left);
            if new_from_i < li0 {
                let (area, _) =
                    xy_integration(&data.x[new_from_i..=ri0], &data.y[new_from_i..=ri0]);
                p.from = data.x[new_from_i];
                p.integral = area;
                p.np = (ri0 - new_from_i + 1) as i32;
            }
        }

        let li1 = closest_index(&data.x, p.from);
        let ri1 = closest_index(&data.x, p.to);

        if ri1 < guard_right {
            let mut minv = ys0[ri1];
            let mut minpos = ri1;
            let mut below_run = 0usize;
            let mut rise_run = 0usize;
            let mut stop_at = ri1;
            let mut j = ri1;
            while j + 1 <= guard_right {
                let v = ys0[j];
                if v < minv {
                    minv = v;
                    minpos = j;
                }
                if v <= floor {
                    below_run += 1;
                } else {
                    below_run = 0;
                }
                if below_run >= need_below {
                    stop_at = j;
                    break;
                }
                let s1 = dy1[j];
                let c2 = d2[j];
                let rising = s1 > slope_gate && c2 > curv_gate && (v - minv) > rise_amp;
                if rising {
                    rise_run += 1;
                } else {
                    rise_run = 0;
                }
                if rise_run >= need_rise {
                    stop_at = minpos;
                    break;
                }
                stop_at = j + 1;
                j += 1;
            }
            let new_to_i = stop_at.min(guard_right);
            if new_to_i > li1 {
                let (area, _) = xy_integration(&data.x[li1..=new_to_i], &data.y[li1..=new_to_i]);
                p.to = data.x[new_to_i];
                p.integral = area;
                p.np = (new_to_i - li1 + 1) as i32;
            }
        }

        peaks[i] = p;
    }

    peaks
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

fn coalesce_tiny_tail_bumps(data: &DataXY, peaks: Vec<Peak>, noise: f64) -> Vec<Peak> {
    if peaks.len() <= 1 {
        return peaks;
    }

    let n = data.x.len();
    let dx_avg = if n > 1 {
        ((data.x[n - 1] - data.x[0]).abs()) / ((n as f64) - 1.0)
    } else {
        1.0
    };
    let near_gap = 2.0 * dx_avg;

    let mut ps = peaks.clone();
    ps.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal));

    let mut out: Vec<Peak> = Vec::with_capacity(ps.len());
    let mut cur = ps[0].clone();

    for k in 1..ps.len() {
        let p = &ps[k];

        let gap = (p.from - cur.to).abs();
        let cur_w = (cur.to - cur.from).abs().max(dx_avg);
        let p_w = (p.to - p.from).abs().max(dx_avg);

        let rel = if cur.intensity > 0.0 {
            p.intensity / cur.intensity
        } else {
            1.0
        };
        let tiny = rel < 0.02 || (noise > 0.0 && p.intensity < 6.0 * noise);
        let narrow = p_w < 0.80 * cur_w;

        if gap <= near_gap && tiny && narrow {
            let mut merged = cur.clone();
            merged.to = merged.to.max(p.to);
            merged.integral += p.integral;
            merged.np += p.np;
            if p.intensity > merged.intensity {
                merged.intensity = p.intensity;
                merged.rt = p.rt;
            }
            cur = merged;
        } else {
            out.push(cur);
            cur = p.clone();
        }
    }
    out.push(cur);
    out
}
