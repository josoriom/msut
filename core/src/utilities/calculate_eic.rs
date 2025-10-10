use std::{cmp::Ordering, sync::Arc};

use crate::utilities::{
    parse::{decode::decode, parse_mzml::MzML},
    structs::{FromTo, Peak},
};

#[derive(Clone, Copy)]
pub struct EicOptions {
    pub ppm_tolerance: f64,
    pub mz_tolerance: f64,
}

impl Default for EicOptions {
    fn default() -> Self {
        Self {
            ppm_tolerance: 20.0,
            mz_tolerance: 0.005,
        }
    }
}

pub struct Eic {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

pub fn calculate_eic_from_bin1(
    bin1: &[u8],
    target_mass: &f64,
    from_to: FromTo,
    options: EicOptions,
) -> Result<Eic, &'static str> {
    let mzml = decode(bin1).map_err(|_| "decode BIN1 failed")?;
    calculate_eic_from_mzml(&mzml, target_mass, from_to, options)
}

pub fn calculate_eic_from_mzml(
    mzml: &MzML,
    target_mass: &f64,
    from_to: FromTo,
    options: EicOptions,
) -> Result<Eic, &'static str> {
    let (times, scans) = collect_ms1_scans(mzml, from_to);
    if scans.is_empty() || times.is_empty() {
        return Ok(Eic {
            x: Vec::new(),
            y: Vec::new(),
        });
    }
    let y = compute_eic_for_mz(&scans, times.len(), target_mass, options);
    Ok(Eic { x: times, y })
}

#[derive(Clone)]
pub struct CentroidScan {
    pub rt: f64,
    pub mz: Arc<[f64]>,
    pub intensity: Arc<[f64]>,
}

pub fn compute_eic_for_mz(
    scans: &[CentroidScan],
    rt_len: usize,
    center: &f64,
    opts: EicOptions,
) -> Vec<f64> {
    let tol_ppm = if opts.ppm_tolerance > 0.0 {
        (opts.ppm_tolerance * 1e-6) * center
    } else {
        0.0
    };
    let tol = tol_ppm.max(opts.mz_tolerance.max(0.0));
    if !(tol.is_finite()) || tol <= 0.0 {
        panic!("[panic] invalid EIC tol for center={}", center);
    }
    let lo = center - tol;
    let hi = center + tol;

    let mut y = vec![0.0f64; rt_len];
    for (i, s) in scans.iter().enumerate() {
        let mzs = &s.mz;
        let ints = &s.intensity;
        let mut acc = 0.0f64;
        let mut j = lower_bound(mzs, lo);
        let mut guard = 0usize;
        while j < mzs.len() {
            let v = mzs[j];
            if v > hi {
                break;
            }
            acc += ints[j] as f64;
            j += 1;
            guard += 1;
            if guard > 5_000_000 {
                panic!(
                    "[panic] compute_eic_for_mz inner loop too long at rt index {}",
                    i
                );
            }
        }
        y[i] = acc;
    }
    y
}

pub fn collect_ms1_scans(mzml: &MzML, time_window: FromTo) -> (Vec<f64>, Vec<CentroidScan>) {
    let mut scans = Vec::new();
    let mut total_points: usize = 0;
    let mut dropped_points: usize = 0;

    if let Some(run) = &mzml.run {
        for s in &run.spectra {
            let is_ms1 = matches!(s.ms_level, Some(1));
            let ok_rt = matches!(s.retention_time, Some(rt) if rt >= time_window.from && rt <= time_window.to);
            let has_arrays = s.mz_array.as_ref().map(|v| !v.is_empty()).unwrap_or(false)
                && s.intensity_array
                    .as_ref()
                    .map(|v| !v.is_empty())
                    .unwrap_or(false);
            if is_ms1 && ok_rt && has_arrays {
                let rt = s.retention_time.unwrap_or_default();
                let mzs_src = s.mz_array.clone().unwrap_or_default();
                let ints_src = s.intensity_array.clone().unwrap_or_default();
                let len = mzs_src.len().min(ints_src.len());
                let mut mzs = Vec::with_capacity(len);
                let mut ints = Vec::with_capacity(len);
                for i in 0..len {
                    let m = mzs_src[i];
                    let it = ints_src[i];
                    total_points += 1;
                    if m.is_finite() && it.is_finite() {
                        mzs.push(m);
                        ints.push(it);
                    } else {
                        dropped_points += 1;
                    }
                }
                if !mzs.is_empty() {
                    scans.push(CentroidScan {
                        rt,
                        mz: Arc::from(mzs),
                        intensity: Arc::from(ints),
                    });
                }
            }
        }
    }
    if dropped_points > 0 {
        eprintln!(
            "[warn] dropped {} non-finite points out of {}",
            dropped_points, total_points
        );
    }
    if scans.is_empty() {
        eprintln!("[warn] no MS1 scans found in time window");
    }
    scans.sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap_or(Ordering::Equal));
    let rt = scans.iter().map(|s| s.rt).collect::<Vec<_>>();
    (rt, scans)
}

fn max_in_range(rt: &[f64], y: &[f64], from_rt: f64, to_rt: f64) -> f64 {
    let i0 = lower_bound(rt, from_rt);
    let mut i1 = upper_bound(rt, to_rt);
    if i0 >= y.len() {
        return 0.0;
    }
    if i1 > y.len() {
        i1 = y.len();
    }
    if i1 <= i0 {
        return 0.0;
    }
    let mut m = y[i0];
    let mut i = i0 + 1;
    while i < i1 {
        let v = y[i];
        if v > m {
            m = v;
        }
        i += 1;
    }
    m
}

pub fn with_eic_apex_intensity(rt: &[f64], y: &[f64], mut p: Peak) -> Peak {
    let a = max_in_range(rt, y, p.from, p.to);
    if a.is_finite() && a > 0.0 {
        p.intensity = a;
    }
    p
}

#[inline]
pub fn lower_bound(a: &[f64], x: f64) -> usize {
    let mut lo = 0usize;
    let mut hi = a.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] < x {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

#[inline]
pub fn upper_bound(a: &[f64], x: f64) -> usize {
    let mut lo = 0usize;
    let mut hi = a.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] <= x {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}
