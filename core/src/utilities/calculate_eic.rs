use crate::utilities::{
    parse::{decode::decode, parse_mzml::MzML},
    structs::FromTo,
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
    pub y: Vec<f32>,
}

pub fn calculate_eic_from_bin1(
    bin1: &[u8],
    target_masses: &str,
    from_to: FromTo,
    options: EicOptions,
) -> Result<Eic, &'static str> {
    let mzml = decode(bin1).map_err(|_| "decode BIN1 failed")?;
    calculate_eic_from_mzml(&mzml, target_masses, from_to, options)
}

pub fn calculate_eic_from_mzml(
    mzml: &MzML,
    target_masses: &str,
    from_to: FromTo,
    options: EicOptions,
) -> Result<Eic, &'static str> {
    let run = mzml.run.as_ref().ok_or("no run")?;

    let mut ms1_idx = Vec::with_capacity(run.spectra.len());
    let mut times = Vec::with_capacity(run.spectra.len());
    for (i, s) in run.spectra.iter().enumerate() {
        if s.ms_level.unwrap_or(0) == 1 {
            if let Some(rt) = s.retention_time {
                if rt.is_finite() {
                    ms1_idx.push(i);
                    times.push(rt);
                }
            }
        }
    }
    if ms1_idx.is_empty() {
        return Err("no MS1 spectra");
    }
    if times.is_empty() {
        return Err("no valid retention_time");
    }

    let start = lower_bound(&times, from_to.from);
    let end = upper_bound(&times, from_to.to);
    if start >= end {
        return Ok(Eic {
            x: Vec::new(),
            y: Vec::new(),
        });
    }
    let n = end - start;

    let masses = parse_target_masses(target_masses).ok_or("invalid masses")?;
    if masses.is_empty() {
        return Err("no valid masses");
    }
    let mass_windows: Vec<(f64, f64)> = masses
        .iter()
        .map(|&m| {
            let tol = ((options.ppm_tolerance / 1e6) * m).max(options.mz_tolerance);
            (m - tol, m + tol)
        })
        .collect();

    let mut x = Vec::with_capacity(n);
    let mut y = vec![0f32; n];

    for (j, &spec_i) in ms1_idx[start..end].iter().enumerate() {
        let s = &run.spectra[spec_i];
        x.push(s.retention_time.unwrap_or(f64::NAN));

        let (mz_opt, intens_opt) = (s.mz_array.as_deref(), s.intensity_array.as_deref());
        let (mz, intens) = match (mz_opt, intens_opt) {
            (Some(mz), Some(intens)) if !mz.is_empty() && mz.len() == intens.len() => (mz, intens),
            _ => continue,
        };

        let mut acc = 0.0f64;
        for &(min_mz, max_mz) in &mass_windows {
            let lo = lower_bound(mz, min_mz);
            let hi = upper_bound(mz, max_mz);
            if lo < hi {
                acc += intens[lo..hi]
                    .iter()
                    .copied()
                    .map(|v| v as f64)
                    .sum::<f64>();
            }
        }
        y[j] = acc as f32;
    }

    Ok(Eic { x, y })
}

fn parse_target_masses(s: &str) -> Option<Vec<f64>> {
    let mut out = Vec::new();
    for tok in s.split(|c: char| c.is_whitespace() || c == ',' || c == ';') {
        if tok.is_empty() {
            continue;
        }
        match tok.parse::<f64>() {
            Ok(v) if v.is_finite() => out.push(v),
            _ => return None,
        }
    }
    Some(out)
}

#[inline]
fn lower_bound(a: &[f64], x: f64) -> usize {
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
fn upper_bound(a: &[f64], x: f64) -> usize {
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
