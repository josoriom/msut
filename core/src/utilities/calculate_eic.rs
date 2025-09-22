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
            ms1_idx.push(i);
            times.push(s.retention_time.unwrap_or(f64::NAN));
        }
    }
    if ms1_idx.is_empty() {
        return Err("no MS1 spectra");
    }
    if !times.iter().any(|t| t.is_finite()) {
        return Err("no valid retention_time");
    }

    let from_idx = lower_bound(&times, from_to.from);
    let to_idx = lower_bound(&times, from_to.to);
    let (start, end) = if from_idx <= to_idx {
        (from_idx, to_idx)
    } else {
        (to_idx, from_idx)
    };
    let n = end.saturating_sub(start);
    if n == 0 {
        return Ok(Eic {
            x: Vec::new(),
            y: Vec::new(),
        });
    }

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

        for &(min_mz, max_mz) in &mass_windows {
            let lo = lower_bound(mz, min_mz);
            let hi = upper_bound(mz, max_mz);
            if lo < hi {
                y[j] = intens[hi - 1];
            }
        }
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
    let mut lo = 0;
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
    let mut lo = 0;
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
