use rayon::{ThreadPoolBuilder, prelude::*};

use crate::utilities::{
    find_peaks::FindPeaksOptions,
    get_peak::get_peak,
    parse::{
        calculate_eic::{EicOptions, calculate_eic_from_mzml},
        decode,
        parse_mzml::MzML,
    },
    structs::{DataXY, EicRoi, FromTo, Peak, Roi},
};

pub fn get_peaks_from_eic(
    bytes: &[u8],
    from_to: FromTo,
    items: &[EicRoi],
    options: Option<FindPeaksOptions>,
    cores: usize,
) -> Option<Vec<(String, f64, f64, Peak)>> {
    let mzml = decode(bytes).ok()?;
    if cores <= 1 || items.len() < 2 {
        let mut out = Vec::with_capacity(items.len());
        for it in items {
            out.push(compute_one(&mzml, from_to, it, &options));
        }
        return Some(out);
    }
    let pool = ThreadPoolBuilder::new().num_threads(cores).build().ok()?;
    Some(pool.install(|| {
        items
            .par_iter()
            .map(|it| compute_one(&mzml, from_to, it, &options))
            .collect()
    }))
}

#[inline]
fn compute_one(
    mzml: &MzML,
    from_to: FromTo,
    roi: &EicRoi,
    options: &Option<FindPeaksOptions>,
) -> (String, f64, f64, Peak) {
    let l = roi.rt - from_to.from;
    let r = roi.rt + from_to.to;
    let (rt_lo, rt_hi) = if l <= r { (l, r) } else { (r, l) };
    let mz_str = roi.mz.to_string();
    let eic = match calculate_eic_from_mzml(
        mzml,
        &mz_str,
        (rt_lo, rt_hi),
        EicOptions {
            ppm_tolerance: 20.0,
            mz_tolerance: 0.005,
        },
    ) {
        Ok(v) => v,
        Err(_) => return (roi.id.clone(), roi.rt, roi.mz, Peak::default()),
    };
    if eic.x.len() < 3 || eic.x.len() != eic.y.len() {
        return (roi.id.clone(), roi.rt, roi.mz, Peak::default());
    }
    let pk = match get_peak(
        &DataXY { x: eic.x, y: eic.y },
        Roi {
            rt: roi.rt,
            window: roi.window,
        },
        options.clone(),
    ) {
        Some(p) => p,
        None => Peak::default(),
    };
    (roi.id.clone(), roi.rt, roi.mz, pk)
}
