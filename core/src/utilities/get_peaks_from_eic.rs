use rayon::{ThreadPoolBuilder, prelude::*};

use crate::utilities::{
    EicOptions, calculate_eic_from_mzml,
    find_peaks::FindPeaksOptions,
    get_peak::get_peak,
    parse::{decode, parse_mzml::MzML},
    structs::{DataXY, EicRoi, FromTo, Peak, Roi},
};

pub fn get_peaks_from_eic(
    bytes: &[u8],
    from_to: FromTo,
    rois: &[EicRoi],
    options: Option<FindPeaksOptions>,
    cores: usize,
) -> Option<Vec<(String, f64, f64, Peak)>> {
    let mzml = decode(bytes).ok()?;
    if cores <= 1 || rois.len() < 2 {
        let mut out: Vec<(String, f64, f64, Peak)> = Vec::with_capacity(rois.len());
        for roi in rois {
            out.push(compute_one(&mzml, from_to, roi, &options));
        }
        return Some(out);
    }
    let pool = ThreadPoolBuilder::new().num_threads(cores).build().ok()?;
    Some(pool.install(|| {
        rois.par_iter()
            .map(|roi| compute_one(&mzml, from_to, roi, &options))
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
    let mz_str = roi.mz.to_string();
    let eic = match calculate_eic_from_mzml(
        mzml,
        &mz_str,
        from_to,
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
