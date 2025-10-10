use rayon::{ThreadPoolBuilder, prelude::*};

use crate::utilities::{
    find_peaks::FindPeaksOptions,
    get_peak::get_peak,
    parse::parse_mzml::MzML,
    structs::{ChromRoi, DataXY, Roi},
};

pub fn get_peaks_from_chrom(
    mzml: &MzML,
    items: &[ChromRoi],
    options: Option<FindPeaksOptions>,
    cores: usize,
) -> Option<Vec<(usize, String, f64, f64, f64, f64, f64, f64)>> {
    let chroms = &mzml.run.as_ref()?.chromatograms;
    let f = |roi: &ChromRoi| {
        if roi.window <= 0.0 || !roi.rt.is_finite() {
            return (roi.idx, roi.id.clone(), roi.rt, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        let i = roi.idx;
        if i >= chroms.len() {
            return (i, roi.id.clone(), roi.rt, 0.0, 0.0, 0.0, 0.0, 0.0);
        }
        let ch = &chroms[i];
        let (x, y) = match (&ch.time_array, &ch.intensity_array) {
            (Some(t), Some(ints)) => (t.iter().copied().collect(), ints.iter().copied().collect()),
            _ => (Vec::new(), Vec::new()),
        };
        compute_one(ch.index, &ch.id, x, y, roi, &options)
    };
    if cores <= 1 || items.len() < 2 {
        Some(items.iter().map(f).collect())
    } else {
        let pool = ThreadPoolBuilder::new().num_threads(cores).build().ok()?;
        Some(pool.install(|| items.par_iter().map(f).collect()))
    }
}

#[inline]
fn compute_one(
    ch_index: usize,
    ch_id: &str,
    x: Vec<f64>,
    y: Vec<f64>,
    roi: &ChromRoi,
    options: &Option<FindPeaksOptions>,
) -> (usize, String, f64, f64, f64, f64, f64, f64) {
    let mut rt = 0.0_f64;
    let mut from = 0.0_f64;
    let mut to = 0.0_f64;
    let mut intensity_val = 0.0_f64;
    let mut integral = 0.0_f64;
    if x.len() >= 3 && x.len() == y.len() {
        if let Some(p) = get_peak(
            &DataXY { x, y },
            Roi {
                rt: roi.rt,
                window: roi.window,
            },
            options.clone(),
        ) {
            rt = p.rt;
            from = p.from;
            to = p.to;
            intensity_val = p.intensity;
            integral = p.integral;
        }
    }
    (
        ch_index,
        ch_id.trim_end().to_string(),
        roi.rt,
        rt,
        from,
        to,
        intensity_val,
        integral,
    )
}
