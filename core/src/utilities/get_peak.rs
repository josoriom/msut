use crate::utilities::find_peaks::{FindPeaksOptions, find_peaks};
use crate::utilities::structs::{DataXY, Peak, Roi};

pub fn get_peak(data: &DataXY, roi: Roi, options: Option<FindPeaksOptions>) -> Option<Peak> {
    let peaks = find_peaks(data, options);
    if peaks.is_empty() {
        return None;
    }

    let target = roi.rt;
    let w = if roi.window.is_finite() && roi.window > 0.0 {
        roi.window
    } else {
        0.0
    };

    let mut best: Option<&Peak> = None;

    if w > 0.0 {
        let lo = target - w;
        let hi = target + w;
        for p in &peaks {
            if p.rt.is_finite() && p.rt >= lo && p.rt <= hi {
                match best {
                    None => best = Some(p),
                    Some(b) => {
                        let db = (b.rt - target).abs();
                        let dp = (p.rt - target).abs();
                        if dp < db || ((dp - db).abs() <= f64::EPSILON && p.intensity > b.intensity)
                        {
                            best = Some(p);
                        }
                    }
                }
            }
        }
    }

    if best.is_none() {
        for p in &peaks {
            match best {
                None => best = Some(p),
                Some(b) => {
                    let db = (b.rt - target).abs();
                    let dp = (p.rt - target).abs();
                    if dp < db || ((dp - db).abs() <= f64::EPSILON && p.intensity > b.intensity) {
                        best = Some(p);
                    }
                }
            }
        }
    }

    best.cloned()
}
