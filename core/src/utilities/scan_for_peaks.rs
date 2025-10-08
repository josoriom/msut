use crate::utilities::find_noise_level::find_noise_level;
use crate::utilities::sgg::{SggOptions, sgg};
use crate::utilities::structs::DataXY;
use crate::utilities::{
    Point, closest_index, kmeans, mean_step, min_positive_step, min_sep, odd_in_range, quad_peak,
};

#[derive(Clone, Copy, Debug)]
pub struct ScanPeaksOptions {
    pub epsilon: f64,
    pub window_size: usize,
}
impl Default for ScanPeaksOptions {
    fn default() -> Self {
        Self {
            epsilon: 1e-5,
            window_size: 11,
        }
    }
}

const DEFAULT_WINDOW_SIZES: &[usize] = &[5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33];

pub fn scan_for_peaks(data: &DataXY, options: Option<ScanPeaksOptions>) -> Vec<f64> {
    scan_for_peaks_across_windows(data, options, None)
}

pub fn scan_for_peaks_across_windows(
    data: &DataXY,
    options: Option<ScanPeaksOptions>,
    window_sizes: Option<&[usize]>,
) -> Vec<f64> {
    let n = data.x.len();
    if n < 3 || n != data.y.len() {
        return Vec::new();
    }

    let opts = options.unwrap_or_default();
    let epsilon = opts.epsilon as f64;

    let mut windows = Vec::new();
    for &w in window_sizes.unwrap_or(DEFAULT_WINDOW_SIZES) {
        if let Some(ow) = odd_in_range(w, n) {
            windows.push(ow);
        }
    }
    if windows.is_empty() {
        return Vec::new();
    }

    let mut xs = Vec::<f64>::new();
    let mut hs = Vec::<f64>::new();
    for w in windows.iter().copied() {
        for (rt, h) in scan_one_window(data, w, epsilon) {
            xs.push(rt);
            hs.push(h);
        }
    }
    if xs.is_empty() {
        return Vec::new();
    }

    let step = mean_step(&data.x);
    let mut rt_tol = (8.0 * step).max(0.02);
    let max_w = *windows.iter().max().unwrap() as f64;
    let bump = 0.10 * max_w * step;
    if bump > rt_tol {
        rt_tol = bump;
    }

    let (centers, votes, best_h) = group_peak_positions(&xs, &hs, rt_tol);

    let need = if windows.len() <= 3 {
        2
    } else {
        ((windows.len() as f64) * 0.30).ceil() as usize
    };
    let mut h_max = 0.0f64;
    for &h in &best_h {
        if h > h_max {
            h_max = h;
        }
    }
    let strong_gate = 0.60 * h_max;

    let mut keep = Vec::<(f64, f64)>::new();
    for i in 0..centers.len() {
        if votes[i] >= need || best_h[i] >= strong_gate {
            keep.push((centers[i], best_h[i]));
        }
    }
    if keep.is_empty() {
        return Vec::new();
    }
    keep.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let noise = find_noise_level(&data.y);
    let w_largest = *windows.iter().max().unwrap();
    let y64: Vec<f64> = data.y.iter().map(|&v| v as f64).collect();

    let y_smooth = if w_largest <= n {
        sgg(
            &y64,
            &data.x,
            SggOptions {
                window_size: w_largest,
                derivative: 0,
                polynomial: 3,
            },
        )
    } else {
        y64.clone()
    };

    let mut merged = Vec::<(f64, f64)>::new();
    let mut i = 0usize;
    while i < keep.len() {
        let (xa, ha) = keep[i];
        if i + 1 >= keep.len() {
            merged.push((xa, ha));
            break;
        }
        let (xb, hb) = keep[i + 1];

        let ia = closest_index(&data.x, xa);
        let ib = closest_index(&data.x, xb);
        let l = ia.min(ib);
        let r = ia.max(ib);

        let mut valley = f64::INFINITY;
        for j in l..=r {
            if y_smooth[j] < valley {
                valley = y_smooth[j];
            }
        }

        let h_top = if ha > hb { ha } else { hb } as f64;
        let rel_drop = if h_top > 0.0 {
            (h_top - valley) / h_top
        } else {
            0.0
        };
        let valley_gate = noise.max(0.92_f64 * (h_top as f64));

        if (valley as f64) > valley_gate || rel_drop < 0.08_f64 {
            merged.push(if ha >= hb { (xa, ha) } else { (xb, hb) });
            i += 2;
        } else {
            merged.push((xa, ha));
            i += 1;
        }
    }
    if merged.is_empty() {
        return Vec::new();
    }

    let sep = min_sep(&data.x, 5);
    let mut out = Vec::<f64>::new();
    let mut last_x = f64::NEG_INFINITY;
    let mut last_h = f64::NEG_INFINITY;
    for (x, h) in merged {
        if !last_x.is_finite() || x - last_x >= sep {
            out.push(x);
            last_x = x;
            last_h = h;
        } else if h > last_h {
            if let Some(t) = out.last_mut() {
                *t = x;
            }
            last_x = x;
            last_h = h;
        }
    }
    out
}

fn scan_one_window(data: &DataXY, window: usize, epsilon: f64) -> Vec<(f64, f64)> {
    let n = data.x.len();
    if n < 3 {
        return Vec::new();
    }
    let w = window.max(3);

    let sm0 = SggOptions {
        window_size: w,
        derivative: 0,
        polynomial: 3,
    };
    let sm1 = SggOptions {
        window_size: w,
        derivative: 1,
        polynomial: 3,
    };
    let sm2 = SggOptions {
        window_size: w,
        derivative: 2,
        polynomial: 3,
    };

    let y_f64: Vec<f64> = data.y.iter().copied().map(f64::from).collect();

    let y0 = sgg(&y_f64, &data.x, sm0);
    let y1 = sgg(&y_f64, &data.x, sm1);
    let y2 = sgg(&y_f64, &data.x, sm2);

    let mut cand = Vec::<(f64, f64, usize)>::with_capacity(n / 3);

    for k in 0..(n - 1) {
        let a = if y1[k] > epsilon {
            1
        } else if y1[k] < -epsilon {
            -1
        } else {
            0
        };
        let b = if y1[k + 1] > epsilon {
            1
        } else if y1[k + 1] < -epsilon {
            -1
        } else {
            0
        };
        if (a > 0 && b <= 0) || (a >= 0 && b < 0) {
            let denom = y1[k] - y1[k + 1];
            let xp = if denom.abs() > f64::EPSILON {
                data.x[k] + (data.x[k + 1] - data.x[k]) * (y1[k] / denom) as f64
            } else {
                0.5 * (data.x[k] + data.x[k + 1])
            };
            let i = if (xp - data.x[k]).abs() <= (data.x[k + 1] - xp).abs() {
                k
            } else {
                k + 1
            };
            if y2[i] < 0.0 {
                cand.push((xp, y0[i], i));
            }
        }
    }

    let mut i = 0usize;
    while i < n {
        if y1[i].abs() <= epsilon {
            let a = i;
            while i + 1 < n && y1[i + 1].abs() <= epsilon {
                i += 1;
            }
            let b = i;
            let left_ok = a == 0 || y0[a] >= y0[a - 1];
            let right_ok = b + 1 >= n || y0[b] >= y0[b + 1];
            if left_ok && right_ok {
                let mut im = a;
                let mut ym = y0[a];
                for j in (a + 1)..=b {
                    if y0[j] > ym {
                        ym = y0[j];
                        im = j;
                    }
                }
                if y2[im] < 0.0 {
                    let xp = if im > 0 && im + 1 < n {
                        quad_peak(&data.x, &y0, im)
                    } else {
                        data.x[im]
                    };
                    cand.push((xp, ym, im));
                }
            }
        }
        i += 1;
    }

    if cand.is_empty() {
        return Vec::new();
    }
    cand.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let mut out = Vec::<(f64, f64)>::with_capacity(cand.len());
    let min_step = min_positive_step(&data.x)
        .unwrap_or_else(|| mean_step(&data.x))
        .max(f64::EPSILON);
    let min_sep = (1.2 * min_step).max(0.01);

    let mut last_x = f64::NEG_INFINITY;
    let mut last_y = f64::NEG_INFINITY;
    for (x, y, _) in cand {
        if !last_x.is_finite() || x - last_x >= min_sep {
            out.push((x, y));
            last_x = x;
            last_y = y;
        } else if y > last_y {
            if let Some(t) = out.last_mut() {
                *t = (x, y);
            }
            last_x = x;
            last_y = y;
        }
    }
    out
}

fn group_peak_positions(xs: &[f64], hs: &[f64], tol: f64) -> (Vec<f64>, Vec<usize>, Vec<f64>) {
    if xs.is_empty() || xs.len() != hs.len() {
        return (Vec::new(), Vec::new(), Vec::new());
    }

    let mut xmin = xs[0];
    let mut xmax = xs[0];
    for &x in xs {
        if x < xmin {
            xmin = x;
        }
        if x > xmax {
            xmax = x;
        }
    }
    let range = (xmax - xmin).abs();
    let mut k = if tol > 0.0 && range > 0.0 {
        (range / tol).ceil() as usize
    } else {
        1
    };
    if k == 0 {
        k = 1;
    }
    if k > xs.len() {
        k = xs.len();
    }

    let mut order: Vec<usize> = (0..xs.len()).collect();
    order.sort_by(|a, b| xs[*a].partial_cmp(&xs[*b]).unwrap());

    let mut init = Vec::<Point>::with_capacity(k);
    if k == 1 {
        init.push(vec![xs[order[xs.len() / 2]]]);
    } else {
        for s in 0..k {
            let pos = s * (xs.len() - 1) / (k - 1);
            init.push(vec![xs[order[pos]]]);
        }
    }

    let pts: Vec<Point> = xs.iter().map(|&x| vec![x]).collect();
    let cents = kmeans(&pts, init);

    let mut centers = Vec::<f64>::with_capacity(cents.len());
    for c in &cents {
        centers.push(c[0]);
    }
    let kfin = centers.len();

    let mut votes = vec![0usize; kfin];
    let mut best_h = vec![f64::NEG_INFINITY; kfin];

    for i in 0..xs.len() {
        let x = xs[i];
        let h = hs[i];
        let mut idx = 0usize;
        let mut best = (x - centers[0]).abs();
        for j in 1..kfin {
            let d = (x - centers[j]).abs();
            if d < best {
                best = d;
                idx = j;
            }
        }
        votes[idx] += 1;
        if h > best_h[idx] {
            best_h[idx] = h;
        }
    }

    (centers, votes, best_h)
}
