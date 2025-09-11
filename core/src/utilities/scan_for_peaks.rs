use crate::utilities::sgg::{SggOptions, sgg};
use crate::utilities::structs::DataXY;

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

pub fn scan_for_peaks(data: &DataXY, options: Option<ScanPeaksOptions>) -> Vec<f64> {
    let n = data.x.len();
    if n < 3 || n != data.y.len() {
        return Vec::new();
    }

    let opts = options.unwrap_or_default();
    let ws = opts.window_size.max(3);
    let s0 = SggOptions {
        window_size: ws,
        derivative: 0,
        polynomial: 3,
    };
    let s1 = SggOptions {
        window_size: ws,
        derivative: 1,
        polynomial: 3,
    };
    let ys_sm: Vec<f32> = sgg(&data.y, &data.x, s0);
    let dy: Vec<f32> = sgg(&data.y, &data.x, s1);
    let eps: f32 = opts.epsilon as f32;
    if dy.iter().all(|&v| v.abs() <= eps) {
        return Vec::new();
    }

    let min_sep = min_separation(&data.x, ws);

    let mut cands: Vec<(f64, f32, usize)> = Vec::with_capacity(n / 3);
    for k in 0..(n - 1) {
        let a = sign_eps(dy[k], eps);
        let b = sign_eps(dy[k + 1], eps);
        if (a > 0 && b <= 0) || (a >= 0 && b < 0) {
            let xp = refine_zero_cross(data.x[k], data.x[k + 1], dy[k], dy[k + 1]);
            let i = if (xp - data.x[k]).abs() <= (data.x[k + 1] - xp).abs() {
                k
            } else {
                k + 1
            };
            cands.push((xp, ys_sm[i], i));
        }
    }

    let mut i = 0usize;
    while i < n {
        if dy[i].abs() <= eps {
            let a = i;
            while i + 1 < n && dy[i + 1].abs() <= eps {
                i += 1
            }
            let b = i;
            if !(a == 0 && b + 1 >= n) {
                let left_ok = a == 0 || ys_sm[a] >= ys_sm[a - 1];
                let right_ok = b + 1 >= n || ys_sm[b] >= ys_sm[b + 1];
                if left_ok && right_ok {
                    let mut im = a;
                    let mut ym = ys_sm[a];
                    for j in (a + 1)..=b {
                        let yj = ys_sm[j];
                        if yj > ym {
                            ym = yj;
                            im = j;
                        }
                    }
                    let xp = if im > 0 && im + 1 < n {
                        quad_vertex(
                            data.x[im - 1],
                            data.x[im],
                            data.x[im + 1],
                            ys_sm[im - 1] as f64,
                            ys_sm[im] as f64,
                            ys_sm[im + 1] as f64,
                        )
                    } else {
                        data.x[im]
                    };
                    cands.push((xp, ym, im));
                }
            }
        }
        i += 1;
    }

    if cands.is_empty() {
        return Vec::new();
    }
    cands.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut out: Vec<f64> = Vec::with_capacity(cands.len());
    let mut last_x = f64::NEG_INFINITY;
    let mut last_i = usize::MAX;
    let mut last_y = f32::MIN;

    for (x, y, idx) in cands {
        if last_i == usize::MAX || x - last_x >= min_sep {
            out.push(x);
            last_x = x;
            last_i = idx;
            last_y = y;
        } else {
            let l = last_i.min(idx);
            let r = last_i.max(idx);
            let mut valley = f32::INFINITY;
            for j in l..=r {
                let v = ys_sm[j];
                if v < valley {
                    valley = v;
                }
            }
            let separated = valley <= 0.8 * last_y.min(y);
            if separated {
                out.push(x);
                last_x = x;
                last_i = idx;
                last_y = y;
            } else if y > last_y {
                if let Some(t) = out.last_mut() {
                    *t = x;
                }
                last_x = x;
                last_i = idx;
                last_y = y;
            }
        }
    }
    out
}

#[inline(always)]
fn sign_eps(v: f32, eps: f32) -> i8 {
    if v > eps {
        1
    } else if v < -eps {
        -1
    } else {
        0
    }
}

#[inline(always)]
fn refine_zero_cross(x0: f64, x1: f64, d0: f32, d1: f32) -> f64 {
    let denom = d0 - d1;
    if denom.abs() > f32::EPSILON {
        x0 + (x1 - x0) * (d0 / denom) as f64
    } else {
        0.5 * (x0 + x1)
    }
}

#[inline(always)]
fn quad_vertex(xm1: f64, x0: f64, xp1: f64, ym1: f64, y0: f64, yp1: f64) -> f64 {
    let a0 = xm1 - x0;
    let a1 = xp1 - x0;
    let dy0 = ym1 - y0;
    let dy1 = yp1 - y0;
    let denom = a0 * a1 * (a0 - a1);
    if denom == 0.0 {
        return x0;
    }
    let a = (dy0 * a1 - dy1 * a0) / denom;
    if a == 0.0 {
        return x0;
    }
    let b = (dy1 * a0 * a0 - dy0 * a1 * a1) / denom;
    x0 - b / (2.0 * a)
}

#[inline(always)]
fn min_separation(x: &[f64], window_size: usize) -> f64 {
    let n = x.len();
    let dx_avg = ((x[n - 1] - x[0]).abs()) / ((n as f64) - 1.0);
    let mut dx_min = f64::INFINITY;
    for w in x.windows(2) {
        let d = w[1] - w[0];
        if d > 0.0 && d < dx_min {
            dx_min = d;
        }
    }
    let dx_min = if dx_min.is_finite() {
        dx_min
    } else {
        dx_avg.max(f64::EPSILON)
    };
    let sep_ws = 0.25 * (window_size as f64) * dx_avg;
    let sep_floor = 1.5 * dx_min;
    sep_ws.max(sep_floor)
}
