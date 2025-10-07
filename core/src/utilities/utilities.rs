#[inline]
pub fn is_finite_positive(x: f32) -> bool {
    x.is_finite() && x > 0.0
}

#[inline]
pub fn closest_index(xs: &[f64], v: f64) -> usize {
    if xs.is_empty() {
        return 0;
    }
    let (mut lo, mut hi) = (0usize, xs.len());
    // lower_bound(v)
    while lo < hi {
        let mid = (lo + hi) / 2;
        if xs[mid] < v {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo == 0 {
        0
    } else if lo >= xs.len() {
        xs.len() - 1
    } else if (v - xs[lo - 1]).abs() <= (xs[lo] - v).abs() {
        lo - 1
    } else {
        lo
    }
}

#[inline]
pub fn odd_in_range(w: usize, n: usize) -> Option<usize> {
    if n < 5 || w < 5 || w > n {
        return None;
    }
    let o = if w % 2 == 0 { w - 1 } else { w };
    if o >= 5 && o <= n { Some(o) } else { None }
}

#[inline]
pub fn mean_step(xs: &[f64]) -> f64 {
    let mut sum = 0.0;
    let mut c = 0usize;
    for w in xs.windows(2) {
        let d = (w[1] - w[0]).abs();
        if d.is_finite() && d > 0.0 {
            sum += d;
            c += 1;
        }
    }
    if c == 0 {
        f64::EPSILON.max(0.01)
    } else {
        (sum / c as f64).max(f64::EPSILON)
    }
}

#[inline]
pub fn min_positive_step(xs: &[f64]) -> Option<f64> {
    let mut m = f64::INFINITY;
    for w in xs.windows(2) {
        let d = w[1] - w[0];
        if d > 0.0 && d < m {
            m = d;
        }
    }
    if m.is_finite() { Some(m) } else { None }
}

#[inline]
pub fn min_sep(xs: &[f64], window_size: usize) -> f64 {
    let avg = mean_step(xs);
    let base = min_positive_step(xs).unwrap_or(avg).max(f64::EPSILON);
    let by_window = 0.15 * (window_size as f64) * avg;
    let floor = 1.2 * base;
    if by_window > floor { by_window } else { floor }
}

#[inline]
pub fn quad_peak(xs: &[f64], ys: &[f64], i: usize) -> f64 {
    debug_assert!(i > 0 && i + 1 < xs.len());
    let xm1 = xs[i - 1];
    let x0 = xs[i];
    let xp1 = xs[i + 1];
    let ym1 = ys[i - 1] as f64;
    let y0 = ys[i] as f64;
    let yp1 = ys[i + 1] as f64;

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
        x0
    } else {
        let b = (dy1 * a0 * a0 - dy0 * a1 * a1) / denom;
        x0 - b / (2.0 * a)
    }
}

pub fn xy_integration(x: &[f64], y: &[f32]) -> (f64, f64) {
    let n = x.len();
    if n == 0 || n != y.len() {
        return (0.0, f64::NEG_INFINITY);
    }
    if n == 1 {
        return (0.0, y[0] as f64);
    }
    let mut s = 0.0f64;
    let mut m = y[0];
    for i in 0..(n - 1) {
        let dx = x[i + 1] - x[i];
        let yi = y[i] as f64;
        let yj = y[i + 1] as f64;
        s += dx * (yi + yj) * 0.5;
        if y[i + 1] > m {
            m = y[i + 1];
        }
    }
    (s, m as f64)
}
