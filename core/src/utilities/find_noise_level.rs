use std::borrow::Cow;
use std::f64;

use crate::utilities::kmeans;

pub fn find_noise_level<'a, Y>(y: Y) -> f64
where
    Y: IntoF64Slice<'a>,
{
    let y = y.into_f64_slice();
    let points = scan_points(y.as_ref());
    if points.is_empty() {
        return 0.0;
    }

    let params = ClusterParams {
        log_n: true,
        log_i: true,
        w_n: 3.0,
        w_i: 1.0,
    };
    let labels = cluster_intensities(&points, params);

    let mut max_in_noise = f64::NEG_INFINITY;
    for (i, &(_n, intensity)) in points.iter().enumerate() {
        if labels[i] == 0 && intensity.is_finite() && intensity > max_in_noise {
            max_in_noise = intensity;
        }
    }

    if max_in_noise.is_finite() {
        max_in_noise
    } else {
        0.0
    }
}

pub trait IntoF64Slice<'a> {
    fn into_f64_slice(self) -> Cow<'a, [f64]>;
}

impl<'a> IntoF64Slice<'a> for &'a [f64] {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        Cow::Borrowed(self)
    }
}
impl<'a> IntoF64Slice<'a> for &'a Vec<f64> {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        Cow::Borrowed(self.as_slice())
    }
}
impl<'a> IntoF64Slice<'a> for Vec<f64> {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        Cow::Owned(self)
    }
}

impl<'a> IntoF64Slice<'a> for &'a [f32] {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        let mut v = Vec::with_capacity(self.len());
        for &x in self {
            v.push(x as f64);
        }
        Cow::Owned(v)
    }
}
impl<'a> IntoF64Slice<'a> for &'a Vec<f32> {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        let mut v = Vec::with_capacity(self.len());
        for &x in self {
            v.push(x as f64);
        }
        Cow::Owned(v)
    }
}
impl<'a> IntoF64Slice<'a> for Vec<f32> {
    fn into_f64_slice(self) -> Cow<'a, [f64]> {
        let mut v = Vec::with_capacity(self.len());
        for x in self {
            v.push(x as f64);
        }
        Cow::Owned(v)
    }
}

#[derive(Clone, Copy)]
struct ClusterParams {
    log_n: bool,
    log_i: bool,
    w_n: f64,
    w_i: f64,
}

fn scan_points(y: &[f64]) -> Vec<(usize, f64)> {
    let mut runs: Vec<(usize, f64)> = Vec::new();
    if y.is_empty() {
        return runs;
    }
    let min = x_min_value(y);

    let mut in_run = false;
    let mut count: usize = 0;
    let mut max_val = f64::NEG_INFINITY;

    for &v in y {
        let is_non_zero = v.is_finite() && v.abs() > min;
        if is_non_zero {
            if !in_run {
                in_run = true;
                count = 1;
                max_val = v;
            } else {
                count += 1;
                if v > max_val {
                    max_val = v;
                }
            }
        } else {
            if in_run {
                runs.push((count, max_val));
                in_run = false;
                count = 0;
                max_val = f64::NEG_INFINITY;
            }
        }
    }

    if in_run {
        runs.push((count, max_val));
    }

    runs
}

fn cluster_intensities(points: &[(usize, f64)], params: ClusterParams) -> Vec<usize> {
    if points.is_empty() {
        return Vec::new();
    }

    let mut tx: Vec<f64> = Vec::with_capacity(points.len());
    let mut ty: Vec<f64> = Vec::with_capacity(points.len());
    for &(n, i) in points.iter() {
        let a = if params.log_n {
            log1p_safe(n as f64)
        } else {
            n as f64
        };
        let b = if params.log_i { log1p_safe(i) } else { i };
        tx.push(a);
        ty.push(b);
    }

    let mx = x_mean(&tx);
    let my = x_mean(&ty);
    let sx = x_stddev(&tx);
    let sy = x_stddev(&ty);
    let sx = if sx == 0.0 { 1.0 } else { sx };
    let sy = if sy == 0.0 { 1.0 } else { sy };

    let wx = params.w_n.sqrt();
    let wy = params.w_i.sqrt();

    let mut weighted: Vec<[f64; 2]> = Vec::with_capacity(points.len());
    for i in 0..points.len() {
        let a = (tx[i] - mx) / sx;
        let b = (ty[i] - my) / sy;
        weighted.push([a * wx, b * wy]);
    }

    let seeds = farthest_seeds(&weighted);

    let mut km_points: Vec<Vec<f64>> = Vec::with_capacity(weighted.len());
    for p in &weighted {
        km_points.push(vec![p[0], p[1]]);
    }
    let centroids = kmeans(&km_points, seeds);

    let mut centroids_sorted = centroids.clone();
    centroids_sorted.sort_by(|a, b| a[1].partial_cmp(&b[1]).unwrap());

    // println!("{centroids:#?}");
    // println!("{centroids_sorted:#?}");

    let mut result: Vec<usize> = Vec::with_capacity(weighted.len());
    for p in weighted.iter() {
        let mut idx = 0usize;
        let mut best = {
            let dx = p[0] - centroids_sorted[0][0];
            let dy = p[1] - centroids_sorted[0][1];
            (dx * dx + dy * dy).sqrt()
        };
        for j in 1..centroids_sorted.len() {
            let dx = p[0] - centroids_sorted[j][0];
            let dy = p[1] - centroids_sorted[j][1];
            let d = (dx * dx + dy * dy).sqrt();
            if d < best {
                best = d;
                idx = j;
            }
        }
        result.push(idx);
    }

    result
}

fn farthest_seeds(weighted: &[[f64; 2]]) -> Vec<Vec<f64>> {
    let mut min = weighted[0];
    let mut max = weighted[0];
    for &q in weighted.iter() {
        if q[1] < min[1] {
            min = q;
        }
        if q[1] > max[1] {
            max = q;
        }
    }
    vec![vec![min[0], min[1]], vec![max[0], max[1]]]
}

fn x_min_value(a: &[f64]) -> f64 {
    let mut m = f64::INFINITY;
    for &v in a {
        if v.is_finite() && v < m {
            m = v;
        }
    }
    if m.is_finite() { m } else { 0.0 }
}

fn x_mean(a: &[f64]) -> f64 {
    if a.is_empty() {
        return 0.0;
    }
    let mut s = 0.0;
    let mut n = 0usize;
    for &v in a {
        if v.is_finite() {
            s += v;
            n += 1;
        }
    }
    if n == 0 { 0.0 } else { s / n as f64 }
}

fn x_stddev(a: &[f64]) -> f64 {
    if a.is_empty() {
        return 0.0;
    }
    let mu = x_mean(a);
    let mut s = 0.0;
    let mut n = 0usize;
    for &v in a {
        if v.is_finite() {
            let d = v - mu;
            s += d * d;
            n += 1;
        }
    }
    if n == 0 { 0.0 } else { (s / n as f64).sqrt() }
}

fn log1p_safe(v: f64) -> f64 {
    (0.0_f64.max(v)).ln_1p()
}
