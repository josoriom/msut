use crate::utilities::kmeans::{Point, kmeans};
use crate::utilities::utilities::is_finite_positive;

pub fn find_noise_level(intensities: &[f32]) -> f32 {
    let n = intensities.len();
    if n == 0 {
        return f32::INFINITY;
    }
    let est = noise_from_series(intensities);
    if est.is_finite() && est > 0.0 {
        est
    } else {
        f32::INFINITY
    }
}

fn noise_from_series(y: &[f32]) -> f32 {
    let n = y.len();
    if n < 128 {
        return f32::INFINITY;
    }
    let (w, s) = window_plan(n);
    let (bas, caps, spans) = window_low_quantiles(y, w, s, 0.20, 0.30);
    if bas.len() < 4 {
        let mut lows = Vec::<f32>::new();
        for (i, (a, b)) in spans.iter().enumerate() {
            let cap = caps[i];
            for &v in &y[*a..*b] {
                if is_finite_positive(v) && v <= cap {
                    lows.push(v);
                }
            }
        }
        if lows.is_empty() {
            let mut v: Vec<f32> = y
                .iter()
                .copied()
                .filter(|x| is_finite_positive(*x))
                .collect();
            if v.is_empty() {
                return f32::INFINITY;
            }
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            return vec_quantile_sorted(&v, 0.80);
        }
        lows.sort_by(|a, b| a.partial_cmp(b).unwrap());
        return vec_quantile_sorted(&lows, 0.995);
    }

    let mut logs = Vec::<f32>::with_capacity(bas.len());
    for &b in &bas {
        logs.push(b.log10());
    }

    let sep = kmeans_1d_logs(&logs, 30);
    let thr_log = if let Some((m0, m1)) = sep {
        let d = (m1 - m0).abs();
        if d < 0.12 {
            vec_quantile_sorted(
                &{
                    let mut t = logs.clone();
                    t.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    t
                },
                0.75,
            )
        } else {
            0.5 * (m0 + m1)
        }
    } else {
        vec_quantile_sorted(
            &{
                let mut t = logs.clone();
                t.sort_by(|a, b| a.partial_cmp(b).unwrap());
                t
            },
            0.75,
        )
    };

    let thr = 10f32.powf(thr_log);

    let mut sel = vec![false; bas.len()];
    let mut any = 0usize;
    for i in 0..bas.len() {
        if bas[i] >= thr {
            sel[i] = true;
            any += 1;
        }
    }

    let mut pool = if any > 0 {
        pool_from_windows(y, &spans, &caps, &sel)
    } else {
        Vec::new()
    };

    if pool.len() < y.len() / 200 {
        let mut idx: Vec<usize> = (0..bas.len()).collect();
        idx.sort_by(|&i, &j| bas[i].partial_cmp(&bas[j]).unwrap());
        let keep = (bas.len() as f32 * 0.25).ceil() as usize;
        for k in 0..idx.len() {
            sel[idx[k]] = false;
        }
        for k in (idx.len().saturating_sub(keep))..idx.len() {
            sel[idx[k]] = true;
        }
        pool = pool_from_windows(y, &spans, &caps, &sel);
    }

    if pool.is_empty() {
        let mut lows = Vec::<f32>::new();
        for (i, (a, b)) in spans.iter().enumerate() {
            let cap = caps[i];
            for &v in &y[*a..*b] {
                if is_finite_positive(v) && v <= cap {
                    lows.push(v);
                }
            }
        }
        if lows.is_empty() {
            let mut v: Vec<f32> = y
                .iter()
                .copied()
                .filter(|x| is_finite_positive(*x))
                .collect();
            if v.is_empty() {
                return f32::INFINITY;
            }
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            return vec_quantile_sorted(&v, 0.80);
        }
        lows.sort_by(|a, b| a.partial_cmp(b).unwrap());
        return vec_quantile_sorted(&lows, 0.995);
    }

    pool.sort_by(|a, b| a.partial_cmp(b).unwrap());
    vec_quantile_sorted(&pool, 0.995)
}

fn kmeans_1d_logs(x: &[f32], _iters: usize) -> Option<(f32, f32)> {
    if x.len() < 2 {
        return None;
    }

    let mut minv = f32::INFINITY;
    let mut maxv = f32::NEG_INFINITY;
    for &v in x {
        if v < minv {
            minv = v;
        }
        if v > maxv {
            maxv = v;
        }
    }
    if !minv.is_finite() || !maxv.is_finite() || maxv <= minv {
        return None;
    }

    let points: Vec<Point> = x.iter().map(|&v| vec![v as f64]).collect();
    let centroids_init: Vec<Point> = vec![vec![minv as f64], vec![maxv as f64]];

    let res = kmeans(&points, centroids_init);
    if res.len() < 2 {
        return None;
    }

    let c0 = res[0][0] as f32;
    let c1 = res[1][0] as f32;
    if c0.is_finite() && c1.is_finite() {
        if c0 < c1 {
            Some((c0, c1))
        } else {
            Some((c1, c0))
        }
    } else {
        None
    }
}

#[inline]
fn vec_quantile_sorted(v: &[f32], q: f32) -> f32 {
    if v.is_empty() {
        return f32::NAN;
    }
    if v.len() == 1 {
        return v[0];
    }
    let r = (q.clamp(0.0, 1.0) * (v.len() as f32 - 1.0)) as f32;
    let i = r.floor() as usize;
    let j = r.ceil() as usize;
    if i == j {
        v[i]
    } else {
        let w = r - i as f32;
        v[i] * (1.0 - w) + v[j] * w
    }
}

fn window_plan(n: usize) -> (usize, usize) {
    let mut w = (n / 64).max(256).min(4096);
    if w > n {
        w = n;
    }
    let s = (w / 4).max(64);
    (w, s)
}

fn window_low_quantiles(
    y: &[f32],
    w: usize,
    s: usize,
    q_low: f32,
    q_cap: f32,
) -> (Vec<f32>, Vec<f32>, Vec<(usize, usize)>) {
    let mut bas = Vec::<f32>::new();
    let mut caps = Vec::<f32>::new();
    let mut spans = Vec::<(usize, usize)>::new();
    let n = y.len();
    if n == 0 {
        return (bas, caps, spans);
    }
    let mut i = 0usize;
    while i < n {
        let a = i;
        let b = (i + w).min(n);
        let mut v: Vec<f32> = y[a..b]
            .iter()
            .copied()
            .filter(|x| is_finite_positive(*x))
            .collect();
        if v.len() >= 32 {
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let p20 = vec_quantile_sorted(&v, q_low);
            let p30 = vec_quantile_sorted(&v, q_cap).max(p20);
            if p20.is_finite() && p20 > 0.0 {
                bas.push(p20);
                caps.push(p30);
                spans.push((a, b));
            }
        }
        if b == n {
            break;
        }
        i = a + s;
    }
    (bas, caps, spans)
}

fn pool_from_windows(
    y: &[f32],
    spans: &[(usize, usize)],
    caps: &[f32],
    selector: &[bool],
) -> Vec<f32> {
    let mut pool = Vec::<f32>::new();
    for (idx, &(a, b)) in spans.iter().enumerate() {
        if selector[idx] {
            let cap = caps[idx];
            for &v in &y[a..b] {
                if is_finite_positive(v) && v <= cap {
                    pool.push(v);
                }
            }
        }
    }
    pool
}
