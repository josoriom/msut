#[derive(Clone, Copy, Debug)]
pub struct SggOptions {
    pub window_size: usize,
    pub derivative: usize,
    pub polynomial: usize,
}

impl Default for SggOptions {
    fn default() -> Self {
        Self {
            window_size: 9,
            derivative: 0,
            polynomial: 3,
        }
    }
}

pub fn sgg(ys: &[f32], xs: &[f64], opts: SggOptions) -> Vec<f32> {
    let window_size = opts.window_size;
    let derivative = opts.derivative;
    let polynomial = opts.polynomial;

    if window_size % 2 == 0 || window_size < 5 {
        panic!("Invalid window size (should be odd and at least 5 integer number)");
    }
    if ys.is_empty() {
        panic!("Y values must be an array");
    }
    if xs.is_empty() {
        panic!("X must be defined");
    }
    if window_size > ys.len() {
        panic!(
            "Window size is higher than the data length {}>{}",
            window_size,
            ys.len()
        );
    }
    if polynomial < 1 {
        panic!("Polynomial should be a positive integer");
    }
    if polynomial >= 6 {
        eprintln!(
            "You should not use polynomial grade higher than 5 if you are not sure that your data arises from such a model. Possible polynomial oscillation problems"
        );
    }

    let half = window_size / 2;
    let n = ys.len();

    let mut y64 = vec![0.0f64; n];
    for i in 0..n {
        y64[i] = ys[i] as f64;
    }

    let hs = get_hs(xs, half, derivative);
    let weights = full_weights(window_size, polynomial, derivative);

    let mut ans = vec![0.0f32; n];

    for i in 0..half {
        let wl = &weights[half - i - 1];
        let wr = &weights[half + i + 1];

        let mut dl = 0.0f64;
        let mut dr = 0.0f64;

        for l in 0..window_size {
            dl += wl[l] * y64[l];
            dr += wr[l] * y64[n - window_size + l];
        }

        let idx_l = half - i - 1;
        let idx_r = n - half + i;

        ans[idx_l] = (dl / hs[idx_l]) as f32;
        ans[idx_r] = (dr / hs[idx_r]) as f32;
    }

    let wc = &weights[half];
    for i in window_size..=n {
        let mut d = 0.0f64;
        for l in 0..window_size {
            d += wc[l] * y64[l + i - window_size];
        }
        let idx = i - half - 1;
        ans[idx] = (d / hs[idx]) as f32;
    }

    ans
}

fn get_hs(xs: &[f64], half: usize, derivative: usize) -> Vec<f64> {
    let n = xs.len();
    if derivative == 0 || n < 2 {
        let v = vec![1.0f64; n];
        return v;
    }

    let mut pref = vec![0.0f64; n];
    for i in 0..(n - 1) {
        pref[i + 1] = pref[i] + (xs[i + 1] - xs[i]);
    }

    let mut hs = vec![1.0f64; n];
    for c in 0..n {
        let start = if c >= half { c - half } else { 0 };
        let end_excl = {
            let e = c + half;
            if e < n - 1 { e } else { n - 1 }
        };
        let count = if end_excl > start {
            end_excl - start
        } else {
            0
        };
        let avg = if count > 0 {
            (pref[end_excl] - pref[start]) / (count as f64)
        } else {
            1.0
        };
        hs[c] = avg.powi(derivative as i32);
    }
    hs
}

fn full_weights(m: usize, n: usize, s: usize) -> Vec<Vec<f64>> {
    let half = (m / 2) as i32;
    let n_i = n as i32;
    let s_i = s as i32;

    let mut gi: Vec<Vec<f64>> = Vec::with_capacity(m);
    for idx in 0..m {
        let i_off = idx as i32 - half;
        let tbl = gram_table(i_off, half, n_i, 0);
        let mut row = vec![0.0f64; (n_i as usize) + 1];
        for k in 0..=n_i as usize {
            row[k] = tbl[k][0];
        }
        gi.push(row);
    }

    let mut gt: Vec<Vec<f64>> = Vec::with_capacity(m);
    for idx in 0..m {
        let t_off = idx as i32 - half;
        let tbl = gram_table(t_off, half, n_i, s_i);
        let mut row = vec![0.0f64; (n_i as usize) + 1];
        for k in 0..=n_i as usize {
            row[k] = tbl[k][s_i as usize];
        }
        gt.push(row);
    }

    let two_m = 2 * half;
    let mut coef = vec![0.0f64; (n_i as usize) + 1];
    for k in 0..=n_i {
        let num = gen_fact(two_m, k);
        let den = gen_fact(two_m + k + 1, k + 1);
        coef[k as usize] = (2 * k + 1) as f64 * (num / den);
    }

    let mut w = vec![vec![0.0f64; m]; m];
    for t_idx in 0..m {
        for j_idx in 0..m {
            let mut sum = 0.0f64;
            for k in 0..=n_i as usize {
                sum += coef[k] * gi[j_idx][k] * gt[t_idx][k];
            }
            w[t_idx][j_idx] = sum;
        }
    }
    w
}

fn gram_table(i: i32, m: i32, n_max: i32, s_max: i32) -> Vec<Vec<f64>> {
    let nm = (n_max as usize) + 1;
    let sm = (s_max as usize) + 1;
    let mut g = vec![vec![0.0f64; sm]; nm];
    g[0][0] = 1.0;

    let mut k = 1;
    while k <= n_max {
        let kf = k as f64;
        let denom = kf * (2 * m - k + 1) as f64;
        let a = (4 * k - 2) as f64 / denom;
        let b = ((k - 1) as f64 * (2 * m + k) as f64) / denom;

        let mut s = 0;
        while s <= s_max {
            let s_usize = s as usize;
            let term1 = (i as f64) * g[(k - 1) as usize][s_usize];
            let term2 = if s > 0 {
                (s as f64) * g[(k - 1) as usize][(s - 1) as usize]
            } else {
                0.0
            };
            let term3 = if k >= 2 {
                g[(k - 2) as usize][s_usize]
            } else {
                0.0
            };
            g[k as usize][s_usize] = a * (term1 + term2) - b * term3;
            s += 1;
        }
        k += 1;
    }
    g
}

fn gen_fact(a: i32, b: i32) -> f64 {
    if a >= b {
        let start = a - b + 1;
        let mut acc = 1.0f64;
        let mut j = start;
        while j <= a {
            acc *= j as f64;
            j += 1;
        }
        acc
    } else {
        1.0
    }
}
