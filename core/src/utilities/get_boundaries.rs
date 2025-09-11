use crate::utilities::{
    sgg::{SggOptions, sgg},
    structs::DataXY,
};

#[derive(Clone, Copy, Debug)]
pub struct Boundary {
    pub index: Option<usize>,
    pub value: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct Boundaries {
    pub from: Boundary,
    pub to: Boundary,
}

#[derive(Clone, Copy, Debug)]
pub struct BoundariesOptions {
    pub epsilon: f64,
    pub window_size: usize,
}
impl Default for BoundariesOptions {
    fn default() -> Self {
        Self {
            epsilon: 1e-5,
            window_size: 17,
        }
    }
}

fn sg_is_valid(sg: SggOptions) -> bool {
    sg.window_size >= 5
        && sg.window_size % 2 == 1
        && sg.polynomial >= 1
        && sg.polynomial < sg.window_size
}

pub fn get_boundaries(
    data: &DataXY,
    peak_position: f64,
    options: Option<BoundariesOptions>,
) -> Boundaries {
    let n = data.x.len();
    if n < 2 || n != data.y.len() {
        return Boundaries {
            from: Boundary {
                index: None,
                value: None,
            },
            to: Boundary {
                index: None,
                value: None,
            },
        };
    }

    let opt = options.unwrap_or_default();
    let sg = SggOptions {
        window_size: opt.window_size,
        derivative: 0,
        polynomial: 3,
    };

    let idx = find_closest_index(&data.x, peak_position);

    let maybe_smooth = if sg_is_valid(sg) && sg.window_size <= n {
        Some(sgg(&data.y, &data.x, sg))
    } else {
        None
    };
    let y_ref: &[f32] = match &maybe_smooth {
        Some(v) => v.as_slice(),
        None => &data.y,
    };

    Boundaries {
        from: get_boundary(&data.x, y_ref, idx, -1, opt.epsilon),
        to: get_boundary(&data.x, y_ref, idx, 1, opt.epsilon),
    }
}

fn get_boundary(x: &[f64], y: &[f32], index: usize, inc: isize, epsilon: f64) -> Boundary {
    let n = x.len();
    if n < 2 {
        return Boundary {
            index: None,
            value: None,
        };
    }

    let n_isize = n as isize;
    let mut pos = index as isize + inc;

    while pos >= 0 && pos < n_isize {
        let j = pos + inc;
        if j < 0 || j >= n_isize {
            break;
        }

        let i = pos as usize;
        let jj = j as usize;

        let dx = x[jj] - x[i];
        let dx_eff = if dx == 0.0 {
            if inc > 0 { epsilon } else { -epsilon }
        } else {
            dx
        };

        let dy = (y[jj] as f64) - (y[i] as f64);
        let s = dy * dx_eff;

        if (inc > 0 && s >= 0.0) || (inc < 0 && s <= 0.0) {
            return Boundary {
                index: Some(i),
                value: Some(x[i]),
            };
        }
        pos += inc;
    }

    Boundary {
        index: None,
        value: None,
    }
}

fn find_closest_index(xs: &[f64], v: f64) -> usize {
    if xs.is_empty() {
        return 0;
    }
    let mut lo = 0usize;
    let mut hi = xs.len();
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
