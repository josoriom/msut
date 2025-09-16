use crate::utilities::{
    closest_index,
    find_noise_level::find_noise_level,
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

const DEFAULT_WINDOW_SIZES: &[usize] = &[5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33];

pub fn get_boundaries(data: &DataXY, peak_x: f64, opt: Option<BoundariesOptions>) -> Boundaries {
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

    let o = opt.unwrap_or_default();
    let idx = closest_index(&data.x, peak_x);

    let mut win_sizes = Vec::new();
    for &win in DEFAULT_WINDOW_SIZES {
        if win >= 5 && win % 2 == 1 && win <= n {
            win_sizes.push(win);
        }
    }

    let mut from_votes = vec![0u32; n];
    let mut to_votes = vec![0u32; n];
    let mut has_any_votes = false;
    let streak_limit = 2usize;
    let noise_level = find_noise_level(&data.y);

    for win in win_sizes {
        let smooth_y = sgg(
            &data.y,
            &data.x,
            SggOptions {
                window_size: win,
                derivative: 0,
                polynomial: 3,
            },
        );

        let peak_level = smooth_y[idx].abs().max(1.0);
        let mut tolerance = 0.01 * peak_level;
        if tolerance < 3.0 {
            tolerance = 3.0;
        }
        if tolerance > 250.0 {
            tolerance = 250.0;
        }

        let noise_floor = noise_level.max(0.035 * peak_level);

        let left_idx = walk_left(&smooth_y, idx, tolerance, streak_limit, noise_floor);
        let right_idx = walk_right(&smooth_y, idx, tolerance, streak_limit, noise_floor);

        if left_idx < right_idx {
            from_votes[left_idx] = from_votes[left_idx].saturating_add(1);
            to_votes[right_idx] = to_votes[right_idx].saturating_add(1);
            has_any_votes = true;
        }
    }

    if has_any_votes {
        let from_idx = pick_mode(&from_votes, false);
        let to_idx = pick_mode(&to_votes, true);
        if let (Some(f), Some(t)) = (from_idx, to_idx) {
            if f < t {
                return Boundaries {
                    from: Boundary {
                        index: Some(f),
                        value: Some(data.x[f]),
                    },
                    to: Boundary {
                        index: Some(t),
                        value: Some(data.x[t]),
                    },
                };
            }
        }
    }

    let ws = o.window_size;
    let y_buf;
    let y_slice: &[f32] = if ws >= 5 && ws % 2 == 1 && ws <= n {
        y_buf = sgg(
            &data.y,
            &data.x,
            SggOptions {
                window_size: ws,
                derivative: 0,
                polynomial: 3,
            },
        );
        &y_buf
    } else {
        &data.y
    };

    Boundaries {
        from: step_boundary(&data.x, y_slice, idx, -1, o.epsilon),
        to: step_boundary(&data.x, y_slice, idx, 1, o.epsilon),
    }
}

fn walk_left(
    y: &[f32],
    start: usize,
    tolerance: f32,
    streak_limit: usize,
    noise_floor: f32,
) -> usize {
    if y.is_empty() {
        return 0;
    }
    let mut i = if start > 0 { start - 1 } else { 0 };
    let mut streak = 0usize;
    while i > 0 {
        let current = y[i];
        if current <= noise_floor {
            break;
        }
        let previous = y[i - 1];
        if previous + tolerance <= current {
            streak = 0;
        } else {
            streak += 1;
            if streak >= streak_limit {
                break;
            }
        }
        i -= 1;
    }
    i
}

fn walk_right(
    y: &[f32],
    start: usize,
    tolerance: f32,
    streak_limit: usize,
    noise_floor: f32,
) -> usize {
    let n = y.len();
    if n == 0 {
        return 0;
    }
    let mut i = if start + 1 < n { start + 1 } else { n - 1 };
    let mut streak = 0usize;
    let mut min_value = y[i];
    while i + 1 < n {
        let current = y[i];
        if current <= noise_floor {
            break;
        }
        let next_value = y[i + 1];
        if next_value + tolerance <= current {
            streak = 0;
        } else {
            streak += 1;
            if streak >= streak_limit {
                if next_value > min_value + tolerance {
                    break;
                }
                streak = 0;
            }
        }
        if next_value < min_value {
            min_value = next_value;
        }
        i += 1;
    }
    i
}

fn pick_mode(votes: &[u32], right_bias: bool) -> Option<usize> {
    let mut best = None::<usize>;
    let mut best_val = 0u32;
    for i in 0..votes.len() {
        let v = votes[i];
        if v == 0 {
            continue;
        }
        match best {
            None => {
                best = Some(i);
                best_val = v;
            }
            Some(b) => {
                let tie_right = right_bias && i > b;
                let tie_left = !right_bias && i < b;
                if v > best_val || (v == best_val && (tie_right || tie_left)) {
                    best = Some(i);
                    best_val = v;
                }
            }
        }
    }
    best
}

fn step_boundary(x: &[f64], y: &[f32], idx: usize, inc: isize, epsilon: f64) -> Boundary {
    let n = x.len() as isize;
    if n < 2 {
        return Boundary {
            index: None,
            value: None,
        };
    }
    let mut p = idx as isize + inc;
    while p >= 0 && p < n {
        let q = p + inc;
        if q < 0 || q >= n {
            break;
        }
        let i = p as usize;
        let j = q as usize;
        let dx = x[j] - x[i];
        let dxn = if dx == 0.0 {
            if inc > 0 { epsilon } else { -epsilon }
        } else {
            dx
        };
        let dy = (y[j] as f64) - (y[i] as f64);
        let s = dy * dxn;
        if (inc > 0 && s >= 0.0) || (inc < 0 && s <= 0.0) {
            return Boundary {
                index: Some(i),
                value: Some(x[i]),
            };
        }
        p += inc;
    }
    Boundary {
        index: None,
        value: None,
    }
}
