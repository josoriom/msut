use crate::utilities::{closest_index, structs::DataXY};

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
    pub noise: f64,
    pub n_steps: usize,
    pub baseline_run: usize,
}
impl Default for BoundariesOptions {
    fn default() -> Self {
        Self {
            epsilon: 1e-5,
            noise: 0.0,
            n_steps: 3,
            baseline_run: 5,
        }
    }
}

pub fn get_boundaries(
    data: &DataXY,
    peak_x: f64,
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

    let opts = options.unwrap_or_default();
    let idx = closest_index(&data.x, peak_x);

    let global_min = data.y.iter().fold(f64::INFINITY, |m, &v| m.min(v as f64));

    let from = walk(
        &data.x,
        &data.y,
        idx,
        -1,
        opts.epsilon,
        opts.noise,
        opts.n_steps,
        opts.baseline_run,
        global_min,
    );

    let to = walk(
        &data.x,
        &data.y,
        idx,
        1,
        opts.epsilon,
        opts.noise,
        opts.n_steps,
        opts.baseline_run,
        global_min,
    );

    Boundaries { from, to }
}

fn walk(
    x: &[f64],
    y: &[f32],
    start: usize,
    direction: isize,
    epsilon: f64,
    noise: f64,
    n_steps: usize,
    baseline_run: usize,
    global_min: f64,
) -> Boundary {
    let n = x.len() as isize;
    if n < 2 {
        return Boundary {
            index: None,
            value: None,
        };
    }

    let mut current = start as isize + if direction > 0 { 1 } else { -1 };
    let min_stop = global_min + epsilon;

    let mut baseline_count: usize = 0;
    let mut is_tail = false;
    let mut is_checking = false;
    let mut check_steps: usize = 0;
    let mut check_start_idx: Option<usize> = None;
    let mut check_start_val: f64 = 0.0;

    while current >= 0 && current < n - 1 {
        let next = current + if direction > 0 { 1 } else { -1 };
        if next < 0 || next >= n {
            break;
        }

        let i = current as usize;
        let j = next as usize;

        if (y[i] as f64) <= min_stop {
            return Boundary {
                index: Some(i),
                value: Some(x[i]),
            };
        }

        let dx = x[j] - x[i];
        let denom = if dx.abs() < epsilon {
            if direction > 0 { epsilon } else { -epsilon }
        } else {
            dx
        };
        let dy = (y[j] as f64) - (y[i] as f64);
        let slope_dir = (dy / denom) * if direction > 0 { 1.0 } else { -1.0 };
        let is_neg = slope_dir < 0.0;

        if !is_tail && is_neg {
            is_tail = true;
        }

        if is_tail {
            if !is_checking {
                if !is_neg {
                    is_checking = true;
                    check_steps = 1;
                    check_start_idx = Some(i);
                    check_start_val = y[i] as f64;
                }
            } else {
                if is_neg {
                    is_checking = false;
                    check_steps = 0;
                    check_start_idx = None;
                } else {
                    check_steps += 1;
                    if check_steps >= n_steps {
                        if let Some(cs) = check_start_idx {
                            let rise = ((y[j] as f64) - check_start_val).max(0.0);
                            if rise >= noise {
                                return Boundary {
                                    index: Some(cs),
                                    value: Some(x[cs]),
                                };
                            } else {
                                is_checking = false;
                                check_steps = 0;
                                check_start_idx = None;
                            }
                        }
                    }
                }
            }
        }

        if (y[i] as f64) <= noise {
            baseline_count += 1;
            if baseline_count >= baseline_run.max(1) {
                let first = if direction > 0 {
                    current - (baseline_run as isize - 1)
                } else {
                    current + (baseline_run as isize - 1)
                };
                let clamped = first.clamp(0, n - 1) as usize;
                return Boundary {
                    index: Some(clamped),
                    value: Some(x[clamped]),
                };
            }
        } else {
            baseline_count = 0;
        }

        current += direction;
    }

    Boundary {
        index: None,
        value: None,
    }
}
