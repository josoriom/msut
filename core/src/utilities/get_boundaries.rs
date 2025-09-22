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

    let mut tail_started = false;
    let mut checking = false;
    let mut steps_up: usize = 0;

    let mut local_min_idx: isize = start as isize;
    let mut local_min_val: f64 = y[start] as f64;

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
        let is_desc = slope_dir < 0.0;
        let is_asc_or_flat = !is_desc;

        if is_desc {
            tail_started = true;
        }

        if tail_started {
            if (y[j] as f64) < local_min_val {
                local_min_val = y[j] as f64;
                local_min_idx = next;
            }

            if !checking {
                if is_asc_or_flat {
                    checking = true;
                    steps_up = 1;
                }
            } else {
                if is_desc {
                    checking = false;
                    steps_up = 0;
                } else {
                    steps_up += 1;
                    if steps_up >= n_steps {
                        let rise = (y[j] as f64) - local_min_val;
                        if rise >= noise {
                            let k = local_min_idx.clamp(0, n - 1) as usize;
                            return Boundary {
                                index: Some(k),
                                value: Some(x[k]),
                            };
                        } else {
                            checking = false;
                            steps_up = 0;
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
                let k = first.clamp(0, n - 1) as usize;
                return Boundary {
                    index: Some(k),
                    value: Some(x[k]),
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
