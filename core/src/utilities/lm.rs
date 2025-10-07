//! Levenberg–Marquardt (LM) — nonlinear least-squares curve fitting
//!
//! Rust implementation of JavaScript version from:
//! https://github.com/mljs/levenberg-marquardt
//!
//! # References
//! * Kenneth M. Brown & J. E. Dennis Jr. (1971).
//!   “Derivative-free analogues of the Levenberg–Marquardt and Gauss algorithms
//!   for nonlinear least-squares approximation.” DOI: https://doi.org/10.1007/BF01404679
//! * Henri P. Gavin (2024).
//!   “The Levenberg–Marquardt algorithm for nonlinear least squares curve-fitting problems.”
//! * MLJS `levenberg-marquardt` repository.

use std::time::{Duration, Instant};

#[derive(Clone, Debug)]
pub struct Data2D {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

#[derive(Clone, Debug)]
pub enum Weights {
    Scalar(f64),
    Array(Vec<f64>),
}

#[derive(Clone, Debug)]
pub enum GradientDifference {
    Scalar(f64),
    Array(Vec<f64>),
}

#[derive(Clone, Debug)]
pub struct LevenbergMarquardtOptions {
    pub timeout: Option<f64>,
    pub initial_values: Vec<f64>,
    pub weights: Option<Weights>,
    pub damping: Option<f64>,
    pub damping_step_up: Option<f64>,
    pub damping_step_down: Option<f64>,
    pub max_iterations: Option<usize>,
    pub error_tolerance: Option<f64>,
    pub central_difference: Option<bool>,
    pub gradient_difference: Option<GradientDifference>,
    pub improvement_threshold: Option<f64>,
    pub min_values: Option<Vec<f64>>,
    pub max_values: Option<Vec<f64>>,
}

impl Default for LevenbergMarquardtOptions {
    fn default() -> Self {
        Self {
            timeout: Some(0.01),
            initial_values: Vec::new(),
            weights: Some(Weights::Scalar(1.0)),
            damping: Some(1e-2),
            damping_step_up: Some(11.0),
            damping_step_down: Some(9.0),
            max_iterations: None,
            error_tolerance: Some(1e-10),
            central_difference: Some(false),
            gradient_difference: None,
            improvement_threshold: Some(1e-3),
            min_values: None,
            max_values: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct CheckedOptions {
    pub check_timeout_deadline: Option<Instant>,
    pub min_values: Vec<f64>,
    pub max_values: Vec<f64>,
    pub parameters: Vec<f64>,
    pub weight_square: Vec<f64>,
    pub damping: f64,
    pub damping_step_up: f64,
    pub damping_step_down: f64,
    pub max_iterations: usize,
    pub error_tolerance: f64,
    pub central_difference: bool,
    pub gradient_difference: Vec<f64>,
    pub improvement_threshold: f64,
}

#[derive(Clone, Debug)]
pub struct LevenbergMarquardtReturn {
    pub parameter_values: Vec<f64>,
    pub parameter_error: f64,
    pub iterations: usize,
}

pub type ParameterizedFunction = dyn Fn(&[f64]) -> Box<dyn Fn(f64) -> f64>;

pub fn lm(
    data: &Data2D,
    parameterized_function: &ParameterizedFunction,
    options: &LevenbergMarquardtOptions,
) -> Result<LevenbergMarquardtReturn, String> {
    let checked_options = check_options(data, options)?;
    let mut damping = checked_options.damping;

    let check_timeout = |deadline: &Option<Instant>| -> bool {
        if let Some(t) = deadline {
            Instant::now() > *t
        } else {
            false
        }
    };

    let min_values = checked_options.min_values.clone();
    let max_values = checked_options.max_values.clone();
    let mut parameters = checked_options.parameters.clone();
    let weight_square = checked_options.weight_square.clone();
    let damping_step_up = checked_options.damping_step_up;
    let damping_step_down = checked_options.damping_step_down;
    let max_iterations = checked_options.max_iterations;
    let error_tolerance = checked_options.error_tolerance;
    let central_difference = checked_options.central_difference;
    let gradient_difference = checked_options.gradient_difference.clone();
    let improvement_threshold = checked_options.improvement_threshold;

    let mut error = error_calculation(data, &parameters, parameterized_function, &weight_square);
    let mut optimal_error = error;
    let mut optimal_parameters = parameters.clone();

    let mut converged = error <= error_tolerance;

    let mut iteration: usize = 0;
    while iteration < max_iterations && !converged {
        let previous_error = error;

        let step_out = step(
            data,
            &parameters,
            damping,
            &gradient_difference,
            parameterized_function,
            central_difference,
            Some(&weight_square),
        );

        let perturbations = step_out.perturbations;
        let jacobian_weight_residual_error = step_out.jacobian_weight_residual_error;

        for k in 0..parameters.len() {
            let new_val = parameters[k] - perturbations[k];
            parameters[k] = new_val.max(min_values[k]).min(max_values[k]);
        }

        error = error_calculation(data, &parameters, parameterized_function, &weight_square);

        if error.is_nan() {
            break;
        }

        if error < optimal_error - error_tolerance {
            optimal_error = error;
            optimal_parameters = parameters.clone();
        }

        let denom: f64 = perturbations
            .iter()
            .zip(jacobian_weight_residual_error.iter())
            .map(|(p, jwr)| p * (damping * p + jwr))
            .sum();

        let improvement_metric = if denom.abs() > 1e-18 && denom.is_finite() {
            (previous_error - error) / denom
        } else {
            0.0
        };

        if improvement_metric > improvement_threshold {
            damping = (damping / damping_step_down).max(1e-7);
        } else {
            damping = (damping * damping_step_up).min(1e7);
        }

        if check_timeout(&checked_options.check_timeout_deadline) {
            let secs = options.timeout.unwrap_or(0.0);
            return Err(format!("The execution time is over to {} seconds", secs));
        }

        converged = error <= error_tolerance;
        iteration += 1;
    }

    Ok(LevenbergMarquardtReturn {
        parameter_values: optimal_parameters,
        parameter_error: optimal_error,
        iterations: iteration,
    })
}

pub struct StepReturn {
    pub perturbations: Vec<f64>,
    pub jacobian_weight_residual_error: Vec<f64>,
}

pub fn step(
    data: &Data2D,
    params: &[f64],
    damping: f64,
    gradient_difference: &[f64],
    parameterized_function: &ParameterizedFunction,
    central_difference: bool,
    weights: Option<&[f64]>,
) -> StepReturn {
    let n_params = params.len();
    let eye_l = eye(n_params, damping.max(1e-12));

    let func = parameterized_function(params);
    let mut evaluated_data = vec![0.0_f64; data.x.len()];
    for i in 0..data.x.len() {
        evaluated_data[i] = func(data.x[i]);
    }

    let j = gradient_function(
        data,
        &evaluated_data,
        params,
        gradient_difference,
        parameterized_function,
        central_difference,
    );

    let residual = matrix_function(data, &evaluated_data);

    let j_t = transpose(&j);
    let j_t_w = scale_rows(&j_t, weights);

    let a = add(&eye_l, &mmul(&j, &j_t_w));

    let wr = scale_vec_rows(&residual, weights);
    let b = mat_vec_mul(&j, &wr);

    let perturbations = match solve_sym_posdef_or_svd(&a, &b) {
        Some(p) => p,
        None => vec![0.0; n_params],
    };

    StepReturn {
        perturbations,
        jacobian_weight_residual_error: b,
    }
}

fn solve_sym_posdef_or_svd(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let mut a_bumped = a.to_vec();
    let n = a_bumped.len();
    for i in 0..n {
        a_bumped[i][i] += 1e-12;
    }
    if let Some(x) = solve_linear(&a_bumped, b) {
        return Some(x);
    }
    if let Ok(inv) = invert(&a_bumped) {
        return Some(mat_vec_mul(&inv, b));
    }
    None
}

fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    let mut aug = vec![vec![0.0_f64; n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j];
        }
        aug[i][n] = b[i];
    }
    for i in 0..n {
        let mut pivot = i;
        let mut maxv = aug[i][i].abs();
        for r in (i + 1)..n {
            let v = aug[r][i].abs();
            if v > maxv {
                maxv = v;
                pivot = r;
            }
        }
        if maxv <= 0.0 || !maxv.is_finite() {
            return None;
        }
        if pivot != i {
            aug.swap(i, pivot);
        }
        let diag = aug[i][i];
        for j in i..=n {
            aug[i][j] /= diag;
        }
        for r in 0..n {
            if r == i {
                continue;
            }
            let factor = aug[r][i];
            for j in i..=n {
                aug[r][j] -= factor * aug[i][j];
            }
        }
    }
    Some((0..n).map(|i| aug[i][n]).collect())
}

pub fn error_calculation(
    data: &Data2D,
    parameters: &[f64],
    parameterized_function: &ParameterizedFunction,
    weight_square: &[f64],
) -> f64 {
    let func = parameterized_function(parameters);
    let mut error = 0.0_f64;
    for i in 0..data.x.len() {
        let r = data.y[i] - func(data.x[i]);
        error += (r * r) / weight_square[i];
    }
    error
}

pub fn gradient_function(
    data: &Data2D,
    evaluated_data: &[f64],
    params: &[f64],
    gradient_difference: &[f64],
    param_function: &ParameterizedFunction,
    central_difference: bool,
) -> Vec<Vec<f64>> {
    let nb_params = params.len();
    let nb_points = data.x.len();
    let mut ans = vec![vec![0.0_f64; nb_points]; nb_params];

    let mut row_index = 0usize;
    for param in 0..nb_params {
        if gradient_difference[param] == 0.0 {
            continue;
        }
        let mut delta = gradient_difference[param];
        let mut aux_params = params.to_vec();
        aux_params[param] += delta;
        let func_param = param_function(&aux_params);
        if !central_difference {
            for point in 0..nb_points {
                ans[row_index][point] = (evaluated_data[point] - func_param(data.x[point])) / delta;
            }
        } else {
            aux_params = params.to_vec();
            aux_params[param] -= delta;
            delta *= 2.0;
            let func_param2 = param_function(&aux_params);
            for point in 0..nb_points {
                ans[row_index][point] =
                    (func_param2(data.x[point]) - func_param(data.x[point])) / delta;
            }
        }
        row_index += 1;
    }

    ans
}

pub fn check_options(
    data: &Data2D,
    options: &LevenbergMarquardtOptions,
) -> Result<CheckedOptions, String> {
    let damping = options.damping.unwrap_or(1e-2);
    let damping_step_up = options.damping_step_up.unwrap_or(11.0);
    let damping_step_down = options.damping_step_down.unwrap_or(9.0);
    let max_iterations = options.max_iterations.unwrap_or(100);
    let error_tolerance = options.error_tolerance.unwrap_or(1e-7);
    let central_difference = options.central_difference.unwrap_or(false);
    let improvement_threshold = options.improvement_threshold.unwrap_or(1e-3);

    if damping <= 0.0 {
        return Err("The damping option must be a positive number".into());
    } else if data.x.is_empty() || data.y.is_empty() {
        return Err("The data parameter must have x and y elements".into());
    } else if data.x.len() < 2 || data.y.len() < 2 {
        return Err("The data parameter elements must be an array with more than 2 points".into());
    } else if data.x.len() != data.y.len() {
        return Err("The data parameter elements must have the same size".into());
    }

    if options.initial_values.is_empty() {
        return Err("The initialValues option is mandatory and must be an array".into());
    }
    let parameters = options.initial_values.clone();
    let par_len = parameters.len();

    let mut max_values = options
        .max_values
        .clone()
        .unwrap_or_else(|| vec![f64::MAX; par_len]);
    let mut min_values = options
        .min_values
        .clone()
        .unwrap_or_else(|| vec![f64::MIN; par_len]);

    if max_values.len() != min_values.len() {
        return Err("minValues and maxValues must be the same size".into());
    }

    let gradient_difference_array = get_gradient_difference_array(
        options
            .gradient_difference
            .clone()
            .unwrap_or(GradientDifference::Scalar(1.0e-1)),
        &parameters,
    );

    let check_timeout_deadline = get_check_timeout(options.timeout);

    let filler = get_filler(options.weights.clone(), data.x.len())?;
    let mut weight_square = vec![0.0_f64; data.x.len()];
    for i in 0..data.x.len() {
        weight_square[i] = filler(i);
    }

    if max_values.len() != par_len {
        max_values = vec![f64::MAX; par_len];
    }
    if min_values.len() != par_len {
        min_values = vec![f64::MIN; par_len];
    }

    Ok(CheckedOptions {
        check_timeout_deadline,
        min_values,
        max_values,
        parameters,
        weight_square,
        damping,
        damping_step_up,
        damping_step_down,
        max_iterations,
        error_tolerance,
        central_difference,
        gradient_difference: gradient_difference_array,
        improvement_threshold,
    })
}

fn matrix_function(data: &Data2D, evaluated_data: &[f64]) -> Vec<f64> {
    let m = data.x.len();
    let mut ans = vec![0.0_f64; m];
    for i in 0..m {
        ans[i] = data.y[i] - evaluated_data[i];
    }
    ans
}

fn eye(n: usize, value: f64) -> Vec<Vec<f64>> {
    let mut m = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        m[i][i] = value;
    }
    m
}

fn add(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let r = a.len();
    let c = a[0].len();
    let mut out = vec![vec![0.0_f64; c]; r];
    for i in 0..r {
        for j in 0..c {
            out[i][j] = a[i][j] + b[i][j];
        }
    }
    out
}

fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let r = a.len();
    let c = a[0].len();
    let mut out = vec![vec![0.0_f64; r]; c];
    for i in 0..r {
        for j in 0..c {
            out[j][i] = a[i][j];
        }
    }
    out
}

fn scale_rows(m: &[Vec<f64>], weights: Option<&[f64]>) -> Vec<Vec<f64>> {
    let r = m.len();
    let c = m[0].len();
    let mut out = vec![vec![0.0_f64; c]; r];
    if let Some(w) = weights {
        for i in 0..r {
            let wv = w[i];
            for j in 0..c {
                out[i][j] = m[i][j] * wv;
            }
        }
    } else {
        for i in 0..r {
            for j in 0..c {
                out[i][j] = m[i][j];
            }
        }
    }
    out
}

fn scale_vec_rows(v: &[f64], weights: Option<&[f64]>) -> Vec<f64> {
    let mut out = vec![0.0_f64; v.len()];
    if let Some(w) = weights {
        for i in 0..v.len() {
            out[i] = v[i] * w[i];
        }
    } else {
        out.copy_from_slice(v);
    }
    out
}

fn mmul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let r = a.len();
    let m = a[0].len();
    let c = b[0].len();
    let mut out = vec![vec![0.0_f64; c]; r];
    for i in 0..r {
        for k in 0..m {
            let aik = a[i][k];
            for j in 0..c {
                out[i][j] += aik * b[k][j];
            }
        }
    }
    out
}

fn mat_vec_mul(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    let r = a.len();
    let c = a[0].len();
    let mut out = vec![0.0_f64; r];
    for i in 0..r {
        let mut s = 0.0_f64;
        for j in 0..c {
            s += a[i][j] * x[j];
        }
        out[i] = s;
    }
    out
}

fn invert(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    let mut aug = vec![vec![0.0_f64; 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j];
        }
        aug[i][n + i] = 1.0;
    }

    for i in 0..n {
        let mut pivot = i;
        let mut maxv = aug[i][i].abs();
        for r in (i + 1)..n {
            let v = aug[r][i].abs();
            if v > maxv {
                maxv = v;
                pivot = r;
            }
        }
        if maxv == 0.0 {
            return Err("singular matrix".into());
        }
        if pivot != i {
            aug.swap(i, pivot);
        }

        let diag = aug[i][i];
        for j in 0..(2 * n) {
            aug[i][j] /= diag;
        }

        for r in 0..n {
            if r == i {
                continue;
            }
            let factor = aug[r][i];
            for j in 0..(2 * n) {
                aug[r][j] -= factor * aug[i][j];
            }
        }
    }

    let mut inv = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            inv[i][j] = aug[i][n + j];
        }
    }
    Ok(inv)
}

fn get_gradient_difference_array(
    gradient_difference: GradientDifference,
    parameters: &[f64],
) -> Vec<f64> {
    match gradient_difference {
        GradientDifference::Scalar(v) => vec![v; parameters.len()],
        GradientDifference::Array(arr) => {
            if arr.len() != parameters.len() {
                vec![arr[0]; parameters.len()]
            } else {
                arr
            }
        }
    }
}

fn get_filler(
    weights: Option<Weights>,
    data_length: usize,
) -> Result<Box<dyn Fn(usize) -> f64>, String> {
    match weights {
        Some(Weights::Scalar(w)) => {
            let value = 1.0 / (w * w);
            Ok(Box::new(move |_| value))
        }
        Some(Weights::Array(arr)) => {
            if arr.len() < data_length {
                let value = 1.0 / (arr[0] * arr[0]);
                Ok(Box::new(move |_| value))
            } else {
                Ok(Box::new(move |i: usize| 1.0 / (arr[i] * arr[i])))
            }
        }
        None => {
            let value = 1.0 / (1.0 * 1.0);
            Ok(Box::new(move |_| value))
        }
    }
}

fn get_check_timeout(timeout: Option<f64>) -> Option<Instant> {
    if let Some(secs) = timeout {
        Some(Instant::now() + Duration::from_secs_f64(secs.max(0.0)))
    } else {
        None
    }
}
