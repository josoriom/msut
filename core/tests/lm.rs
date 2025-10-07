use msut::utilities::lm::{
    Data2D, GradientDifference, LevenbergMarquardtOptions, Weights, error_calculation, lm,
};

fn assert_close_decimals(a: f64, b: f64, decimals: usize) {
    let tol = 10f64.powi(-(decimals as i32));
    let diff = (a - b).abs();
    assert!(
        diff <= tol,
        "got {a}, expected {b}, diff {diff} > tol {tol}"
    );
}

fn assert_vec_close_decimals(a: &[f64], b: &[f64], decimals: usize) {
    assert_eq!(
        a.len(),
        b.len(),
        "length mismatch ({} vs {})",
        a.len(),
        b.len()
    );
    for (i, (ai, bi)) in a.iter().zip(b.iter()).enumerate() {
        let tol = 10f64.powi(-(decimals as i32));
        let diff = (ai - bi).abs();
        assert!(
            diff <= tol,
            "index {i}: got {ai}, expected {bi}, diff {diff} > tol {tol}"
        );
    }
}

fn linspace(n: usize, start: f64, end: f64) -> Vec<f64> {
    if n == 1 {
        return vec![start];
    }
    let step = (end - start) / ((n - 1) as f64);
    (0..n).map(|i| start + (i as f64) * step).collect()
}

fn sin_function(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
    let a = params[0];
    let b = params[1];
    Box::new(move |t: f64| a * (b * t).sin())
}

#[test]
fn linear_regression() {
    fn line(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
        let a = params[0];
        let b = params[1];
        Box::new(move |x: f64| a * x + b)
    }

    let x = vec![0., 1., 2., 3., 4., 5., 6.];
    let y = vec![-2., 0., 2., 4., 6., 8., 10.];
    let data = Data2D { x, y };

    let opts = LevenbergMarquardtOptions {
        initial_values: vec![1., 0.],
        timeout: None,
        weights: None,
        damping: None,
        damping_step_up: None,
        damping_step_down: None,
        max_iterations: None,
        error_tolerance: None,
        central_difference: None,
        gradient_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };

    let res = lm(&data, &line, &opts).unwrap();
    assert_vec_close_decimals(&res.parameter_values, &[2., -2.], 5);
    assert_close_decimals(res.parameter_error, 0.0, 9);
}

#[test]
fn contrived_bennet5() {
    fn f(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
        let b1 = params[0];
        let b2 = params[1];
        let b3 = params[2];
        Box::new(move |t: f64| b1 * (t + b2).powf(-1.0 / b3))
    }

    let n = 154usize;
    let x = linspace(n, -2.6581, 49.6526);
    let problem_parameters = vec![2., 3., 5.];
    let y: Vec<f64> = x.iter().map(|&t| f(&problem_parameters)(t)).collect();
    let data = Data2D { x, y };

    let opts = LevenbergMarquardtOptions {
        damping: Some(0.00001),
        max_iterations: Some(1000),
        error_tolerance: Some(1e-7),
        min_values: Some(vec![1., 2.7, 1.]),
        max_values: Some(vec![11., 11., 11.]),
        initial_values: vec![3.5, 3.8, 4.],
        timeout: None,
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        central_difference: None,
        gradient_difference: None,
        improvement_threshold: None,
    };

    let res = lm(&data, &f, &opts).unwrap();
    assert_vec_close_decimals(&res.parameter_values, &problem_parameters, 1);
    assert_close_decimals(res.parameter_error, 0.0, 2);
}

#[test]
fn contrived_sin() {
    let n = 20usize;
    let x = linspace(n, 0.0, 19.0);
    let problem_parameters = vec![2., 2.];
    let y: Vec<f64> = x
        .iter()
        .map(|&t| sin_function(&problem_parameters)(t))
        .collect();
    let data = Data2D { x, y };

    let opts = LevenbergMarquardtOptions {
        max_iterations: Some(100),
        gradient_difference: Some(GradientDifference::Scalar(1.0e-1)),
        damping: Some(0.1),
        damping_step_down: Some(1.0),
        damping_step_up: Some(1.0),
        initial_values: vec![3., 3.],
        timeout: None,
        weights: None,
        error_tolerance: None,
        central_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };

    let res = lm(&data, &sin_function, &opts).unwrap();
    assert_vec_close_decimals(&res.parameter_values, &problem_parameters, 1);
    assert_close_decimals(res.parameter_error, 0.0, 2);
}

#[test]
fn contrived_sigmoid() {
    fn sigmoid(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
        let a = params[0];
        let b = params[1];
        let c = params[2];
        Box::new(move |t: f64| a / (b + (-t * c).exp()))
    }

    let n = 20usize;
    let x = linspace(n, 0.0, 19.0);
    let problem_parameters = vec![2., 2., 2.];
    let y: Vec<f64> = x.iter().map(|&t| sigmoid(&problem_parameters)(t)).collect();
    let data = Data2D { x, y };

    let opts = LevenbergMarquardtOptions {
        damping: Some(0.1),
        initial_values: vec![3., 3., 3.],
        max_iterations: Some(200),
        timeout: None,
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        error_tolerance: None,
        central_difference: None,
        gradient_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };

    let res = lm(&data, &sigmoid, &opts).unwrap();
    assert_vec_close_decimals(&res.parameter_values, &problem_parameters, 1);
    assert_close_decimals(res.parameter_error, 0.0, 2);
}

fn lorentzians_model(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
    let p = params.to_vec();
    Box::new(move |t: f64| {
        let mut s = 0.0;
        let n = p.len();
        let mut i = 0usize;
        while i < n {
            let p2 = (p[i + 2] / 2.0).powi(2);
            let factor = p[i + 1] * p2;
            s += factor / (((t - p[i]).powi(2)) + p2);
            i += 3;
        }
        s
    })
}

#[test]
fn contrived_sum_of_lorentzians() {
    let n = 100usize;
    let x = linspace(n, 0.0, 99.0);
    let problem_parameters = vec![1.05, 0.1, 0.3, 4.0, 0.15, 0.3];
    let y: Vec<f64> = x
        .iter()
        .map(|&t| lorentzians_model(&problem_parameters)(t))
        .collect();
    let data = Data2D { x, y };

    let opts = LevenbergMarquardtOptions {
        damping: Some(0.01),
        gradient_difference: Some(GradientDifference::Array(vec![
            0.01, 0.0001, 0.0001, 0.01, 0.0001, 0.0,
        ])),
        initial_values: vec![1.1, 0.15, 0.29, 4.05, 0.17, 0.3],
        max_iterations: Some(500),
        timeout: None,
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        error_tolerance: None,
        central_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };

    let res = lm(&data, &lorentzians_model, &opts).unwrap();
    assert_vec_close_decimals(&res.parameter_values, &problem_parameters, 1);
    assert_close_decimals(res.parameter_error, 0.0, 2);
}

#[test]
fn contrived_sum_of_lorentzians_central() {
    let n = 100usize;
    let x = linspace(n, 0.0, 99.0);
    let problem_parameters = vec![1.0, 0.1, 0.3, 4.0, 0.15, 0.3];
    let y: Vec<f64> = x
        .iter()
        .map(|&t| lorentzians_model(&problem_parameters)(t))
        .collect();
    let data = Data2D { x, y };

    let opts = LevenbergMarquardtOptions {
        damping: Some(0.01),
        gradient_difference: Some(GradientDifference::Array(vec![
            0.01, 0.0001, 0.0001, 0.01, 0.0001,
        ])),
        central_difference: Some(true),
        initial_values: vec![1.1, 0.15, 0.29, 4.05, 0.17, 0.28],
        max_iterations: Some(500),
        error_tolerance: Some(1.0e-7),
        timeout: None,
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };

    let res = lm(&data, &lorentzians_model, &opts).unwrap();
    assert_vec_close_decimals(&res.parameter_values, &problem_parameters, 1);
    assert_close_decimals(res.parameter_error, 0.0, 2);
}

#[test]
fn should_return_solution_with_lowest_error() {
    let data = Data2D {
        x: vec![
            0.0,
            0.6283185307179586,
            1.2566370614359172,
            1.8849555921538759,
            2.5132741228718345,
            3.141592653589793,
            3.7699111843077517,
            4.39822971502571,
            5.026548245743669,
            5.654866776461628,
        ],
        y: vec![
            0.0,
            1.902113032590307,
            1.1755705045849465,
            -1.175570504584946,
            -1.9021130325903073,
            -4.898587196589413e-16,
            1.902113032590307,
            1.1755705045849467,
            -1.1755705045849456,
            -1.9021130325903075,
        ],
    };

    let opts = LevenbergMarquardtOptions {
        damping: Some(1.5),
        initial_values: vec![0.594398586701882, 0.3506424963635226],
        gradient_difference: Some(GradientDifference::Scalar(1e-2)),
        max_iterations: Some(100),
        error_tolerance: Some(1e-2),
        timeout: None,
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        central_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };

    let res = lm(&data, &sin_function, &opts).unwrap();

    let manual_error: f64 = data
        .x
        .iter()
        .zip(data.y.iter())
        .map(|(&x, &y)| {
            let y_hat = sin_function(&res.parameter_values)(x);
            (y - y_hat).powi(2)
        })
        .sum();

    assert_close_decimals(res.parameter_error, manual_error, 2);
    assert_close_decimals(res.parameter_error, 15.5, 1);
}

#[test]
fn parameter_error_linear_exact_fit() {
    fn linear_function(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
        let slope = params[0];
        let intercept = params[1];
        Box::new(move |x: f64| slope * x + intercept)
    }

    let sample_parameters = vec![1., 1.];
    let n = 10usize;
    let w = vec![1.0; n];
    let xs: Vec<f64> = (0..n).map(|i| i as f64).collect();
    let ys: Vec<f64> = xs
        .iter()
        .map(|&x| linear_function(&sample_parameters)(x))
        .collect();
    let data = Data2D { x: xs, y: ys };

    let err = error_calculation(&data, &sample_parameters, &linear_function, &w);
    assert_close_decimals(err, 0.0, 3);
}

#[test]
fn parameter_error_linear_offset_by_one() {
    fn linear_function(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
        let slope = params[0];
        let intercept = params[1];
        Box::new(move |x: f64| slope * x + intercept)
    }

    let sample_parameters = vec![1., 1.];
    let n = 10usize;
    let w = vec![1.0; n];
    let xs: Vec<f64> = (0..n).map(|i| i as f64).collect();
    let fct = linear_function(&sample_parameters);
    let ys: Vec<f64> = (0..n).map(|i| fct(i as f64)).collect();
    let data = Data2D { x: xs, y: ys };

    let mut parameters = sample_parameters.clone();
    parameters[1] += 1.0;

    let err = error_calculation(&data, &parameters, &linear_function, &w);
    assert_close_decimals(err, n as f64, 3);
}

#[test]
fn invalid_options_negative_damping() {
    let data = Data2D {
        x: vec![],
        y: vec![],
    };
    let opts = LevenbergMarquardtOptions {
        damping: Some(-1.0),
        initial_values: vec![],
        timeout: None,
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        max_iterations: None,
        error_tolerance: None,
        central_difference: None,
        gradient_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };
    let err = lm(&data, &|_| Box::new(|_| 1.0), &opts).err().unwrap();
    assert!(
        err.contains("The damping option must be a positive number"),
        "{}",
        err
    );
}

#[test]
fn invalid_options_initial_values_missing() {
    let input_data = Data2D {
        x: vec![1., 2.],
        y: vec![1., 2.],
    };
    let opts = LevenbergMarquardtOptions {
        damping: Some(0.1),
        initial_values: vec![],
        timeout: None,
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        max_iterations: None,
        error_tolerance: None,
        central_difference: None,
        gradient_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };
    let err = lm(&input_data, &sin_function, &opts).err().unwrap();
    assert!(
        err.contains("The initialValues option is mandatory and must be an array"),
        "{}",
        err
    );
}

#[test]
fn invalid_options_min_max_mismatch() {
    let data = Data2D {
        x: vec![1., 2.],
        y: vec![1., 2.],
    };
    let opts = LevenbergMarquardtOptions {
        damping: Some(0.1),
        min_values: Some(vec![1., 2., 3.]),
        max_values: Some(vec![1., 2.]),
        initial_values: vec![1., 1.],
        timeout: None,
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        max_iterations: None,
        error_tolerance: None,
        central_difference: None,
        gradient_difference: None,
        improvement_threshold: None,
    };
    let err = lm(&data, &sin_function, &opts).err().unwrap();
    assert!(
        err.contains("minValues and maxValues must be the same size"),
        "{}",
        err
    );
}

#[test]
fn invalid_data_too_short_or_mismatch() {
    let opts = LevenbergMarquardtOptions {
        damping: Some(0.1),
        initial_values: vec![3., 3.],
        timeout: None,
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        max_iterations: None,
        error_tolerance: None,
        central_difference: None,
        gradient_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };

    let data_short = Data2D {
        x: vec![1.0],
        y: vec![1.0],
    };
    let err1 = lm(&data_short, &sin_function, &opts).err().unwrap();
    assert!(
        err1.contains("must be an array with more than 2 points"),
        "{}",
        err1
    );

    let data_mismatch = Data2D {
        x: vec![1.0, 2.0],
        y: vec![1.0, 2.0, 3.0],
    };
    let err2 = lm(&data_mismatch, &sin_function, &opts).err().unwrap();
    assert!(err2.contains("must have the same size"), "{}", err2);
}

fn four_param_eq(params: &[f64]) -> Box<dyn Fn(f64) -> f64> {
    let a = params[0];
    let b = params[1];
    let c = params[2];
    let d = params[3];
    Box::new(move |t: f64| a + (b - a) / (1.0 + c.powf(d) * t.powf(-d)))
}

#[test]
fn ill_behaved_returns_initial_on_nan_after_start() {
    let data = Data2D {
        x: vec![
            9.22e-12, 5.53e-11, 3.32e-10, 1.99e-9, 1.19e-8, 7.17e-8, 4.3e-7, 0.00000258, 0.0000155,
            0.0000929,
        ],
        y: vec![
            7.807, -3.74, 21.119, 2.382, 4.269, 41.57, 73.401, 98.535, 97.059, 92.147,
        ],
    };

    let opts = LevenbergMarquardtOptions {
        damping: Some(0.01),
        max_iterations: Some(200),
        initial_values: vec![0., 100., 1., 0.1],
        timeout: None,
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        error_tolerance: None,
        central_difference: None,
        gradient_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };

    let res = lm(&data, &four_param_eq, &opts).unwrap();
    assert_eq!(res.iterations, 0);
    assert_vec_close_decimals(&res.parameter_values, &[0., 100., 1., 0.1], 3);
    assert_close_decimals(res.parameter_error, 19289.706, 3);
}

#[test]
fn ill_behaved_timeout() {
    let data = Data2D {
        x: vec![
            9.22e-12, 5.53e-11, 3.32e-10, 1.99e-9, 1.19e-8, 7.17e-8, 4.3e-7, 0.00000258, 0.0000155,
            0.0000929,
        ],
        y: vec![
            7.807, -3.74, 21.119, 2.382, 4.269, 41.57, 73.401, 98.535, 97.059, 92.147,
        ],
    };

    let opts = LevenbergMarquardtOptions {
        timeout: Some(0.0),
        damping: Some(0.00001),
        max_iterations: Some(200),
        initial_values: vec![0., 100., 1., 0.1],
        weights: None,
        damping_step_up: None,
        damping_step_down: None,
        error_tolerance: None,
        central_difference: None,
        gradient_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };

    let err = lm(&data, &four_param_eq, &opts).err().unwrap();
    assert!(
        err.contains("The execution time is over to 0 seconds"),
        "{}",
        err
    );
}

#[test]
fn real_world_problem_four_param_eq() {
    let data = Data2D {
        x: vec![
            9.22e-12, 5.53e-11, 3.32e-10, 1.99e-9, 1.19e-8, 7.17e-8, 4.3e-7, 0.00000258, 0.0000155,
            0.0000929,
        ],
        y: vec![
            7.807, -3.74, 21.119, 2.382, 4.269, 41.57, 73.401, 98.535, 97.059, 92.147,
        ],
    };

    let expected_iterations = 200usize;
    let expected_parameter_error = 16398.0009709;
    let expected_param_values = vec![-16.7697, 43.4549, 1018.8938, -4.3514];

    let opts = LevenbergMarquardtOptions {
        damping: Some(0.00001),
        max_iterations: Some(200),
        weights: Some(Weights::Scalar(1.0)),
        initial_values: vec![0., 100., 1., 0.1],
        timeout: None,
        damping_step_up: None,
        damping_step_down: None,
        error_tolerance: None,
        central_difference: None,
        gradient_difference: None,
        improvement_threshold: None,
        min_values: None,
        max_values: None,
    };

    let actual = lm(&data, &four_param_eq, &opts).unwrap();
    assert_eq!(actual.iterations, expected_iterations);
    assert_close_decimals(actual.parameter_error, expected_parameter_error, 3);
    assert_vec_close_decimals(&actual.parameter_values, &expected_param_values, 3);
}
