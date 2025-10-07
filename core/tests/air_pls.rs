use msut::utilities::air_pls::{AirPlsOptions, air_pls};

#[cfg(test)]
#[test]
fn small_vector_finds_flat_baseline() {
    let y = vec![1.0, 1.0, 1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0, 1.0];
    let x = (0..10).map(|i| i as f64).collect::<Vec<_>>();

    let res = air_pls(
        &x,
        &y,
        AirPlsOptions {
            lambda: 10.0,
            max_iterations: 100,
            tolerance: 1e-3,
            ..AirPlsOptions::default()
        },
    );

    let expected_baseline = vec![1.0; 10];
    let expected_corrected = vec![0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0];

    assert_vec_close(&res.baseline, &expected_baseline, 1e-2);
    assert_vec_close(&res.corrected, &expected_corrected, 1e-2);
}

fn assert_vec_close(a: &[f64], b: &[f64], tol: f64) {
    assert_eq!(
        a.len(),
        b.len(),
        "length mismatch ({} vs {})",
        a.len(),
        b.len()
    );
    for (i, (ai, bi)) in a.iter().zip(b.iter()).enumerate() {
        let diff = (ai - bi).abs();
        assert!(
            diff <= tol,
            "index {i}: got {ai}, expected {bi}, diff {diff} > tol {tol}"
        );
    }
}
