use crate::utilities::{lm::ParameterizedFunction, utilities::erfc_approx};

use std::f64::consts::{LN_2, PI, SQRT_2};

/// ## Gaussian
/// Wikipedia https://en.wikipedia.org/wiki/Gaussian_function
/// ### Function
/// f(rt) = h * exp( - (rt-μ)^2 / (2σ^2))
pub fn gaussian_fn(rt: f64, h: f64, mu: f64, sigma: f64) -> f64 {
    let sigma = sigma.max(1e-12);
    let d = rt - mu;
    h * (-(d * d) / (2.0 * sigma * sigma)).exp()
}

/// ## Exponentially modified Gaussian distribution
/// https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
/// ### Function
/// f(rt; h, μ, σ, τ) = (h σ / τ) * sqrt(π/2) * exp(0.5*(σ/τ)^2 - (rt - μ)/τ) * erfc( (1/√2) * (σ/τ - (rt - μ)/σ ))
pub fn emg_fn(rt: f64, h: f64, mu: f64, sigma: f64, tau: f64) -> f64 {
    let sigma = sigma.max(1e-12);
    let tau = tau.max(1e-12);
    let a = sigma / tau;
    let z = (a - (rt - mu) / sigma) / SQRT_2;
    let g = (-0.5 * ((rt - mu) / sigma) * ((rt - mu) / sigma)).exp();
    h * a * (PI / 2.0).sqrt() * g * erfcx(z)
}

#[inline]
fn erfcx(z: f64) -> f64 {
    (z * z).exp() * erfc_approx(z)
}

#[inline]
pub fn sigma_from_fwhm(fwhm: f64) -> f64 {
    // fwhm = 2 * sqrt(2 ln 2) * sigma
    fwhm / (2.0 * (2.0 * LN_2).sqrt())
}

pub fn gaussian(mu: f64) -> Box<ParameterizedFunction> {
    Box::new(move |params: &[f64]| {
        let h = params.get(0).copied().unwrap_or(1.0);
        let fwhm = params.get(1).copied().unwrap_or(1.0).abs().max(1e-12);
        let sigma = sigma_from_fwhm(fwhm);
        Box::new(move |rt: f64| gaussian_fn(rt, h, mu, sigma))
    })
}

pub fn emg(mu: f64) -> Box<ParameterizedFunction> {
    Box::new(move |params: &[f64]| {
        let h_at_mu = params.get(0).copied().unwrap_or(1.0);
        let fwhm = params.get(1).copied().unwrap_or(1.0).abs().max(1e-12);
        let tau = params.get(2).copied().unwrap_or(1.0).abs().max(1e-12);

        let sigma = sigma_from_fwhm(fwhm);
        let a = sigma / tau;

        let unit_at_mu = a * (PI / 2.0).sqrt() * erfcx(a / SQRT_2);
        let h_amp = h_at_mu / unit_at_mu.max(1e-24);

        Box::new(move |rt: f64| emg_fn(rt, h_amp, mu, sigma, tau))
    })
}
