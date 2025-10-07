//! SAN plot (Signal–Artifact–Noise) — spectral quality analysis
//!
//! Rust implementation of JavaScript version from:
//! https://github.com/mljs/spectra-processing/blob/main/src/x/xNoiseSanPlot.ts
//!
//! # References
//! * Kirill F. Sheberstov, Eduard Sistaré Guardiola, Marion Pupier, Damien Jeannerat (2019).
//! “SAN plot: A graphical representation of the signal, noise, and artifacts content of spectra.”
//!   https://doi.org/10.1002/mrc.4882
//! * MLJS `xNoiseSanPlot.ts` `spectra-processing` repository.

use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct DataXY {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct FromTo {
    pub from: usize,
    pub to: usize,
}

#[derive(Clone, Debug)]
pub struct ConsiderList {
    pub from: f64,
    pub step: f64,
    pub to: f64,
}

#[derive(Clone, Debug)]
pub struct NoiseSanPlotOptions {
    pub mask: Option<Vec<f64>>,
    pub cut_off: Option<f64>,
    pub refine: bool,
    pub magnitude_mode: bool,
    pub scale_factor: f64,
    pub factor_std: f64,
    pub fix_offset: bool,
    pub log_base_y: f64,
    pub consider_list: Option<ConsiderList>,
    pub from_to: Option<HashMap<String, FromTo>>,
}

impl Default for NoiseSanPlotOptions {
    fn default() -> Self {
        Self {
            mask: None,
            cut_off: None,
            refine: true,
            magnitude_mode: false,
            scale_factor: 1.0,
            factor_std: 5.0,
            fix_offset: true,
            log_base_y: 2.0,
            consider_list: None,
            from_to: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct NoiseSanPlotResult {
    pub positive: f64,
    pub negative: f64,
    pub snr: f64,
    pub sanplot: HashMap<String, DataXY>,
}

pub fn noise_san_plot(array: &[f64], options: NoiseSanPlotOptions) -> NoiseSanPlotResult {
    let mask = options.mask.as_ref();
    let cut_off_opt = options.cut_off;
    let refine = options.refine;
    let magnitude_mode = options.magnitude_mode;
    let scale_factor = options.scale_factor;
    let factor_std = options.factor_std;
    let fix_offset = options.fix_offset;
    let log_base_y = options.log_base_y;
    let consider_list = options.consider_list.unwrap_or(ConsiderList {
        from: 0.5,
        step: 0.1,
        to: 0.9,
    });

    let mut input = prepare_data(array, scale_factor, mask);

    if fix_offset && !magnitude_mode {
        let m = input.len() / 2;
        let median = if input.len() % 2 == 0 {
            0.5 * (input[m - 1] + input[m])
        } else {
            input[m]
        };
        for v in &mut input {
            *v -= median;
        }
    }

    let first_negative_index = if *input.last().unwrap_or(&0.0) >= 0.0 {
        input.len()
    } else {
        match input.iter().position(|e| *e < 0.0) {
            Some(idx) => idx,
            None => input.len(),
        }
    };

    let mut last_positive_index: isize = first_negative_index as isize - 1;
    if last_positive_index >= 0 {
        let mut i = last_positive_index;
        while i >= 0 {
            if input[i as usize] > 0.0 {
                last_positive_index = i;
                break;
            }
            i -= 1;
        }
    }

    let sign_positive = if last_positive_index >= 0 {
        input[0..=last_positive_index as usize].to_vec()
    } else {
        Vec::new()
    };
    let sign_negative = if first_negative_index < input.len() {
        input[first_negative_index..].to_vec()
    } else {
        Vec::new()
    };

    let cut_off_dist = if let Some(c) = cut_off_opt {
        c
    } else {
        determine_cut_off(&sign_positive, magnitude_mode, &consider_list)
    };

    let p_index = (sign_positive.len() as f64 * cut_off_dist).floor() as usize;
    let mut noise_level_positive = *sign_positive
        .get(p_index.min(sign_positive.len().saturating_sub(1)))
        .unwrap_or(&0.0);

    let sky_point = *sign_positive.get(0).unwrap_or(&0.0);

    let mut initial_noise_level_negative = 0.0;
    if !sign_negative.is_empty() {
        let n_index = (sign_negative.len() as f64 * (1.0 - cut_off_dist)).floor() as usize;
        let pick = *sign_negative
            .get(n_index.min(sign_negative.len().saturating_sub(1)))
            .unwrap_or(&0.0);
        initial_noise_level_negative = -pick;
    }
    let mut noise_level_negative = initial_noise_level_negative;

    let mut clone_positive = sign_positive.clone();
    let mut clone_negative = sign_negative.clone();

    let mut cut_off_signals_index_plus: isize = 0;
    let mut cut_off_signals_index_neg: isize = 2;

    if refine {
        let cut_off_signals_pos = noise_level_positive * factor_std;
        let idx_plus = sign_positive
            .iter()
            .position(|e| *e < cut_off_signals_pos)
            .map(|v| v as isize)
            .unwrap_or(-1);
        cut_off_signals_index_plus = idx_plus;

        if idx_plus > -1 {
            let start = idx_plus as usize;
            clone_positive = sign_positive[start..].to_vec();
            let p2 = (clone_positive.len() as f64 * cut_off_dist).floor() as usize;
            noise_level_positive = *clone_positive
                .get(p2.min(clone_positive.len().saturating_sub(1)))
                .unwrap_or(&noise_level_positive);
        }

        let cut_off_signals_neg = noise_level_negative * factor_std;
        let idx_neg = sign_negative
            .iter()
            .position(|e| *e < cut_off_signals_neg)
            .map(|v| v as isize)
            .unwrap_or(-1);
        cut_off_signals_index_neg = idx_neg;

        if idx_neg > -1 {
            let start = idx_neg as usize;
            clone_negative = sign_negative[start..].to_vec();
            let n2 = (clone_negative.len() as f64 * (1.0 - cut_off_dist)).floor() as usize;
            if !clone_positive.is_empty() {
                let pick = *clone_positive
                    .get(n2.min(clone_positive.len().saturating_sub(1)))
                    .unwrap_or(&noise_level_negative);
                noise_level_negative = pick;
            }
        }
    }

    let correction_factor = -simple_norm_inv_number(cut_off_dist / 2.0, magnitude_mode);
    let mut effective_cut_off_dist;
    let mut refined_correction_factor;

    if refine && cut_off_signals_index_plus > -1 {
        effective_cut_off_dist = (cut_off_dist * clone_positive.len() as f64
            + cut_off_signals_index_plus as f64)
            / (clone_positive.len() as f64 + cut_off_signals_index_plus as f64);
        refined_correction_factor =
            -simple_norm_inv_number(effective_cut_off_dist / 2.0, magnitude_mode);
        noise_level_positive /= refined_correction_factor;

        if cut_off_signals_index_neg > -1 {
            effective_cut_off_dist = (cut_off_dist * clone_negative.len() as f64
                + cut_off_signals_index_neg as f64)
                / (clone_negative.len() as f64 + cut_off_signals_index_neg as f64);
            refined_correction_factor =
                -simple_norm_inv_number(effective_cut_off_dist / 2.0, magnitude_mode);
            if noise_level_negative != 0.0 {
                noise_level_negative /= refined_correction_factor;
            }
        }
    } else {
        noise_level_positive /= correction_factor;
        noise_level_negative /= correction_factor;
    }

    let mut sanplot = HashMap::new();
    let mut from_to = HashMap::new();
    from_to.insert(
        "positive".to_string(),
        FromTo {
            from: 0,
            to: if last_positive_index >= 0 {
                last_positive_index as usize
            } else {
                0
            },
        },
    );
    from_to.insert(
        "negative".to_string(),
        FromTo {
            from: first_negative_index,
            to: input.len(),
        },
    );
    let scaled = generate_san_plot(&input, Some(from_to), log_base_y);
    for (k, v) in scaled {
        sanplot.insert(k, v);
    }

    NoiseSanPlotResult {
        positive: noise_level_positive,
        negative: noise_level_negative,
        snr: if noise_level_positive == 0.0 {
            0.0
        } else {
            sky_point / noise_level_positive
        },
        sanplot,
    }
}

fn prepare_data(array: &[f64], scale_factor: f64, mask: Option<&Vec<f64>>) -> Vec<f64> {
    let mut v = Vec::with_capacity(array.len());
    if let Some(m) = mask {
        if m.len() == array.len() {
            for (i, &val) in array.iter().enumerate() {
                if m[i] == 0.0 {
                    v.push(val);
                }
            }
        } else {
            v.extend_from_slice(array);
        }
    } else {
        v.extend_from_slice(array);
    }
    if scale_factor > 1.0 {
        for x in &mut v {
            *x *= scale_factor;
        }
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v.reverse();
    v
}

fn determine_cut_off(sign_positive: &[f64], magnitude_mode: bool, consider: &ConsiderList) -> f64 {
    let mut cut_off_pairs: Vec<(f64, f64)> = Vec::new();
    let index_max = sign_positive.len().saturating_sub(1) as f64;
    let mut i = 0.01;
    while i <= 0.99 + 1e-12 {
        let idx = (index_max * i).round() as usize;
        let val = -sign_positive
            .get(idx.min(sign_positive.len().saturating_sub(1)))
            .copied()
            .unwrap_or(0.0)
            / simple_norm_inv_number(i / 2.0, magnitude_mode);
        cut_off_pairs.push((i, val));
        i += 0.01;
    }

    let mut min_ki = f64::MAX;
    let delta = consider.step / 2.0;
    let mut where_to_cut = 0.5;
    let mut t = consider.from;
    while t <= consider.to + 1e-12 {
        let floor = t - delta;
        let top = t + delta;
        let window: Vec<f64> = cut_off_pairs
            .iter()
            .filter(|(p, _)| *p < top && *p > floor)
            .map(|(_, v)| *v)
            .collect();
        let mut average_value = 0.0;
        for val in &window {
            average_value += val.abs();
        }
        let mut ki_sqrt = 0.0;
        for val in &window {
            ki_sqrt += (val - average_value) * (val - average_value);
        }
        if ki_sqrt < min_ki {
            min_ki = ki_sqrt;
            where_to_cut = t;
        }
        t += consider.step;
    }
    where_to_cut
}

fn generate_san_plot(
    array: &[f64],
    from_to: Option<HashMap<String, FromTo>>,
    log_base_y: f64,
) -> HashMap<String, DataXY> {
    let mut out = HashMap::new();
    if let Some(ranges) = from_to {
        for (key, rng) in ranges {
            let segment = if rng.from != rng.to {
                scale_slice(&array[rng.from..rng.to], log_base_y)
            } else {
                DataXY {
                    x: Vec::new(),
                    y: Vec::new(),
                }
            };
            if key == "negative" {
                let mut flipped = segment.clone();
                flipped.y.reverse();
                out.insert(key, flipped);
            } else {
                out.insert(key, segment);
            }
        }
    }
    out
}

fn scale_slice(slice: &[f64], log_base_y: f64) -> DataXY {
    let mut y = slice.to_vec();
    if log_base_y != 0.0 {
        let log_of_base = log_base_y.log10();
        for v in &mut y {
            *v = v.abs().log10() / log_of_base;
        }
    }
    let mut x = Vec::with_capacity(y.len());
    for i in 0..y.len() {
        x.push(i as f64);
    }
    DataXY { x, y }
}

pub fn simple_norm_inv(data: &[f64], magnitude_mode: bool) -> Vec<f64> {
    let mut out = vec![0.0; data.len()];
    if magnitude_mode {
        for (i, &v) in data.iter().enumerate() {
            out[i] = -(-2.0 * (1.0 - v).ln()).sqrt();
        }
    } else {
        for (i, &v) in data.iter().enumerate() {
            out[i] = -std::f64::consts::SQRT_2 * erfcinv(2.0 * v);
        }
    }
    out
}

pub fn simple_norm_inv_number(x: f64, magnitude_mode: bool) -> f64 {
    simple_norm_inv(&[x], magnitude_mode)[0]
}

fn polyval(c: &[f64], x: f64) -> f64 {
    let mut p = 0.0;
    for &coef in c {
        p = p * x + coef;
    }
    p
}

fn calc(x: f64, v: f64, p: &[f64], q: &[f64], y: f64) -> f64 {
    let s = x - v;
    let r = polyval(p, s) / polyval(q, s);
    y * x + r * x
}

pub fn erfcinv(mut x: f64) -> f64 {
    const Y1: f64 = 8.913_147_449_493_408e-2;
    const P1: [f64; 8] = [
        -5.387_729_650_712_429e-3,
        8.226_878_746_769_158e-3,
        2.198_786_811_111_689e-2,
        -3.656_379_714_117_626e-2,
        -1.269_261_476_629_740e-2,
        3.348_066_254_097_446e-2,
        -8.368_748_197_417_368e-3,
        -5.087_819_496_582_806e-4,
    ];
    const Q1: [f64; 9] = [
        8.862_163_904_564_247e-4,
        -2.333_937_593_741_900e-3,
        7.952_836_873_415_717e-2,
        -5.273_963_823_400_997e-2,
        -7.122_890_234_154_285e-1,
        6.623_288_404_720_030e-1,
        1.562_215_583_984_230,
        -1.565_745_582_341_758,
        -9.700_050_433_032_906e-1,
    ];

    const Y2: f64 = 2.249_481_201_171_875;
    const P2: [f64; 9] = [
        -3.671_922_547_077_293,
        2.112_946_554_483_405e1,
        1.744_538_598_557_086e1,
        -4.463_823_244_417_870e1,
        -1.885_106_480_587_142e1,
        1.764_472_984_083_740e1,
        8.370_503_283_431_199,
        1.052_646_806_993_917e-1,
        -2.024_335_083_559_388e-1,
    ];
    const Q2: [f64; 9] = [
        1.721_147_657_612_003,
        -2.264_369_334_131_397e1,
        1.082_686_673_554_602e1,
        4.856_092_131_087_399e1,
        -2.014_326_346_804_852e1,
        -2.866_081_804_998_000e1,
        3.971_343_795_334_387,
        6.242_641_248_542_475,
        1.0,
    ];

    const Y3: f64 = 8.072_204_589_843_75e-1;
    const P3: [f64; 11] = [
        -6.811_499_568_537_770e-10,
        2.852_253_317_822_170e-8,
        -6.794_655_751_811_264e-7,
        2.145_589_953_888_053e-3,
        2.901_579_100_053_291e-2,
        1.428_695_344_081_571e-1,
        3.377_855_389_120_359e-1,
        3.870_797_389_726_043e-1,
        1.170_301_563_419_953e-1,
        -1.637_940_471_933_171e-1,
        -1.311_027_816_799_519e-1,
    ];
    const Q3: [f64; 8] = [
        1.105_924_229_346_489e-2,
        1.522_643_382_953_318e-1,
        8.488_543_434_579_020e-1,
        2.593_019_216_236_203,
        4.778_465_929_458_438,
        5.381_683_457_070_069,
        3.466_254_072_425_672,
        1.0,
    ];

    const Y4: f64 = 9.399_557_113_647_461e-1;
    const P4: [f64; 9] = [
        2.663_392_274_257_820e-12,
        -2.304_047_769_118_826e-10,
        4.604_698_905_843_180e-6,
        1.575_446_174_249_606e-4,
        1.871_234_928_195_592e-3,
        9.508_047_013_259_196e-3,
        1.855_733_065_142_311e-2,
        -2.224_265_292_134_479e-3,
        -3.503_537_871_831_780e-2,
    ];
    const Q4: [f64; 7] = [
        7.646_752_923_027_945e-5,
        2.638_616_766_570_160e-3,
        3.415_891_436_709_477e-2,
        2.200_911_057_641_312e-1,
        7.620_591_645_536_234e-1,
        1.365_334_981_755_406,
        1.0,
    ];

    const Y5: f64 = 9.836_282_730_102_539e-1;
    const P5: [f64; 9] = [
        9.905_570_997_331_033e-17,
        -2.811_287_356_288_318e-14,
        4.625_961_635_228_786e-9,
        4.496_967_899_277_065e-7,
        1.496_247_837_583_424e-5,
        2.093_863_174_875_881e-4,
        1.056_288_621_524_929e-3,
        -1.129_514_387_455_803e-3,
        -1.674_310_050_766_337e-2,
    ];
    const Q5: [f64; 7] = [
        2.822_431_720_161_080e-7,
        2.753_354_747_647_260e-5,
        9.640_118_070_051_655e-4,
        1.607_460_870_936_765e-2,
        1.381_518_657_490_833e-1,
        5.914_293_448_864_175e-1,
        1.0,
    ];

    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 || x > 2.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::INFINITY;
    }
    if x == 2.0 {
        return f64::NEG_INFINITY;
    }
    if x == 1.0 {
        return 0.0;
    }

    let mut sign = false;
    let mut q;
    if x > 1.0 {
        q = 2.0 - x;
        x = 1.0 - q;
        sign = true;
    } else {
        q = x;
        x = 1.0 - x;
    }

    if x <= 0.5 {
        let g = x * (x + 10.0);
        let r = polyval(&P1, x) / polyval(&Q1, x);
        let val = g * Y1 + g * r;
        return if sign { -val } else { val };
    }

    if q >= 0.25 {
        let g = (-2.0 * q.ln()).sqrt();
        let q2 = q - 0.25;
        let r = polyval(&P2, q2) / polyval(&Q2, q2);
        let val = g / (Y2 + r);
        return if sign { -val } else { val };
    }

    q = -q.ln().sqrt();
    if q < 3.0 {
        return calc(q, 1.125, &P3, &Q3, Y3);
    }
    if q < 6.0 {
        return calc(q, 3.0, &P4, &Q4, Y4);
    }
    calc(q, 6.0, &P5, &Q5, Y5)
}
