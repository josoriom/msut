use core::ffi::c_int;
use serde_json::json;
use std::{
    panic::{AssertUnwindSafe, catch_unwind},
    ptr, slice,
};

pub mod utilities;
use utilities::{
    calculate_eic::{EicOptions, calculate_eic_from_bin1},
    find_noise_level::find_noise_level as find_noise_level_rs,
    find_peaks::{FilterPeaksOptions, FindPeaksOptions, find_peaks as find_peaks_rs},
    get_boundaries::BoundariesOptions,
    get_peak::get_peak as get_peak_rs,
    get_peaks_from_chrom::get_peaks_from_chrom as get_peaks_from_chrom_rs,
    get_peaks_from_eic::get_peaks_from_eic as get_peaks_from_eic_rs,
    parse::{
        decode::{decode, metadata_to_json},
        encode::encode,
        parse_mzml::parse_mzml as parse_mzml_rs,
    },
    scan_for_peaks::ScanPeaksOptions,
    structs::{DataXY, FromTo, Roi},
};

use crate::utilities::structs::{ChromRoi, EicRoi};

const OK: c_int = 0;
const ERR_INVALID_ARGS: c_int = 1;
const ERR_PANIC: c_int = 2;
const ERR_PARSE: c_int = 4;
const EPS: f64 = 1e-5;

#[repr(C)]
pub struct Buf {
    pub ptr: *mut u8,
    pub len: usize,
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CPeakPOptions {
    pub integral_threshold: f64,
    pub intensity_threshold: f64,
    pub width_threshold: c_int,
    pub noise: f64,
    pub auto_noise: c_int,
    pub allow_overlap: c_int,
    pub window_size: c_int,
    pub sn_ratio: f64,
}

#[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
#[link(wasm_import_module = "env")]
unsafe extern "C" {
    fn js_log(ptr: *const u8, len: usize);
}

#[inline]
pub fn log_json<T: serde::Serialize>(v: &T) {
    if let Ok(s) = serde_json::to_string_pretty(v) {
        #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
        unsafe {
            js_log(s.as_ptr(), s.len());
        }

        #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
        eprintln!("{s}");
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn alloc(size: usize) -> *mut u8 {
    if size == 0 {
        return core::ptr::null_mut();
    }
    let mut v = Vec::<u8>::with_capacity(size);
    let p = v.as_mut_ptr();
    core::mem::forget(v);
    p
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn free_(ptr_raw: *mut u8, size: usize) {
    if !ptr_raw.is_null() {
        let _ = unsafe { Vec::<u8>::from_raw_parts(ptr_raw, size, size) };
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_mzml(
    data_ptr: *const u8,
    data_len: usize,
    out_data: *mut Buf,
) -> c_int {
    if data_ptr.is_null() || out_data.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let data = unsafe { slice::from_raw_parts(data_ptr, data_len) };
        let parsed = parse_mzml_rs(data, false).map_err(|_| ERR_PARSE)?;
        let bin = encode(&parsed);
        write_buf(out_data, bin.into_boxed_slice());
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn parse_mzml_to_json(
    data_ptr: *const u8,
    data_len: usize,
    out_json: *mut Buf,
    out_blob: *mut Buf,
) -> c_int {
    if data_ptr.is_null() || out_json.is_null() || out_blob.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let data = unsafe { slice::from_raw_parts(data_ptr, data_len) };
        let parsed = parse_mzml_rs(data, true).map_err(|_| ERR_PARSE)?;
        let meta = metadata_to_json(&parsed).map_err(|_| ERR_PARSE)?;
        let blob = encode(&parsed);

        write_buf(out_json, meta.into_boxed_slice());
        write_buf(out_blob, blob.into_boxed_slice());
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peak(
    x_ptr: *const f64,
    y_ptr: *const f32,
    len: usize,
    rt: f64,
    range: f64,
    options: *const CPeakPOptions,
    out_json: *mut Buf,
) -> c_int {
    if x_ptr.is_null() || y_ptr.is_null() || out_json.is_null() || len < 3 {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let xs = unsafe { slice::from_raw_parts(x_ptr, len) };
        let ys = unsafe { slice::from_raw_parts(y_ptr, len) };
        let data = DataXY {
            x: xs.to_vec(),
            y: ys.to_vec(),
        };

        let fp_opts = build_find_peaks_options(options);
        let roi = Roi { rt, window: range };
        let peak = get_peak_rs(&data, roi, Some(fp_opts));

        let s = match peak {
            Some(p) => serde_json::json!({
                "from": p.from,
                "to": p.to,
                "rt": p.rt,
                "integral": p.integral,
                "intensity": p.intensity,
                "ratio": p.ratio,
                "np": p.np
            })
            .to_string(),
            None => r#"{"from":0,"to":0,"rt":0,"integral":0,"intensity":0,"ratio":0,"np":0}"#
                .to_string(),
        };
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peaks_from_eic(
    bin_ptr: *const u8,
    bin_len: usize,
    rts_ptr: *const f64,
    mzs_ptr: *const f64,
    ranges_ptr: *const f64,
    ids_off_ptr: *const u32,
    ids_len_ptr: *const u32,
    ids_buf_ptr: *const u8,
    ids_buf_len: usize,
    n_items: usize,
    from_left: f64,
    to_right: f64,
    options: *const CPeakPOptions,
    cores: usize,
    out_json: *mut Buf,
) -> i32 {
    if bin_ptr.is_null()
        || rts_ptr.is_null()
        || mzs_ptr.is_null()
        || ranges_ptr.is_null()
        || out_json.is_null()
        || n_items == 0
    {
        return ERR_INVALID_ARGS;
    }
    let run = || -> Result<(), i32> {
        let bytes = unsafe { std::slice::from_raw_parts(bin_ptr, bin_len) };
        let rts = unsafe { std::slice::from_raw_parts(rts_ptr, n_items) };
        let mzs = unsafe { std::slice::from_raw_parts(mzs_ptr, n_items) };
        let ranges = unsafe { std::slice::from_raw_parts(ranges_ptr, n_items) };

        let has_ids = !(ids_off_ptr.is_null()
            || ids_len_ptr.is_null()
            || ids_buf_ptr.is_null()
            || ids_buf_len == 0);
        let (offs, lens, ibuf) = if has_ids {
            (
                unsafe { std::slice::from_raw_parts(ids_off_ptr, n_items) },
                unsafe { std::slice::from_raw_parts(ids_len_ptr, n_items) },
                Some(unsafe { std::slice::from_raw_parts(ids_buf_ptr, ids_buf_len) }),
            )
        } else {
            (&[][..], &[][..], None)
        };

        let mut items: Vec<EicRoi> = Vec::with_capacity(n_items);
        for i in 0..n_items {
            let rt = rts[i];
            let mz = mzs[i];
            let win = ranges[i];
            let ok = rt.is_finite() && mz.is_finite() && win.is_finite() && win > 0.0;

            let id = if let Some(buf) = ibuf {
                if has_ids {
                    let o = offs[i] as usize;
                    let l = lens[i] as usize;
                    if o.checked_add(l).map_or(true, |e| e > buf.len()) {
                        String::new()
                    } else {
                        std::str::from_utf8(&buf[o..o + l])
                            .unwrap_or("")
                            .to_string()
                    }
                } else {
                    String::new()
                }
            } else {
                String::new()
            };

            if ok {
                items.push(EicRoi {
                    id,
                    rt,
                    mz,
                    window: win,
                });
            } else {
                items.push(EicRoi {
                    id: String::new(),
                    rt: 0.0,
                    mz: 0.0,
                    window: 0.0,
                });
            }
        }

        let window = FromTo {
            from: from_left,
            to: to_right,
        };
        let fp = build_find_peaks_options(options);
        let peaks = get_peaks_from_eic_rs(bytes, window, items.as_slice(), Some(fp), cores)
            .ok_or(ERR_PARSE)?;

        let mut arr = Vec::with_capacity(peaks.len());
        for (id, ort, mz, p) in peaks {
            arr.push(serde_json::json!({
                "id": id,
                "mz": mz,
                "ort": ort,
                "rt": p.rt,
                "from": p.from,
                "to": p.to,
                "intensity": p.intensity,
                "integral": p.integral,
                "noise": p.noise
            }));
        }
        let s = serde_json::to_string(&arr).map_err(|_| ERR_PARSE)?;
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    };
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(run)) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn get_peaks_from_chrom(
    bin_ptr: *const u8,
    bin_len: usize,
    idxs_ptr: *const u32,
    rts_ptr: *const f64,
    ranges_ptr: *const f64,
    n_items: usize,
    options: *const CPeakPOptions,
    cores: usize,
    out_json: *mut Buf,
) -> i32 {
    if bin_ptr.is_null()
        || idxs_ptr.is_null()
        || rts_ptr.is_null()
        || ranges_ptr.is_null()
        || out_json.is_null()
        || n_items == 0
    {
        return ERR_INVALID_ARGS;
    }
    let run = || -> Result<(), i32> {
        let bin = unsafe { std::slice::from_raw_parts(bin_ptr, bin_len) };
        let idxs = unsafe { std::slice::from_raw_parts(idxs_ptr, n_items) };
        let rts = unsafe { std::slice::from_raw_parts(rts_ptr, n_items) };
        let wins = unsafe { std::slice::from_raw_parts(ranges_ptr, n_items) };
        let mzml = decode(bin).map_err(|_| ERR_PARSE)?;
        let chroms = &mzml.run.as_ref().ok_or(ERR_PARSE)?.chromatograms;

        let mut items = Vec::with_capacity(n_items);
        for i in 0..n_items {
            let iu = idxs[i];
            if iu == u32::MAX {
                items.push(ChromRoi {
                    id: String::new(),
                    idx: usize::MAX,
                    rt: 0.0,
                    window: 0.0,
                });
                continue;
            }
            let idx = iu as usize;
            if idx >= chroms.len() {
                items.push(ChromRoi {
                    id: String::new(),
                    idx,
                    rt: 0.0,
                    window: 0.0,
                });
                continue;
            }
            let id = chroms[idx].id.clone();
            items.push(ChromRoi {
                id,
                idx,
                rt: rts[i],
                window: wins[i],
            });
        }

        let fp = build_find_peaks_options(options);
        let list =
            get_peaks_from_chrom_rs(&mzml, items.as_slice(), Some(fp), cores).ok_or(ERR_PARSE)?;

        let mut out = Vec::with_capacity(list.len());
        for (index, id, ort, rt, from_, to_, intensity, integral) in list {
            out.push(serde_json::json!({
                "index": index,
                "id": id,
                "ort": ort,
                "rt": rt,
                "from": from_,
                "to": to_,
                "intensity": intensity,
                "integral":  integral
            }));
        }
        let s = serde_json::to_string(&out).map_err(|_| ERR_PARSE)?;
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    };
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(run)) {
        Ok(Ok(())) => OK,
        Ok(Err(c)) => c,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn find_peaks(
    x_ptr: *const f64,
    y_ptr: *const f32,
    len: usize,
    options: *const CPeakPOptions,
    out_json: *mut Buf,
) -> c_int {
    if x_ptr.is_null() || y_ptr.is_null() || out_json.is_null() {
        return ERR_INVALID_ARGS;
    }
    let run = || -> Result<(), c_int> {
        let xs = unsafe { slice::from_raw_parts(x_ptr, len) };
        let ys = unsafe { slice::from_raw_parts(y_ptr, len) };
        if xs.len() != ys.len() || xs.len() < 3 {
            return Err(ERR_INVALID_ARGS);
        }
        let data = DataXY {
            x: xs.to_vec(),
            y: ys.to_vec(),
        };
        let opts = build_find_peaks_options(options);
        let peaks = find_peaks_rs(&data, Some(opts));
        let list: Vec<_> = peaks
            .iter()
            .map(|p| {
                json!({
                    "from": p.from,
                    "to": p.to,
                    "rt": p.rt,
                    "integral": p.integral,
                    "intensity": p.intensity,
                    "ratio": p.ratio,
                    "np": p.np,
                    "noise": p.noise
                })
            })
            .collect();
        let s = serde_json::to_string(&list).map_err(|_| ERR_PARSE)?;
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    };
    match catch_unwind(AssertUnwindSafe(run)) {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn find_noise_level(y_ptr: *const f32, len: usize) -> f32 {
    if y_ptr.is_null() || len == 0 {
        return f32::INFINITY;
    }
    let compute = || {
        let ys = unsafe { slice::from_raw_parts(y_ptr, len) };
        find_noise_level_rs(ys)
    };
    match catch_unwind(AssertUnwindSafe(compute)) {
        Ok(noise) => noise as f32,
        Err(_) => f32::INFINITY,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn bin_to_json(
    bin_ptr: *const u8,
    bin_len: usize,
    out_json: *mut Buf,
) -> c_int {
    if bin_ptr.is_null() || out_json.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let bin = unsafe { slice::from_raw_parts(bin_ptr, bin_len) };
        let mzml = decode(bin);
        let s = serde_json::to_string(&mzml).map_err(|_| ERR_PARSE)?;
        write_buf(out_json, s.into_bytes().into_boxed_slice());
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn calculate_eic(
    bin_ptr: *const u8,
    bin_len: usize,
    target_ptr: *const u8,
    target_len: usize,
    from_time: f64,
    to_time: f64,
    ppm_tolerance: f64,
    mz_tolerance: f64,
    out_x: *mut Buf,
    out_y: *mut Buf,
) -> c_int {
    if bin_ptr.is_null() || target_ptr.is_null() || out_x.is_null() || out_y.is_null() {
        return ERR_INVALID_ARGS;
    }
    let res = catch_unwind(AssertUnwindSafe(|| -> Result<(), c_int> {
        let bin = unsafe { slice::from_raw_parts(bin_ptr, bin_len) };
        let target = std::str::from_utf8(unsafe { slice::from_raw_parts(target_ptr, target_len) })
            .map_err(|_| ERR_PARSE)?;

        let eic = calculate_eic_from_bin1(
            bin,
            target,
            FromTo {
                from: from_time,
                to: to_time,
            },
            EicOptions {
                ppm_tolerance,
                mz_tolerance,
            },
        )
        .map_err(|_| ERR_PARSE)?;

        let x_bytes = f64_slice_to_u8_box(&eic.x);
        let y_bytes = f32_slice_to_u8_box(&eic.y);
        write_buf(out_x, x_bytes);
        write_buf(out_y, y_bytes);
        Ok(())
    }));
    match res {
        Ok(Ok(())) => OK,
        Ok(Err(code)) => code,
        Err(_) => ERR_PANIC,
    }
}

fn f64_slice_to_u8_box(v: &[f64]) -> Box<[u8]> {
    let n = v.len() * 8;
    let mut out = Vec::<u8>::with_capacity(n);
    unsafe {
        out.set_len(n);
        ptr::copy_nonoverlapping(v.as_ptr() as *const u8, out.as_mut_ptr(), n);
    }
    out.into_boxed_slice()
}

fn f32_slice_to_u8_box(v: &[f32]) -> Box<[u8]> {
    let n = v.len() * 4;
    let mut out = Vec::<u8>::with_capacity(n);
    unsafe {
        out.set_len(n);
        ptr::copy_nonoverlapping(v.as_ptr() as *const u8, out.as_mut_ptr(), n);
    }
    out.into_boxed_slice()
}

#[inline]
fn odd_at_least(v: usize, min_: usize, def_: usize) -> usize {
    let v = if v == 0 { def_ } else { v };
    let v = v.max(min_);
    if v % 2 == 0 { v | 1 } else { v }
}

#[inline]
fn pos_usize(raw: c_int, def_: usize) -> usize {
    if raw > 0 { raw as usize } else { def_ }
}

fn write_buf(out: *mut Buf, bytes: Box<[u8]>) {
    let len = bytes.len();
    let ptr_bytes = Box::into_raw(bytes) as *mut u8;
    unsafe {
        ptr::write_unaligned(
            out,
            Buf {
                ptr: ptr_bytes,
                len,
            },
        )
    };
}

fn build_find_peaks_options(options: *const CPeakPOptions) -> FindPeaksOptions {
    if options.is_null() {
        let ws = odd_at_least(17, 5, 17);
        return FindPeaksOptions {
            scan_peaks_options: Some(ScanPeaksOptions {
                epsilon: EPS,
                window_size: ws,
            }),
            get_boundaries_options: Some(BoundariesOptions {
                ..Default::default()
            }),
            filter_peaks_options: None,
        };
    }
    let o = unsafe { *options };
    let ws = odd_at_least(pos_usize(o.window_size, 17), 5, 17);
    let integral = (o.integral_threshold.is_finite() && o.integral_threshold >= 0.0)
        .then_some(o.integral_threshold);
    let intensity = (o.intensity_threshold.is_finite() && o.intensity_threshold >= 0.0)
        .then_some(o.intensity_threshold);
    let width = (o.width_threshold > 0).then_some(o.width_threshold as usize);
    let noise = (o.noise.is_finite() && o.noise > 0.0).then_some(o.noise);
    let auto_noise = Some(o.auto_noise != 0);
    let allow_overlap = Some(o.allow_overlap != 0);
    let sn_ratio = if o.sn_ratio.is_finite() && o.sn_ratio > 0.0 {
        Some(o.sn_ratio)
    } else {
        Some(1.5) // DEefault value for s/n 
    };
    let filter = FilterPeaksOptions {
        integral_threshold: integral,
        intensity_threshold: intensity,
        width_threshold: width,
        noise,
        auto_noise,
        allow_overlap,
        sn_ratio,
    };
    FindPeaksOptions {
        scan_peaks_options: Some(ScanPeaksOptions {
            epsilon: EPS,
            window_size: ws,
        }),
        get_boundaries_options: Some(BoundariesOptions {
            ..Default::default()
        }),
        filter_peaks_options: Some(filter),
    }
}
