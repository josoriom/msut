use crate::utilities::parse::{
    helper::{
        ensure_cap, set_f64_at, set_u32_at, set_u64_at, write_f32_at, write_f32_le, write_f64_at,
        write_f64_le,
    },
    parse_mzml::MzML,
};

pub fn encode(mzml: &MzML) -> Vec<u8> {
    const H: usize = 64;
    const SI: usize = 32;
    const CI: usize = 32;
    const SM: usize = 104;
    const CM: usize = 24;

    #[inline]
    fn a8(x: usize) -> usize {
        (x + 7) & !7
    }

    let run = match mzml.run.as_ref() {
        Some(r) => r,
        None => {
            let mut out = vec![0u8; H];
            out[0..4].copy_from_slice(b"BIN1");
            set_u32_at(&mut out, 4, 0);
            set_u32_at(&mut out, 8, 0);
            out[12] = 2;
            out[13] = 1;
            out[14] = 2;
            out[15] = 1;
            set_u64_at(&mut out, 56, H as u64);
            return out;
        }
    };

    let n_spec = run.spectra.len() as u32;
    let n_ch = run.chromatograms.len() as u32;

    let sb = (n_spec as usize) * SI;
    let cb = (n_ch as usize) * CI;
    let smb = (n_spec as usize) * SM;
    let cmb = (n_ch as usize) * CM;

    let mut plan = H + sb + cb + smb + cmb;

    for s in &run.spectra {
        if let Some(v) = &s.mz_array {
            plan = a8(plan);
            plan += v.len() * 8;
        }
    }
    for s in &run.spectra {
        if let Some(v) = &s.intensity_array {
            plan = a8(plan);
            plan += v.len() * 4;
        }
    }
    for c in &run.chromatograms {
        if let Some(v) = &c.time_array {
            plan = a8(plan);
            plan += v.len() * 8;
        }
    }
    for c in &run.chromatograms {
        if let Some(v) = &c.intensity_array {
            plan = a8(plan);
            plan += v.len() * 4;
        }
    }
    for c in &run.chromatograms {
        let b = c.id.as_bytes();
        if !b.is_empty() {
            plan = a8(plan);
            plan += b.len();
        }
    }

    let mut out: Vec<u8> = Vec::with_capacity(plan);
    unsafe {
        out.set_len(plan);
    }
    let mut cur = H;

    out[0..4].copy_from_slice(b"BIN1");
    set_u32_at(&mut out, 4, n_spec);
    set_u32_at(&mut out, 8, n_ch);
    out[12] = 2;
    out[13] = 1;
    out[14] = 2;
    out[15] = 1;

    let spec_index_off = cur as u64;
    cur += sb;
    let chrom_index_off = cur as u64;
    cur += cb;
    let spec_meta_off = cur as u64;
    cur += smb;
    let chrom_meta_off = cur as u64;
    cur += cmb;
    let data_off = cur as u64;

    let mut sx: Vec<(u64, u32)> = Vec::with_capacity(n_spec as usize);
    for s in &run.spectra {
        let p = match &s.mz_array {
            Some(v) if !v.is_empty() => unsafe { write_f64_le(&mut out, &mut cur, v) },
            _ => (0, 0),
        };
        sx.push(p);
    }
    let mut sy: Vec<(u64, u32)> = Vec::with_capacity(n_spec as usize);
    for s in &run.spectra {
        let p = match &s.intensity_array {
            Some(v) if !v.is_empty() => unsafe { write_f32_le(&mut out, &mut cur, v) },
            _ => (0, 0),
        };
        sy.push(p);
    }
    let mut cx: Vec<(u64, u32)> = Vec::with_capacity(n_ch as usize);
    for c in &run.chromatograms {
        let p = match &c.time_array {
            Some(v) if !v.is_empty() => unsafe { write_f64_le(&mut out, &mut cur, v) },
            _ => (0, 0),
        };
        cx.push(p);
    }
    let mut cy: Vec<(u64, u32)> = Vec::with_capacity(n_ch as usize);
    for c in &run.chromatograms {
        let p = match &c.intensity_array {
            Some(v) if !v.is_empty() => unsafe { write_f32_le(&mut out, &mut cur, v) },
            _ => (0, 0),
        };
        cy.push(p);
    }
    let mut cid: Vec<(u64, u32)> = Vec::with_capacity(n_ch as usize);
    for c in &run.chromatograms {
        let s = c.id.as_bytes();
        if s.is_empty() {
            cid.push((0, 0));
        } else {
            cur = a8(cur);
            let off = cur as u64;
            let len = s.len() as u32;
            ensure_cap(&mut out, cur + s.len());
            unsafe {
                std::ptr::copy_nonoverlapping(s.as_ptr(), out.as_mut_ptr().add(cur), s.len());
            }
            cur += s.len();
            cid.push((off, len));
        }
    }

    for i in 0..(n_spec as usize) {
        let b = spec_index_off as usize + i * SI;
        let (x_off, x_len) = sx[i];
        let (y_off, y_len) = sy[i];
        set_u64_at(&mut out, b + 0, x_off);
        set_u32_at(&mut out, b + 8, x_len);
        set_u64_at(&mut out, b + 12, y_off);
        set_u32_at(&mut out, b + 20, y_len);
        set_u64_at(&mut out, b + 24, 0);
    }
    for i in 0..(n_ch as usize) {
        let b = chrom_index_off as usize + i * CI;
        let (x_off, x_len) = cx[i];
        let (y_off, y_len) = cy[i];
        set_u64_at(&mut out, b + 0, x_off);
        set_u32_at(&mut out, b + 8, x_len);
        set_u64_at(&mut out, b + 12, y_off);
        set_u32_at(&mut out, b + 20, y_len);
        set_u64_at(&mut out, b + 24, 0);
    }

    for (i, s) in run.spectra.iter().enumerate() {
        let b = spec_meta_off as usize + i * SM;
        set_u32_at(&mut out, b + 0, s.index as u32);
        set_u32_at(&mut out, b + 4, s.array_length as u32);
        out[b + 8] = s.ms_level.unwrap_or(255);
        out[b + 9] = s.polarity.unwrap_or(255);
        out[b + 10] = s.spectrum_type.unwrap_or(255);
        out[b + 11] = 0;
        set_f64_at(&mut out, b + 12, s.retention_time.unwrap_or(-1.0));
        set_f64_at(&mut out, b + 20, s.scan_window_lower_limit.unwrap_or(-1.0));
        set_f64_at(&mut out, b + 28, s.scan_window_upper_limit.unwrap_or(-1.0));
        set_f64_at(&mut out, b + 36, s.total_ion_current.unwrap_or(-1.0));
        set_f64_at(&mut out, b + 44, s.base_peak_intensity.unwrap_or(-1.0));
        set_f64_at(&mut out, b + 52, s.base_peak_mz.unwrap_or(-1.0));
        let (tgt, low, up, sel) = match &s.precursor {
            Some(p) => (
                p.isolation_window_target_mz.unwrap_or(-1.0),
                p.isolation_window_lower_offset.unwrap_or(-1.0),
                p.isolation_window_upper_offset.unwrap_or(-1.0),
                p.selected_ion_mz.unwrap_or(-1.0),
            ),
            None => (-1.0, -1.0, -1.0, -1.0),
        };
        set_f64_at(&mut out, b + 60, tgt);
        set_f64_at(&mut out, b + 68, low);
        set_f64_at(&mut out, b + 76, up);
        set_f64_at(&mut out, b + 84, sel);
    }

    for i in 0..(n_ch as usize) {
        let b = chrom_meta_off as usize + i * CM;
        set_u32_at(&mut out, b + 0, run.chromatograms[i].index as u32);
        set_u32_at(&mut out, b + 4, run.chromatograms[i].array_length as u32);
        let (off, len) = cid[i];
        set_u64_at(&mut out, b + 8, off);
        set_u32_at(&mut out, b + 16, len);
        set_u32_at(&mut out, b + 20, 0);
    }

    let total = cur as u64;
    set_u64_at(&mut out, 16, if n_spec > 0 { spec_index_off } else { 0 });
    set_u64_at(&mut out, 24, if n_ch > 0 { chrom_index_off } else { 0 });
    set_u64_at(&mut out, 32, if n_spec > 0 { spec_meta_off } else { 0 });
    set_u64_at(&mut out, 40, if n_ch > 0 { chrom_meta_off } else { 0 });
    set_u64_at(&mut out, 48, data_off);
    set_u64_at(&mut out, 56, total);

    out.truncate(cur);
    out
}

pub fn encode_arrays(mzml: &MzML) -> Vec<u8> {
    const H: usize = 64;
    const I: usize = 32;

    #[inline]
    fn a8(x: usize) -> usize {
        (x + 7) & !7
    }

    let run = match mzml.run.as_ref() {
        Some(r) => r,
        None => {
            let mut out = vec![0u8; H];
            out[0..4].copy_from_slice(b"BINS");
            set_u32_at(&mut out, 4, 1);
            out[12] = 2;
            out[13] = 1;
            out[14] = 2;
            out[15] = 1;
            set_u64_at(&mut out, 48, H as u64);
            set_u64_at(&mut out, 56, H as u64);
            return out;
        }
    };

    let n_spec = run.spectra.len() as u32;
    let n_ch = run.chromatograms.len() as u32;

    let sb = (n_spec as usize) * I;
    let cb = (n_ch as usize) * I;

    let mut cur = H + sb + cb;

    let mut smz: Vec<(u64, u32)> = Vec::with_capacity(n_spec as usize);
    for s in &run.spectra {
        if let Some(v) = &s.mz_array {
            cur = a8(cur);
            smz.push((cur as u64, v.len() as u32));
            cur += v.len() * 8;
        } else {
            smz.push((0, 0));
        }
    }
    let mut sin: Vec<(u64, u32)> = Vec::with_capacity(n_spec as usize);
    for s in &run.spectra {
        if let Some(v) = &s.intensity_array {
            cur = a8(cur);
            sin.push((cur as u64, v.len() as u32));
            cur += v.len() * 4;
        } else {
            sin.push((0, 0));
        }
    }
    let mut ctm: Vec<(u64, u32)> = Vec::with_capacity(n_ch as usize);
    for c in &run.chromatograms {
        if let Some(v) = &c.time_array {
            cur = a8(cur);
            ctm.push((cur as u64, v.len() as u32));
            cur += v.len() * 8;
        } else {
            ctm.push((0, 0));
        }
    }
    let mut cin: Vec<(u64, u32)> = Vec::with_capacity(n_ch as usize);
    for c in &run.chromatograms {
        if let Some(v) = &c.intensity_array {
            cur = a8(cur);
            cin.push((cur as u64, v.len() as u32));
            cur += v.len() * 4;
        } else {
            cin.push((0, 0));
        }
    }

    let total = cur;
    let mut out = vec![0u8; total];

    out[0..4].copy_from_slice(b"BINS");
    set_u32_at(&mut out, 4, n_spec);
    set_u32_at(&mut out, 8, n_ch);
    out[12] = 2;
    out[13] = 1;
    out[14] = 2;
    out[15] = 1;

    let s_idx_off = H as u64;
    let c_idx_off = (H + sb) as u64;
    let data_off = (H + sb + cb) as u64;

    for i in 0..(n_spec as usize) {
        let b = s_idx_off as usize + i * I;
        let (x_off, x_len) = smz[i];
        let (y_off, y_len) = sin[i];
        set_u64_at(&mut out, b + 0, x_off);
        set_u32_at(&mut out, b + 8, x_len);
        set_u64_at(&mut out, b + 12, y_off);
        set_u32_at(&mut out, b + 20, y_len);
    }
    for i in 0..(n_ch as usize) {
        let b = c_idx_off as usize + i * I;
        let (x_off, x_len) = ctm[i];
        let (y_off, y_len) = cin[i];
        set_u64_at(&mut out, b + 0, x_off);
        set_u32_at(&mut out, b + 8, x_len);
        set_u64_at(&mut out, b + 12, y_off);
        set_u32_at(&mut out, b + 20, y_len);
    }

    for (i, s) in run.spectra.iter().enumerate() {
        if let (Some(v), (off, len)) = (&s.mz_array, smz[i]) {
            if off != 0 && len != 0 {
                unsafe { write_f64_at(&mut out, off as usize, v) }
            }
        }
    }
    for (i, s) in run.spectra.iter().enumerate() {
        if let (Some(v), (off, len)) = (&s.intensity_array, sin[i]) {
            if off != 0 && len != 0 {
                unsafe { write_f32_at(&mut out, off as usize, v) }
            }
        }
    }
    for (i, c) in run.chromatograms.iter().enumerate() {
        if let (Some(v), (off, len)) = (&c.time_array, ctm[i]) {
            if off != 0 && len != 0 {
                unsafe { write_f64_at(&mut out, off as usize, v) }
            }
        }
    }
    for (i, c) in run.chromatograms.iter().enumerate() {
        if let (Some(v), (off, len)) = (&c.intensity_array, cin[i]) {
            if off != 0 && len != 0 {
                unsafe { write_f32_at(&mut out, off as usize, v) }
            }
        }
    }

    set_u64_at(&mut out, 16, s_idx_off);
    set_u64_at(&mut out, 24, c_idx_off);
    set_u64_at(&mut out, 32, 0);
    set_u64_at(&mut out, 40, 0);
    set_u64_at(&mut out, 48, data_off);
    set_u64_at(&mut out, 56, total as u64);

    out
}
