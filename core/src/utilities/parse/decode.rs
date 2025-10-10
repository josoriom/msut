use serde_json;

use crate::utilities::parse::{
    helper::{rd_f64, rd_u32, rd_u64, read_array_as_f64},
    parse_mzml::{ChromatogramSummary, MzML, Precursor, Run, SpectrumSummary},
};

pub fn decode(bin: &[u8]) -> Result<MzML, String> {
    if bin.len() < 64 {
        return Err("short header".into());
    }
    let magic = &bin[0..4];

    let n_spec = rd_u32(bin, 4)? as usize;
    let n_ch = rd_u32(bin, 8)? as usize;
    let cx = bin[12];
    let cy = bin[13];
    let sx = bin[14];
    let sy = bin[15];

    let s_idx_off = rd_u64(bin, 16)? as usize;
    let c_idx_off = rd_u64(bin, 24)? as usize;
    let s_meta_off = rd_u64(bin, 32)? as usize;
    let c_meta_off = rd_u64(bin, 40)? as usize;
    let data_off = rd_u64(bin, 48)? as usize;
    let _total = rd_u64(bin, 56)? as usize;
    if data_off > bin.len() {
        return Err("data_off OOB".into());
    }

    let s_idx_len = n_spec.checked_mul(32).ok_or("spec index ovf")?;
    if s_idx_off + s_idx_len > bin.len() {
        return Err("spec index OOB".into());
    }
    let mut sidx: Vec<(u64, u32, u64, u32)> = Vec::with_capacity(n_spec);
    for i in 0..n_spec {
        let b = s_idx_off + i * 32;
        sidx.push((
            rd_u64(bin, b + 0)?,
            rd_u32(bin, b + 8)?,
            rd_u64(bin, b + 12)?,
            rd_u32(bin, b + 20)?,
        ));
    }

    let c_idx_len = n_ch.checked_mul(32).ok_or("chrom index ovf")?;
    if c_idx_off + c_idx_len > bin.len() {
        return Err("chrom index OOB".into());
    }
    let mut cidx: Vec<(u64, u32, u64, u32)> = Vec::with_capacity(n_ch);
    for i in 0..n_ch {
        let b = c_idx_off + i * 32;
        cidx.push((
            rd_u64(bin, b + 0)?,
            rd_u32(bin, b + 8)?,
            rd_u64(bin, b + 12)?,
            rd_u32(bin, b + 20)?,
        ));
    }

    let mut spectra: Vec<SpectrumSummary> = Vec::with_capacity(n_spec);
    let mut chroms: Vec<ChromatogramSummary> = Vec::with_capacity(n_ch);

    if magic == b"BIN1" {
        let s_meta_len = n_spec.checked_mul(104).ok_or("spec meta ovf")?;
        if s_meta_off + s_meta_len > bin.len() {
            return Err("spec meta OOB".into());
        }
        for i in 0..n_spec {
            let b = s_meta_off + i * 104;
            let index = rd_u32(bin, b + 0)? as usize;
            let array_length = rd_u32(bin, b + 4)? as usize;
            let ms_level = {
                let v = bin[b + 8];
                if v == 255 { None } else { Some(v) }
            };
            let polarity = {
                let v = bin[b + 9];
                if v == 255 { None } else { Some(v) }
            };
            let spectrum_type = {
                let v = bin[b + 10];
                if v == 255 { None } else { Some(v) }
            };
            let rt = {
                let v = rd_f64(bin, b + 12)?;
                if v < 0.0 { None } else { Some(v) }
            };
            let swl = {
                let v = rd_f64(bin, b + 20)?;
                if v < 0.0 { None } else { Some(v) }
            };
            let swu = {
                let v = rd_f64(bin, b + 28)?;
                if v < 0.0 { None } else { Some(v) }
            };
            let tic = {
                let v = rd_f64(bin, b + 36)?;
                if v < 0.0 { None } else { Some(v) }
            };
            let bpi = {
                let v = rd_f64(bin, b + 44)?;
                if v < 0.0 { None } else { Some(v) }
            };
            let bpm = {
                let v = rd_f64(bin, b + 52)?;
                if v < 0.0 { None } else { Some(v) }
            };
            let pt = rd_f64(bin, b + 60)?;
            let pl = rd_f64(bin, b + 68)?;
            let pu = rd_f64(bin, b + 76)?;
            let ps = rd_f64(bin, b + 84)?;
            let prec = {
                let a = if pt < 0.0 { None } else { Some(pt) };
                let d = if pl < 0.0 { None } else { Some(pl) };
                let e = if pu < 0.0 { None } else { Some(pu) };
                let f = if ps < 0.0 { None } else { Some(ps) };
                if a.is_none() && d.is_none() && e.is_none() && f.is_none() {
                    None
                } else {
                    Some(Precursor {
                        isolation_window_target_mz: a,
                        isolation_window_lower_offset: d,
                        isolation_window_upper_offset: e,
                        selected_ion_mz: f,
                    })
                }
            };
            spectra.push(SpectrumSummary {
                index,
                array_length,
                ms_level,
                polarity,
                spectrum_type,
                retention_time: rt,
                scan_window_lower_limit: swl,
                scan_window_upper_limit: swu,
                total_ion_current: tic,
                base_peak_intensity: bpi,
                base_peak_mz: bpm,
                mz_array: None,
                intensity_array: None,
                precursor: prec,
            });
        }

        let c_meta_len = n_ch.checked_mul(24).ok_or("chrom meta ovf")?;
        if c_meta_off + c_meta_len > bin.len() {
            return Err("chrom meta OOB".into());
        }
        let mut ids: Vec<(u64, u32)> = Vec::with_capacity(n_ch);
        for i in 0..n_ch {
            let b = c_meta_off + i * 24;
            let index = rd_u32(bin, b + 0)? as usize;
            let array_length = rd_u32(bin, b + 4)? as usize;
            let off = rd_u64(bin, b + 8)?;
            let len = rd_u32(bin, b + 16)?;
            ids.push((off, len));
            chroms.push(ChromatogramSummary {
                index,
                array_length,
                time_array: None,
                intensity_array: None,
                id: String::new(),
            });
        }

        for (i, (x_off, x_len, _, _)) in sidx.iter().enumerate() {
            spectra[i].mz_array = read_array_as_f64(bin, *x_off, *x_len, sx)?;
        }
        for (i, (_, _, y_off, y_len)) in sidx.iter().enumerate() {
            spectra[i].intensity_array = read_array_as_f64(bin, *y_off, *y_len, sy)?;
        }
        for (i, (x_off, x_len, _, _)) in cidx.iter().enumerate() {
            chroms[i].time_array = read_array_as_f64(bin, *x_off, *x_len, cx)?;
        }
        for (i, (_, _, y_off, y_len)) in cidx.iter().enumerate() {
            chroms[i].intensity_array = read_array_as_f64(bin, *y_off, *y_len, cy)?;
        }

        for (i, (off, len)) in ids.into_iter().enumerate() {
            if off == 0 || len == 0 {
                chroms[i].id = String::new();
            } else {
                let o = off as usize;
                let l = len as usize;
                if o + l > bin.len() {
                    return Err("chrom id OOB".into());
                }
                let s = std::str::from_utf8(&bin[o..o + l]).unwrap_or_default();
                chroms[i].id = s.to_owned();
            }
        }
    } else if magic == b"BINS" {
        for i in 0..n_spec {
            spectra.push(SpectrumSummary {
                index: i,
                array_length: 0,
                ms_level: None,
                polarity: None,
                spectrum_type: None,
                retention_time: None,
                scan_window_lower_limit: None,
                scan_window_upper_limit: None,
                total_ion_current: None,
                base_peak_intensity: None,
                base_peak_mz: None,
                mz_array: None,
                intensity_array: None,
                precursor: None,
            });
        }
        for i in 0..n_ch {
            chroms.push(ChromatogramSummary {
                index: i,
                array_length: 0,
                time_array: None,
                intensity_array: None,
                id: String::new(),
            });
        }
        for (i, (x_off, x_len, _, _)) in sidx.iter().enumerate() {
            let a = read_array_as_f64(bin, *x_off, *x_len, sx)?;
            let n = a.as_ref().map(|v| v.len()).unwrap_or(0);
            spectra[i].mz_array = a;
            spectra[i].array_length = n.max(spectra[i].array_length);
        }
        for (i, (_, _, y_off, y_len)) in sidx.iter().enumerate() {
            let a = read_array_as_f64(bin, *y_off, *y_len, sy)?;
            let n = a.as_ref().map(|v| v.len()).unwrap_or(0);
            spectra[i].intensity_array = a;
            spectra[i].array_length = n.max(spectra[i].array_length);
        }
        for (i, (x_off, x_len, _, _)) in cidx.iter().enumerate() {
            let a = read_array_as_f64(bin, *x_off, *x_len, cx)?;
            let n = a.as_ref().map(|v| v.len()).unwrap_or(0);
            chroms[i].time_array = a;
            chroms[i].array_length = n.max(chroms[i].array_length);
        }
        for (i, (_, _, y_off, y_len)) in cidx.iter().enumerate() {
            let a = read_array_as_f64(bin, *y_off, *y_len, cy)?;
            let n = a.as_ref().map(|v| v.len()).unwrap_or(0);
            chroms[i].intensity_array = a;
            chroms[i].array_length = n.max(chroms[i].array_length);
        }
    } else {
        return Err("bad magic".into());
    }

    let run = Run {
        id: String::new(),
        start_time_stamp: None,
        default_instrument_configuration_ref: None,
        spectrum_list_count: Some(n_spec),
        chromatogram_list_count: Some(n_ch),
        spectra,
        chromatograms: chroms,
    };

    Ok(MzML {
        cv_list: Vec::new(),
        file_description: None,
        referenceable_param_groups: Vec::new(),
        sample_list: Vec::new(),
        instrument_configurations: Vec::new(),
        software_list: Vec::new(),
        data_processing_list: Vec::new(),
        acquisition_settings_list: Vec::new(),
        run: Some(run),
        index_list: None,
    })
}

pub fn metadata_to_json(src: &MzML) -> Result<Vec<u8>, serde_json::Error> {
    let run_meta = src.run.as_ref().map(|r| Run {
        id: r.id.clone(),
        start_time_stamp: r.start_time_stamp.clone(),
        default_instrument_configuration_ref: r.default_instrument_configuration_ref.clone(),
        spectrum_list_count: r.spectrum_list_count,
        chromatogram_list_count: r.chromatogram_list_count,
        spectra: r
            .spectra
            .iter()
            .map(|s| SpectrumSummary {
                index: s.index,
                array_length: s.array_length,
                ms_level: s.ms_level,
                polarity: s.polarity,
                spectrum_type: s.spectrum_type,
                retention_time: s.retention_time,
                scan_window_lower_limit: s.scan_window_lower_limit,
                scan_window_upper_limit: s.scan_window_upper_limit,
                total_ion_current: s.total_ion_current,
                base_peak_intensity: s.base_peak_intensity,
                base_peak_mz: s.base_peak_mz,
                mz_array: None,
                intensity_array: None,
                precursor: s.precursor.clone(),
            })
            .collect(),
        chromatograms: r
            .chromatograms
            .iter()
            .map(|c| ChromatogramSummary {
                index: c.index,
                array_length: c.array_length,
                time_array: None,
                intensity_array: None,
                id: c.id.clone(),
            })
            .collect(),
    });

    let meta_only = MzML {
        cv_list: src.cv_list.clone(),
        file_description: src.file_description.clone(),
        referenceable_param_groups: src.referenceable_param_groups.clone(),
        sample_list: src.sample_list.clone(),
        instrument_configurations: src.instrument_configurations.clone(),
        software_list: src.software_list.clone(),
        data_processing_list: src.data_processing_list.clone(),
        acquisition_settings_list: src.acquisition_settings_list.clone(),
        run: run_meta,
        index_list: src.index_list.clone(),
    };

    serde_json::to_vec(&meta_only)
}
