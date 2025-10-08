use base64::Engine;
use base64::engine::general_purpose::STANDARD;
use memchr::{memchr as mc_memchr, memmem};
use miniz_oxide::inflate::decompress_to_vec_zlib;
use serde::{Deserialize, Serialize};
use std::io::{Cursor, Read, Seek, SeekFrom};
use std::str;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectrumSummary {
    pub index: usize,
    pub array_length: usize,
    pub ms_level: Option<u8>,
    pub polarity: Option<u8>,
    pub spectrum_type: Option<u8>,
    pub retention_time: Option<f64>,
    pub scan_window_lower_limit: Option<f64>,
    pub scan_window_upper_limit: Option<f64>,
    pub total_ion_current: Option<f64>,
    pub base_peak_intensity: Option<f64>,
    pub base_peak_mz: Option<f64>,
    pub mz_array: Option<Vec<f64>>,
    pub intensity_array: Option<Vec<f32>>,
    pub precursor: Option<Precursor>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MzML {
    pub cv_list: Vec<CvEntry>,
    pub file_description: Option<FileDescription>,
    pub referenceable_param_groups: Vec<RefParamGroup>,
    pub sample_list: Vec<Sample>,
    pub instrument_configurations: Vec<InstrumentConfiguration>,
    pub software_list: Vec<Software>,
    pub data_processing_list: Vec<DataProcessing>,
    pub acquisition_settings_list: Vec<AcquisitionSettings>,
    pub run: Option<Run>,
    pub index_list: Option<IndexList>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CvEntry {
    pub id: String,
    pub full_name: Option<String>,
    pub version: Option<String>,
    pub uri: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileDescription {
    pub file_content: Vec<CvPair>,
    pub source_files: Vec<SourceFile>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CvPair {
    pub name: String,
    pub value: Option<String>,
    pub accession: Option<String>,
    pub cv_ref: Option<String>,
    pub unit_name: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SourceFile {
    pub id: String,
    pub name: Option<String>,
    pub location: Option<String>,
    pub cv_params: Vec<CvPair>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefParamGroup {
    pub id: String,
    pub cv_params: Vec<CvPair>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sample {
    pub id: String,
    pub name: Option<String>,
    pub cv_params: Vec<CvPair>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Software {
    pub id: String,
    pub version: Option<String>,
    pub cv_params: Vec<CvPair>,
    pub user_params: Vec<(String, Option<String>)>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Component {
    pub kind: String,
    pub order: Option<usize>,
    pub cv_params: Vec<CvPair>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InstrumentConfiguration {
    pub id: String,
    pub ref_param_group: Option<String>,
    pub cv_params: Vec<CvPair>,
    pub components: Vec<Component>,
    pub software_ref: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingMethod {
    pub order: Option<usize>,
    pub software_ref: Option<String>,
    pub cv_params: Vec<CvPair>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DataProcessing {
    pub id: String,
    pub methods: Vec<ProcessingMethod>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AcquisitionSettings {
    pub id: String,
    pub instrument_ref: Option<String>,
    pub cv_params: Vec<CvPair>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromatogramSummary {
    pub index: usize,
    pub array_length: usize,
    pub time_array: Option<Vec<f64>>,
    pub intensity_array: Option<Vec<f32>>,
    pub id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Run {
    pub id: String,
    pub start_time_stamp: Option<String>,
    pub default_instrument_configuration_ref: Option<String>,
    pub spectrum_list_count: Option<usize>,
    pub chromatogram_list_count: Option<usize>,
    pub spectra: Vec<SpectrumSummary>,
    pub chromatograms: Vec<ChromatogramSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexOffset {
    pub id_ref: Option<String>,
    pub offset: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexList {
    pub spectrum: Vec<IndexOffset>,
    pub chromatogram: Vec<IndexOffset>,
    pub index_list_offset: Option<u64>,
    pub file_checksum: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Precursor {
    pub isolation_window_target_mz: Option<f64>,
    pub isolation_window_lower_offset: Option<f64>,
    pub isolation_window_upper_offset: Option<f64>,
    pub selected_ion_mz: Option<f64>,
}

struct Scratch {
    b64_buf: Vec<u8>,
    zlib_buf: Vec<u8>,
}

pub fn parse_mzml(bytes: &[u8], slim: bool) -> Result<MzML, String> {
    if slim {
        let run_header = parse_run_header(bytes);
        let spectra = parse_spectra_internal(bytes)?;
        let chromatograms = parse_chromatograms_linear(bytes)?;
        let run = run_header.map(|mut r| {
            r.spectra = spectra;
            r.chromatograms = chromatograms;
            r
        });

        return Ok(MzML {
            cv_list: Vec::new(),
            file_description: None,
            referenceable_param_groups: Vec::new(),
            sample_list: Vec::new(),
            instrument_configurations: Vec::new(),
            software_list: Vec::new(),
            data_processing_list: Vec::new(),
            acquisition_settings_list: Vec::new(),
            run,
            index_list: None,
        });
    }

    let cv_list = parse_cv_list(bytes);
    let file_description = parse_file_description(bytes);
    let referenceable_param_groups = parse_ref_param_groups(bytes);
    let sample_list = parse_sample_list(bytes);
    let instrument_configurations = parse_instrument_configurations(bytes);
    let software_list = parse_software_list(bytes);
    let data_processing_list = parse_data_processing_list(bytes);
    let acquisition_settings_list = parse_acquisition_settings_list(bytes);
    let index_list = parse_index_list_wrapper(bytes);
    let run_header = parse_run_header(bytes);
    let spectra = parse_spectra_internal(bytes)?;
    let chromatograms = parse_chromatograms_linear(bytes)?;
    let run = run_header.map(|mut r| {
        r.spectra = spectra;
        r.chromatograms = chromatograms;
        r
    });

    Ok(MzML {
        cv_list,
        file_description,
        referenceable_param_groups,
        sample_list,
        instrument_configurations,
        software_list,
        data_processing_list,
        acquisition_settings_list,
        run,
        index_list,
    })
}

fn parse_cv_list(xml: &[u8]) -> Vec<CvEntry> {
    let mut out = Vec::new();
    if let Some(block) = find_section(xml, b"cvList") {
        let mut cur = 0usize;
        let f = memmem::Finder::new(b"<cv ");
        while let Some(p) = f.find(&block[cur..]) {
            let from = cur + p;
            let gt = match mc_memchr(b'>', &block[from..]) {
                Some(x) => from + x,
                None => break,
            };
            let head = &block[from..gt];
            let id = b2s(find_attr_value_in_tag(head, b"id")).unwrap_or_default();
            let full = b2s(find_attr_value_in_tag(head, b"fullName"));
            let ver = b2s(find_attr_value_in_tag(head, b"version"));
            let uri = b2s(find_attr_value_in_tag(head, b"URI"));
            out.push(CvEntry {
                id,
                full_name: full,
                version: ver,
                uri,
            });
            cur = gt + 1;
        }
    }
    out
}

fn parse_file_description(xml: &[u8]) -> Option<FileDescription> {
    let b = find_section(xml, b"fileDescription")?;
    let file_content = if let Some((s, e)) = tag_body(b, b"<fileContent>", b"</fileContent>") {
        collect_cv_params(&b[s..e])
    } else {
        Vec::new()
    };
    let mut source_files = Vec::new();
    if let Some((s, e)) = tag_body(b, b"<sourceFileList", b"</sourceFileList>") {
        let blk = &b[s..e];
        let mut cur = 0usize;
        let open_f = memmem::Finder::new(b"<sourceFile ");
        let close = b"</sourceFile>";
        let close_f = memmem::Finder::new(close);
        while let Some(p) = open_f.find(&blk[cur..]) {
            let from = cur + p;
            let gt = match mc_memchr(b'>', &blk[from..]) {
                Some(x) => from + x,
                None => break,
            };
            let head = &blk[from..gt];
            let id = b2s(find_attr_value_in_tag(head, b"id")).unwrap_or_default();
            let name = b2s(find_attr_value_in_tag(head, b"name"));
            let location = b2s(find_attr_value_in_tag(head, b"location"));
            let end_rel = match close_f.find(&blk[gt..]) {
                Some(v) => v,
                None => break,
            };
            let inner = &blk[gt..gt + end_rel];
            let cv_params = collect_cv_params(inner);
            source_files.push(SourceFile {
                id,
                name,
                location,
                cv_params,
            });
            cur = gt + end_rel + close.len();
        }
    }
    Some(FileDescription {
        file_content,
        source_files,
    })
}

fn parse_ref_param_groups(xml: &[u8]) -> Vec<RefParamGroup> {
    let mut out = Vec::new();
    if let Some(block) = find_section(xml, b"referenceableParamGroupList") {
        let mut cur = 0usize;
        let open_f = memmem::Finder::new(b"<referenceableParamGroup ");
        let close = b"</referenceableParamGroup>";
        let close_f = memmem::Finder::new(close);
        while let Some(p) = open_f.find(&block[cur..]) {
            let from = cur + p;
            let gt = match mc_memchr(b'>', &block[from..]) {
                Some(x) => from + x,
                None => break,
            };
            let head = &block[from..gt];
            let id = b2s(find_attr_value_in_tag(head, b"id")).unwrap_or_default();
            let end_rel = match close_f.find(&block[gt..]) {
                Some(v) => v,
                None => break,
            };
            let inner = &block[gt..gt + end_rel];
            let cv_params = collect_cv_params(inner);
            out.push(RefParamGroup { id, cv_params });
            cur = gt + end_rel + close.len();
        }
    }
    out
}

fn parse_sample_list(xml: &[u8]) -> Vec<Sample> {
    let mut out = Vec::new();
    if let Some(block) = find_section(xml, b"sampleList") {
        let mut cur = 0usize;
        let open_f = memmem::Finder::new(b"<sample ");
        let close = b"</sample>";
        let close_f = memmem::Finder::new(close);
        while let Some(p) = open_f.find(&block[cur..]) {
            let from = cur + p;
            let gt = match mc_memchr(b'>', &block[from..]) {
                Some(x) => from + x,
                None => break,
            };
            let head = &block[from..gt];
            let id = b2s(find_attr_value_in_tag(head, b"id")).unwrap_or_default();
            let name = b2s(find_attr_value_in_tag(head, b"name"));
            let end_rel = match close_f.find(&block[gt..]) {
                Some(v) => v,
                None => break,
            };
            let inner = &block[gt..gt + end_rel];
            let cv_params = collect_cv_params(inner);
            out.push(Sample {
                id,
                name,
                cv_params,
            });
            cur = gt + end_rel + close.len();
        }
    }
    out
}

fn parse_software_list(xml: &[u8]) -> Vec<Software> {
    let mut out = Vec::new();
    if let Some(block) = find_section(xml, b"softwareList") {
        let mut cur = 0usize;
        let open_f = memmem::Finder::new(b"<software ");
        let close = b"</software>";
        let close_f = memmem::Finder::new(close);
        while let Some(p) = open_f.find(&block[cur..]) {
            let from = cur + p;
            let gt = match mc_memchr(b'>', &block[from..]) {
                Some(x) => from + x,
                None => break,
            };
            let head = &block[from..gt];
            let id = b2s(find_attr_value_in_tag(head, b"id")).unwrap_or_default();
            let version = b2s(find_attr_value_in_tag(head, b"version"));
            let end_rel = match close_f.find(&block[gt..]) {
                Some(v) => v,
                None => break,
            };
            let inner = &block[gt..gt + end_rel];
            let cv_params = collect_cv_params(inner);
            let mut user_params = Vec::new();
            let mut pcur = 0usize;
            let up_open = memmem::Finder::new(b"<userParam ");
            while let Some(pp) = up_open.find(&inner[pcur..]) {
                let sf = pcur + pp;
                let sgt = match mc_memchr(b'>', &inner[sf..]) {
                    Some(x) => sf + x,
                    None => break,
                };
                let shead = &inner[sf..sgt];
                let name = b2s(find_attr_value_in_tag(shead, b"name")).unwrap_or_default();
                let val = b2s(find_attr_value_in_tag(shead, b"value"));
                user_params.push((name, val));
                pcur = sgt + 1;
            }
            out.push(Software {
                id,
                version,
                cv_params,
                user_params,
            });
            cur = gt + end_rel + close.len();
        }
    }
    out
}

fn parse_instrument_configurations(xml: &[u8]) -> Vec<InstrumentConfiguration> {
    let mut out = Vec::new();
    if let Some(block) = find_section(xml, b"instrumentConfigurationList") {
        let mut cur = 0usize;
        let open_f = memmem::Finder::new(b"<instrumentConfiguration ");
        let close = b"</instrumentConfiguration>";
        let close_f = memmem::Finder::new(close);
        while let Some(p) = open_f.find(&block[cur..]) {
            let from = cur + p;
            let gt = match mc_memchr(b'>', &block[from..]) {
                Some(x) => from + x,
                None => break,
            };
            let head = &block[from..gt];
            let id = b2s(find_attr_value_in_tag(head, b"id")).unwrap_or_default();
            let end_rel = match close_f.find(&block[gt..]) {
                Some(v) => v,
                None => break,
            };
            let inner = &block[gt..gt + end_rel];
            let ref_group = b2s(find_attr_value_in_tag(inner, b"ref"));
            let cv_params = collect_cv_params(inner);
            let mut components = Vec::new();
            for (kind, open_tag, close_tag) in [
                ("source", b"<source " as &[u8], b"</source>" as &[u8]),
                ("analyzer", b"<analyzer ", b"</analyzer>"),
                ("detector", b"<detector ", b"</detector>"),
            ] {
                let mut kcur = 0usize;
                let of = memmem::Finder::new(open_tag);
                let cf = memmem::Finder::new(close_tag);
                while let Some(pp) = of.find(&inner[kcur..]) {
                    let sf = kcur + pp;
                    let sgt = match mc_memchr(b'>', &inner[sf..]) {
                        Some(x) => sf + x,
                        None => break,
                    };
                    let shead = &inner[sf..sgt];
                    let order = find_attr_usize(shead, kind.as_bytes(), b"order");
                    let send_rel = match cf.find(&inner[sgt..]) {
                        Some(v) => v,
                        None => break,
                    };
                    let sinner = &inner[sgt..sgt + send_rel];
                    let scv = collect_cv_params(sinner);
                    components.push(Component {
                        kind: kind.to_string(),
                        order,
                        cv_params: scv,
                    });
                    kcur = sgt + send_rel + close_tag.len();
                }
            }
            let software_ref = if let Some(pp) = memmem::find(inner, b"<softwareRef ") {
                let sf = pp;
                let sgt = match mc_memchr(b'>', &inner[sf..]) {
                    Some(x) => sf + x,
                    None => sf,
                };
                b2s(find_attr_value_in_tag(&inner[sf..sgt], b"ref"))
            } else {
                None
            };
            out.push(InstrumentConfiguration {
                id,
                ref_param_group: ref_group,
                cv_params,
                components,
                software_ref,
            });
            cur = gt + end_rel + close.len();
        }
    }
    out
}

fn parse_data_processing_list(xml: &[u8]) -> Vec<DataProcessing> {
    let mut out = Vec::new();
    if let Some(block) = find_section(xml, b"dataProcessingList") {
        let mut cur = 0usize;
        let open_f = memmem::Finder::new(b"<dataProcessing ");
        let close = b"</dataProcessing>";
        let close_f = memmem::Finder::new(close);
        while let Some(p) = open_f.find(&block[cur..]) {
            let from = cur + p;
            let gt = match mc_memchr(b'>', &block[from..]) {
                Some(x) => from + x,
                None => break,
            };
            let head = &block[from..gt];
            let id = b2s(find_attr_value_in_tag(head, b"id")).unwrap_or_default();
            let end_rel = match close_f.find(&block[gt..]) {
                Some(v) => v,
                None => break,
            };
            let inner = &block[gt..gt + end_rel];
            let mut methods = Vec::new();
            let mut mcur = 0usize;
            let pm_open = memmem::Finder::new(b"<processingMethod ");
            let pm_close = memmem::Finder::new(b"</processingMethod>");
            while let Some(pp) = pm_open.find(&inner[mcur..]) {
                let sf = mcur + pp;
                let sgt = match mc_memchr(b'>', &inner[sf..]) {
                    Some(x) => sf + x,
                    None => break,
                };
                let shead = &inner[sf..sgt];
                let order = find_attr_usize(shead, b"processingMethod", b"order");
                let swref = b2s(find_attr_value_in_tag(shead, b"softwareRef"));
                let (send_rel, sinner) = if inner.get(sgt) == Some(&b'/') {
                    (1usize, &inner[sgt..sgt])
                } else if let Some(er) = pm_close.find(&inner[sgt..]) {
                    (er + b"</processingMethod>".len(), &inner[sgt..sgt + er])
                } else {
                    (0, &inner[sgt..sgt])
                };
                let cvp = collect_cv_params(sinner);
                methods.push(ProcessingMethod {
                    order,
                    software_ref: swref,
                    cv_params: cvp,
                });
                mcur = sgt + send_rel;
            }
            out.push(DataProcessing { id, methods });
            cur = gt + end_rel + close.len();
        }
    }
    out
}

fn parse_acquisition_settings_list(xml: &[u8]) -> Vec<AcquisitionSettings> {
    let mut out = Vec::new();
    if let Some(block) = find_section(xml, b"acquisitionSettingsList") {
        let mut cur = 0usize;
        let open_f = memmem::Finder::new(b"<acquisitionSettings ");
        let close = b"</acquisitionSettings>";
        let close_f = memmem::Finder::new(close);
        while let Some(p) = open_f.find(&block[cur..]) {
            let from = cur + p;
            let gt = match mc_memchr(b'>', &block[from..]) {
                Some(x) => from + x,
                None => break,
            };
            let head = &block[from..gt];
            let id = b2s(find_attr_value_in_tag(head, b"id")).unwrap_or_default();
            let inst = b2s(find_attr_value_in_tag(head, b"instrumentConfigurationRef"));
            let end_rel = match close_f.find(&block[gt..]) {
                Some(v) => v,
                None => break,
            };
            let inner = &block[gt..gt + end_rel];
            let cv_params = collect_cv_params(inner);
            out.push(AcquisitionSettings {
                id,
                instrument_ref: inst,
                cv_params,
            });
            cur = gt + end_rel + close.len();
        }
    }
    out
}

fn parse_run_header(xml: &[u8]) -> Option<Run> {
    let b = find_section(xml, b"run")?;
    let gt = mc_memchr(b'>', b).unwrap_or(0);
    let head = &b[..gt];
    let id = b2s(find_attr_value_in_tag(head, b"id")).unwrap_or_default();
    let start_time_stamp = b2s(find_attr_value_in_tag(head, b"startTimeStamp"));
    let def_icr = b2s(find_attr_value_in_tag(
        head,
        b"defaultInstrumentConfigurationRef",
    ));
    let spectrum_list_count = find_attr_usize(b, b"spectrumList", b"count");
    let chromatogram_list_count = find_attr_usize(b, b"chromatogramList", b"count");
    Some(Run {
        id,
        start_time_stamp,
        default_instrument_configuration_ref: def_icr,
        spectrum_list_count,
        chromatogram_list_count,
        spectra: Vec::new(),
        chromatograms: Vec::new(),
    })
}

fn parse_index_list_wrapper(xml: &[u8]) -> Option<IndexList> {
    let spectrum = parse_offsets_from_index(xml, b"spectrum");
    let chromatogram = parse_offsets_from_index(xml, b"chromatogram");
    if spectrum.is_empty() && chromatogram.is_empty() {
        return None;
    }
    let index_list_offset = extract_index_list_offset(xml);
    let file_checksum = if let Some((s, e)) = tag_body(xml, b"<fileChecksum>", b"</fileChecksum>") {
        b2s(Some(strip_ws(&xml[s..e])))
    } else {
        None
    };
    Some(IndexList {
        spectrum,
        chromatogram,
        index_list_offset,
        file_checksum,
    })
}

fn parse_offsets_from_index(hay: &[u8], name: &[u8]) -> Vec<IndexOffset> {
    let mut out = Vec::new();
    let mut pat = Vec::with_capacity(20 + name.len());
    pat.extend_from_slice(b"<index name=\"");
    pat.extend_from_slice(name);
    pat.push(b'"');
    if let Some(ix_start) = memmem::find(hay, &pat) {
        if let Some(ix_end_rel) = memmem::find(&hay[ix_start..], b"</index>") {
            let block = &hay[ix_start..ix_start + ix_end_rel];
            let mut cursor = 0usize;
            let off_f = memmem::Finder::new(b"<offset");
            let close_f = memmem::Finder::new(b"</offset>");
            while let Some(p) = off_f.find(&block[cursor..]) {
                let from = cursor + p;
                let g = match mc_memchr(b'>', &block[from..]) {
                    Some(rel) => from + rel + 1,
                    None => break,
                };
                let endp_rel = match close_f.find(&block[g..]) {
                    Some(v) => v,
                    None => break,
                };
                let head = &block[from..g - 1];
                let id_ref = b2s(find_attr_value_in_tag(head, b"idRef"));
                let num = strip_ws(&block[g..g + endp_rel]);
                if let Some(off) = parse_u64_ascii(num) {
                    out.push(IndexOffset {
                        id_ref,
                        offset: off,
                    });
                }
                cursor = g + endp_rel + b"</offset>".len();
            }
        }
    }
    out
}

fn parse_spectra_internal(bytes: &[u8]) -> Result<Vec<SpectrumSummary>, String> {
    let file_len = bytes.len() as u64;
    let mut cursor = Cursor::new(bytes);
    let mut scratch = Scratch {
        b64_buf: Vec::with_capacity(256),
        zlib_buf: Vec::with_capacity(256),
    };
    if let Some(offsets) = read_spectrum_offsets(&mut cursor)? {
        if file_len <= 1_073_741_824 {
            let all = bytes;
            let mut out = Vec::with_capacity(offsets.len());
            for i in 0..offsets.len() {
                let start = offsets[i] as usize;
                let end = if i + 1 < offsets.len() {
                    offsets[i + 1] as usize
                } else {
                    find_spectrum_end_in(all, start)
                        .ok_or_else(|| "no </spectrum> after last offset".to_string())?
                };
                if let Some(sum) = parse_spectrum_block(&all[start..end], &mut scratch) {
                    out.push(sum);
                }
            }
            return Ok(out);
        } else {
            let mut out = Vec::with_capacity(offsets.len());
            for i in 0..offsets.len() {
                let start = offsets[i];
                let next = if i + 1 < offsets.len() {
                    Some(offsets[i + 1])
                } else {
                    None
                };
                if let Some(sum) = read_one_spectrum_span(&mut cursor, start, next, &mut scratch)? {
                    out.push(sum);
                }
            }
            return Ok(out);
        }
    }
    cursor.set_position(0);
    linear_scan_spectra(&mut cursor, &mut scratch)
}

fn extract_index_list_offset(tail: &[u8]) -> Option<u64> {
    let tag = b"<indexListOffset>";
    let endtag = b"</indexListOffset>";
    let pos = memmem::find(tail, tag)?;
    let pos2 = memmem::find(&tail[pos + tag.len()..], endtag)?;
    let num = &tail[pos + tag.len()..pos + tag.len() + pos2];
    parse_u64_ascii(strip_ws(num))
}

fn parse_spectrum_offsets_from_index(buf: &[u8]) -> Vec<u64> {
    let mut out = Vec::new();
    if let Some(ix_start) = memmem::find(buf, br#"<index name="spectrum">"#) {
        if let Some(ix_end_rel) = memmem::find(&buf[ix_start..], b"</index>") {
            let block = &buf[ix_start..ix_start + ix_end_rel];
            let mut cursor = 0usize;
            let tag = b"<offset";
            let endtag = b"</offset>";
            let t_f = memmem::Finder::new(tag);
            let e_f = memmem::Finder::new(endtag);
            while let Some(p) = t_f.find(&block[cursor..]) {
                let from = cursor + p;
                let g = match mc_memchr(b'>', &block[from..]) {
                    Some(rel) => from + rel + 1,
                    None => break,
                };
                let endp_rel = match e_f.find(&block[g..]) {
                    Some(v) => v,
                    None => break,
                };
                let num = &block[g..g + endp_rel];
                if let Some(v) = parse_u64_ascii(strip_ws(num)) {
                    out.push(v);
                }
                cursor = g + endp_rel + endtag.len();
            }
        }
    }
    out
}

fn read_one_spectrum_span<R: Read + Seek>(
    r: &mut R,
    start: u64,
    next: Option<u64>,
    scratch: &mut Scratch,
) -> Result<Option<SpectrumSummary>, String> {
    r.seek(SeekFrom::Start(start))
        .map_err(|e| format!("seek: {e}"))?;
    if let Some(end) = next {
        let len = (end - start) as usize;
        let mut buf = vec![0u8; len];
        r.read_exact(&mut buf)
            .map_err(|e| format!("read span: {e}"))?;
        if let Some(pos) = memmem::find(&buf, b"</spectrum>") {
            buf.truncate(pos + b"</spectrum>".len());
        }
        Ok(parse_spectrum_block(&buf, scratch))
    } else {
        let mut buf = Vec::with_capacity(128 * 1024);
        let mut tmp = [0u8; 128 * 1024];
        let close = b"</spectrum>";
        let close_f = memmem::Finder::new(close);
        let mut search_from = 0usize;
        loop {
            let n = r
                .read(&mut tmp)
                .map_err(|e| format!("read tail spectrum: {e}"))?;
            if n == 0 {
                break;
            }
            buf.extend_from_slice(&tmp[..n]);
            let window_start = search_from.saturating_sub(close.len().saturating_sub(1));
            if let Some(rel) = close_f.find(&buf[window_start..]) {
                let end = window_start + rel + close.len();
                buf.truncate(end);
                break;
            }
            search_from = buf.len();
            if buf.len() > 32 * 1024 * 1024 {
                return Err("spectrum block too large?".into());
            }
        }
        Ok(parse_spectrum_block(&buf, scratch))
    }
}

fn find_spectrum_end_in(hay: &[u8], start: usize) -> Option<usize> {
    let rel = memmem::find(&hay[start..], b"</spectrum>")?;
    Some(start + rel + b"</spectrum>".len())
}

fn linear_scan_spectra<R: Read + Seek>(
    r: &mut R,
    scratch: &mut Scratch,
) -> Result<Vec<SpectrumSummary>, String> {
    r.seek(SeekFrom::Start(0))
        .map_err(|e| format!("seek: {e}"))?;
    let mut file = Vec::new();
    r.read_to_end(&mut file)
        .map_err(|e| format!("read all: {e}"))?;
    let mut out = Vec::new();
    let mut cur = 0usize;
    let open_tag = b"<spectrum ";
    let close_tag = b"</spectrum>";
    let of = memmem::Finder::new(open_tag);
    let cf = memmem::Finder::new(close_tag);
    while let Some(p) = of.find(&file[cur..]) {
        let start = cur + p;
        let end_rel = cf
            .find(&file[start..])
            .ok_or_else(|| "unterminated <spectrum>".to_string())?;
        let end = start + end_rel + close_tag.len();
        if let Some(sum) = parse_spectrum_block(&file[start..end], scratch) {
            out.push(sum);
        }
        cur = end;
    }
    Ok(out)
}

fn first_activation_name(b: &[u8]) -> Option<String> {
    if let Some((s, e)) = tag_body(b, b"<activation>", b"</activation>") {
        let a = &b[s..e];
        let mut cur = 0usize;
        let f = memmem::Finder::new(b"<cvParam");
        while let Some(p) = f.find(&a[cur..]) {
            let from = cur + p;
            let gt = mc_memchr(b'>', &a[from..]).map(|x| from + x)?;
            let head = &a[from..gt];
            if let Some(nm) = find_attr_value_in_tag(head, b"name") {
                return std::str::from_utf8(nm).ok().map(|s| s.to_string());
            }
            cur = gt + 1;
        }
    }
    None
}

fn parse_precursor_from_header(header: &[u8]) -> Option<Precursor> {
    let pblk = find_section(header, b"precursorList")?;
    let isolation_window_target_mz = find_cv_value_f64(pblk, b"isolation window target m/z");
    let isolation_window_lower_offset = find_cv_value_f64(pblk, b"isolation window lower offset");
    let isolation_window_upper_offset = find_cv_value_f64(pblk, b"isolation window upper offset");
    let selected_ion_mz = find_cv_value_f64(pblk, b"selected ion m/z");
    let activation = first_activation_name(pblk);
    if isolation_window_target_mz.is_none()
        && isolation_window_lower_offset.is_none()
        && isolation_window_upper_offset.is_none()
        && selected_ion_mz.is_none()
        && activation.is_none()
    {
        None
    } else {
        Some(Precursor {
            isolation_window_target_mz,
            isolation_window_lower_offset,
            isolation_window_upper_offset,
            selected_ion_mz,
        })
    }
}

fn parse_spectrum_block(block: &[u8], scratch: &mut Scratch) -> Option<SpectrumSummary> {
    let index = find_attr_usize(block, b"spectrum", b"index").unwrap_or(0);
    let array_len = find_attr_usize(block, b"spectrum", b"defaultArrayLength").unwrap_or(0);
    let header_end = memmem::find(block, b"<binaryDataArrayList").unwrap_or(block.len());
    let header = &block[..header_end];
    let ms_level = find_cv_value_u8(header, b"ms level");
    let polarity = if has_cv_name(header, b"positive scan") {
        Some(0)
    } else if has_cv_name(header, b"negative scan") {
        Some(1)
    } else {
        None
    };
    let spectrum_type = if has_cv_name(header, b"profile spectrum") {
        Some(0)
    } else if has_cv_name(header, b"centroid spectrum") {
        Some(1)
    } else {
        None
    };
    let total_ion_current = find_cv_value_f64(header, b"total ion current");
    let base_peak_intensity = find_cv_value_f64(header, b"base peak intensity");
    let base_peak_mz = find_cv_value_f64(header, b"base peak m/z");
    let retention_time = find_scan_start_time_min(header);
    let scan_window_lower_limit = find_cv_value_f64(header, b"scan window lower limit");
    let scan_window_upper_limit = find_cv_value_f64(header, b"scan window upper limit");
    let precursor = parse_precursor_from_header(header);
    let (mz_array_opt, intensity_array_opt) = decode_binary_arrays(block, array_len, scratch);

    let mut tic_out = total_ion_current;
    let mut bpi_out = base_peak_intensity;
    let mut bpmz_out = base_peak_mz;

    if tic_out.is_none() || bpi_out.is_none() || bpmz_out.is_none() {
        if let (Some(mz), Some(inten)) = (&mz_array_opt, &intensity_array_opt) {
            let mut tic_val: f64 = 0.0;
            let mut bpi_val: f64 = 0.0;
            let mut bpmz_val: f64 = 0.0;
            for (mzv, &intv_f32) in mz.iter().zip(inten.iter()) {
                let intv = intv_f32 as f64;
                tic_val += intv;
                if intv > bpi_val {
                    bpi_val = intv;
                    bpmz_val = *mzv;
                }
            }
            if tic_out.is_none() {
                tic_out = Some(tic_val);
            }
            if bpi_out.is_none() {
                bpi_out = Some(bpi_val);
            }
            if bpmz_out.is_none() {
                bpmz_out = Some(bpmz_val);
            }
        }
    }
    Some(SpectrumSummary {
        index,
        array_length: array_len,
        ms_level,
        polarity,
        spectrum_type,
        retention_time,
        scan_window_lower_limit,
        scan_window_upper_limit,
        total_ion_current: tic_out,
        base_peak_intensity: bpi_out,
        base_peak_mz: bpmz_out,
        mz_array: mz_array_opt,
        intensity_array: intensity_array_opt,
        precursor,
    })
}

fn has_cv_name(buf: &[u8], name: &[u8]) -> bool {
    let mut cur = 0usize;
    let f = memmem::Finder::new(b"<cvParam");
    while let Some(p) = f.find(&buf[cur..]) {
        let from = cur + p;
        if let Some(gt_rel) = mc_memchr(b'>', &buf[from..]) {
            let gt = from + gt_rel;
            let head = &buf[from..gt];
            if let Some(v) = find_attr_value_in_tag(head, b"name") {
                if v == name {
                    return true;
                }
            }
            cur = gt + 1;
        } else {
            return false;
        }
    }
    false
}

fn find_cv_value_f64(buf: &[u8], name: &[u8]) -> Option<f64> {
    find_cv_value(buf, name).and_then(|s| str::from_utf8(s).ok()?.parse().ok())
}

fn find_cv_value_u8(buf: &[u8], name: &[u8]) -> Option<u8> {
    find_cv_value(buf, name).and_then(|s| str::from_utf8(s).ok()?.parse().ok())
}

fn find_cv_value<'a>(buf: &'a [u8], name: &[u8]) -> Option<&'a [u8]> {
    let mut cur = 0usize;
    let f = memmem::Finder::new(b"<cvParam");
    while let Some(p) = f.find(&buf[cur..]) {
        let from = cur + p;
        let gt = mc_memchr(b'>', &buf[from..]).map(|x| from + x)?;
        let head = &buf[from..gt];
        if let Some(nm) = find_attr_value_in_tag(head, b"name") {
            if nm == name {
                return find_attr_value_in_tag(head, b"value");
            }
        }
        cur = gt + 1;
    }
    None
}

fn find_scan_start_time_min(buf: &[u8]) -> Option<f64> {
    let mut cur = 0usize;
    let f = memmem::Finder::new(b"<cvParam");
    while let Some(p) = f.find(&buf[cur..]) {
        let from = cur + p;
        let gt = mc_memchr(b'>', &buf[from..]).map(|x| from + x)?;
        let head = &buf[from..gt];
        if let Some(nm) = find_attr_value_in_tag(head, b"name") {
            if nm == b"scan start time" {
                let val = find_attr_value_in_tag(head, b"value")?;
                let mut v: f64 = str::from_utf8(val).ok()?.parse().ok()?;
                if let Some(unit) = find_attr_value_in_tag(head, b"unitName") {
                    if unit == b"second" {
                        v /= 60.0;
                    }
                }
                return Some(v);
            }
        }
        cur = gt + 1;
    }
    None
}

fn bda_flags(b: &[u8]) -> (bool, bool, bool, bool, bool, bool) {
    let stop = memmem::find(b, b"<binary>").unwrap_or(b.len());
    let head = &b[..stop];
    let mut kind_mz = false;
    let mut kind_int = false;
    let mut is_zlib = false;
    let mut is_f64 = false;
    let mut is_f32 = false;
    let mut little = true;
    let mut cur = 0usize;
    let f = memmem::Finder::new(b"<cvParam");
    while let Some(p) = f.find(&head[cur..]) {
        let from = cur + p;
        if let Some(gt_rel) = mc_memchr(b'>', &head[from..]) {
            let gt = from + gt_rel;
            let tag_head = &head[from..gt];
            if let Some(nm) = find_attr_value_in_tag(tag_head, b"name") {
                match nm {
                    b"m/z array" => kind_mz = true,
                    b"intensity array" => kind_int = true,
                    b"zlib compression" => is_zlib = true,
                    b"64-bit float" => is_f64 = true,
                    b"32-bit float" => is_f32 = true,
                    b"little endian" => little = true,
                    b"big endian" => little = false,
                    _ => {}
                }
            }
            cur = gt + 1;
        } else {
            break;
        }
    }
    (kind_mz, kind_int, is_zlib, is_f64, is_f32, little)
}

fn decode_binary_arrays(
    block: &[u8],
    expected_len: usize,
    scratch: &mut Scratch,
) -> (Option<Vec<f64>>, Option<Vec<f32>>) {
    let mut mz: Option<Vec<f64>> = None;
    let mut inten: Option<Vec<f32>> = None;

    let mut cur = 0usize;
    let bda_open = memmem::Finder::new(b"<binaryDataArray");
    let bda_close = b"</binaryDataArray>";
    let bda_close_f = memmem::Finder::new(bda_close);

    while let Some(p) = bda_open.find(&block[cur..]) {
        let start = cur + p;
        let end_rel = match bda_close_f.find(&block[start..]) {
            Some(v) => v,
            None => break,
        };
        let b = &block[start..start + end_rel];

        let (kind_mz, kind_int, is_zlib, is_f64, is_f32, little) = bda_flags(b);

        if let Some((bs, be)) = tag_body(b, b"<binary>", b"</binary>") {
            strip_b64_ws_into(&b[bs..be], &mut scratch.zlib_buf);
            let cap = (scratch.zlib_buf.len() / 4) * 3 + 3;
            scratch.b64_buf.resize(cap, 0);
            let written = match STANDARD.decode_slice(&scratch.zlib_buf, &mut scratch.b64_buf) {
                Ok(n) => n,
                Err(_) => {
                    cur = start + end_rel + bda_close.len();
                    continue;
                }
            };
            scratch.b64_buf.truncate(written);

            let bytes: &[u8] = if is_zlib {
                scratch.zlib_buf.clear();
                match decompress_to_vec_zlib(&scratch.b64_buf) {
                    Ok(v) => {
                        scratch.zlib_buf = v;
                        &scratch.zlib_buf
                    }
                    Err(_) => {
                        cur = start + end_rel + bda_close.len();
                        continue;
                    }
                }
            } else {
                &scratch.b64_buf
            };

            let want = if expected_len > 0 {
                expected_len
            } else if is_f64 {
                bytes.len() / 8
            } else {
                bytes.len() / 4
            };

            if kind_mz {
                let vals = if is_f64 {
                    bytes_to_f64_exact_into(bytes, little, want)
                } else if is_f32 {
                    bytes_to_f32_as_f64_exact_into(bytes, little, want)
                } else {
                    Vec::new()
                };
                mz = Some(vals);
            } else if kind_int {
                let vals = if is_f32 {
                    bytes_to_f32_exact_into(bytes, little, want)
                } else if is_f64 {
                    bytes_to_f64_as_f32_exact_into(bytes, little, want)
                } else {
                    Vec::new()
                };
                inten = Some(vals);
            }
        }

        cur = start + end_rel + bda_close.len();
    }

    (mz, inten)
}

fn bda_flags_chrom(b: &[u8]) -> (bool, bool, bool, bool, bool, bool) {
    let stop = memmem::find(b, b"<binary>").unwrap_or(b.len());
    let head = &b[..stop];
    let mut kind_time = false;
    let mut kind_int = false;
    let mut is_zlib = false;
    let mut is_f64 = false;
    let mut is_f32 = false;
    let mut little = true;
    let mut cur = 0usize;
    let f = memmem::Finder::new(b"<cvParam");
    while let Some(p) = f.find(&head[cur..]) {
        let from = cur + p;
        if let Some(gt_rel) = mc_memchr(b'>', &head[from..]) {
            let gt = from + gt_rel;
            let tag_head = &head[from..gt];
            if let Some(nm) = find_attr_value_in_tag(tag_head, b"name") {
                match nm {
                    b"time array" => kind_time = true,
                    b"intensity array" => kind_int = true,
                    b"zlib compression" => is_zlib = true,
                    b"64-bit float" => is_f64 = true,
                    b"32-bit float" => is_f32 = true,
                    b"little endian" => little = true,
                    b"big endian" => little = false,
                    _ => {}
                }
            }
            cur = gt + 1;
        } else {
            break;
        }
    }
    (kind_time, kind_int, is_zlib, is_f64, is_f32, little)
}

fn decode_chrom_binary_arrays(
    block: &[u8],
    expected_len: usize,
    scratch: &mut Scratch,
) -> (Option<Vec<f64>>, Option<Vec<f32>>) {
    let mut time_arr: Option<Vec<f64>> = None;
    let mut intensity_arr: Option<Vec<f32>> = None;

    let mut cur = 0usize;
    let bda_open = memmem::Finder::new(b"<binaryDataArray");
    let bda_close = b"</binaryDataArray>";
    let bda_close_f = memmem::Finder::new(bda_close);

    while let Some(p) = bda_open.find(&block[cur..]) {
        let start = cur + p;
        let end_rel = match bda_close_f.find(&block[start..]) {
            Some(v) => v,
            None => break,
        };
        let b = &block[start..start + end_rel];

        let (kind_time, kind_int, is_zlib, is_f64, is_f32, little) = bda_flags_chrom(b);

        if let Some((bs, be)) = tag_body(b, b"<binary>", b"</binary>") {
            strip_b64_ws_into(&b[bs..be], &mut scratch.zlib_buf);
            let cap = (scratch.zlib_buf.len() / 4) * 3 + 3;
            scratch.b64_buf.resize(cap, 0);
            let written = match STANDARD.decode_slice(&scratch.zlib_buf, &mut scratch.b64_buf) {
                Ok(n) => n,
                Err(_) => {
                    cur = start + end_rel + bda_close.len();
                    continue;
                }
            };
            scratch.b64_buf.truncate(written);

            let bytes: &[u8] = if is_zlib {
                scratch.zlib_buf.clear();
                match decompress_to_vec_zlib(&scratch.b64_buf) {
                    Ok(v) => {
                        scratch.zlib_buf = v;
                        &scratch.zlib_buf
                    }
                    Err(_) => {
                        cur = start + end_rel + bda_close.len();
                        continue;
                    }
                }
            } else {
                &scratch.b64_buf
            };

            let want = if expected_len > 0 {
                expected_len
            } else if is_f64 {
                bytes.len() / 8
            } else {
                bytes.len() / 4
            };

            if kind_time {
                let vals = if is_f64 {
                    bytes_to_f64_exact_into(bytes, little, want)
                } else if is_f32 {
                    bytes_to_f32_as_f64_exact_into(bytes, little, want)
                } else {
                    Vec::new()
                };
                time_arr = Some(vals);
            } else if kind_int {
                let vals = if is_f32 {
                    bytes_to_f32_exact_into(bytes, little, want)
                } else if is_f64 {
                    bytes_to_f64_as_f32_exact_into(bytes, little, want)
                } else {
                    Vec::new()
                };
                intensity_arr = Some(vals);
            }
        }

        cur = start + end_rel + bda_close.len();
    }

    (time_arr, intensity_arr)
}

fn parse_chromatograms_linear(xml: &[u8]) -> Result<Vec<ChromatogramSummary>, String> {
    let mut scratch = Scratch {
        b64_buf: Vec::with_capacity(256),
        zlib_buf: Vec::with_capacity(256),
    };
    let mut out = Vec::new();
    let mut cur = 0usize;
    const OPEN: &[u8] = b"<chromatogram ";
    const CLOSE: &[u8] = b"</chromatogram>";
    let of = memmem::Finder::new(OPEN);
    let cf = memmem::Finder::new(CLOSE);
    while let Some(p) = of.find(&xml[cur..]) {
        let start = cur + p;
        let end_rel = cf
            .find(&xml[start..])
            .ok_or_else(|| "unterminated <chromatogram>".to_string())?;
        let end = start + end_rel + CLOSE.len();
        if let Some(ch) = parse_chromatogram_block(&xml[start..end], &mut scratch) {
            out.push(ch);
        }
        cur = end;
        if out.len() > 2_000_000 {
            return Err("too many chromatograms?".into());
        }
    }
    Ok(out)
}

fn parse_chromatogram_block(block: &[u8], scratch: &mut Scratch) -> Option<ChromatogramSummary> {
    let index = find_attr_usize(block, b"chromatogram", b"index").unwrap_or(0);
    let array_length = find_attr_usize(block, b"chromatogram", b"defaultArrayLength").unwrap_or(0);
    let (time_array, intensity_array) = decode_chrom_binary_arrays(block, array_length, scratch);
    let id = find_attr_string(block, b"chromatogram", b"id").unwrap_or_default();
    Some(ChromatogramSummary {
        index,
        array_length,
        time_array,
        intensity_array,
        id,
    })
}

fn read_spectrum_offsets<R: Read + Seek>(r: &mut R) -> Result<Option<Vec<u64>>, String> {
    const TAIL: u64 = 64 * 1024;
    let end = r
        .seek(SeekFrom::End(0))
        .map_err(|e| format!("seek end: {e}"))?;
    let start = end.saturating_sub(TAIL);
    r.seek(SeekFrom::Start(start))
        .map_err(|e| format!("seek tail: {e}"))?;
    let mut tail = Vec::with_capacity((end - start) as usize);
    r.take(end - start)
        .read_to_end(&mut tail)
        .map_err(|e| format!("read tail: {e}"))?;
    if let Some(off) = extract_index_list_offset(&tail) {
        r.seek(SeekFrom::Start(off))
            .map_err(|e| format!("seek indexList: {e}"))?;
        let mut buf = vec![0u8; (end - off).min(1_000_000).max(4096) as usize];
        let n = r
            .read(&mut buf)
            .map_err(|e| format!("read indexList: {e}"))?;
        buf.truncate(n);
        let offs = parse_spectrum_offsets_from_index(&buf);
        return Ok(Some(offs));
    }
    Ok(None)
}

fn find_attr_usize(buf: &[u8], tag: &[u8], attr: &[u8]) -> Option<usize> {
    find_attr_ascii(buf, tag, attr).and_then(|s| str::from_utf8(s).ok()?.parse().ok())
}

fn find_attr_string(buf: &[u8], tag: &[u8], attr: &[u8]) -> Option<String> {
    let b = find_attr_ascii(buf, tag, attr)?;
    std::str::from_utf8(b).ok().map(|s| s.to_owned())
}

fn find_attr_ascii<'a>(buf: &'a [u8], tag: &[u8], attr: &[u8]) -> Option<&'a [u8]> {
    let mut pat = Vec::with_capacity(1 + tag.len());
    pat.push(b'<');
    pat.extend_from_slice(tag);
    let start = memmem::find(buf, &pat)?;
    let gt = mc_memchr(b'>', &buf[start..]).map(|x| start + x)?;
    let head = &buf[start..gt];
    find_attr_value_in_tag(head, attr)
}

fn find_attr_value_in_tag<'a>(head: &'a [u8], attr: &[u8]) -> Option<&'a [u8]> {
    let mut pat = Vec::with_capacity(attr.len() + 1);
    pat.extend_from_slice(attr);
    pat.push(b'=');
    let p = memmem::find(head, &pat)?;
    let q = p + pat.len();
    let quote = *head.get(q)?;
    if quote != b'"' && quote != b'\'' {
        return None;
    }
    let rest = &head[q + 1..];
    let end = mc_memchr(quote, rest)?;
    Some(&rest[..end])
}

fn strip_b64_ws_into(src: &[u8], dst: &mut Vec<u8>) {
    dst.clear();
    dst.extend(
        src.iter()
            .copied()
            .filter(|b| !matches!(b, b' ' | b'\n' | b'\r' | b'\t')),
    );
}

fn tag_body(hay: &[u8], open: &[u8], close: &[u8]) -> Option<(usize, usize)> {
    let s = memmem::find(hay, open)?;
    let e_rel = memmem::find(&hay[s + open.len()..], close)?;
    Some((s + open.len(), s + open.len() + e_rel))
}

fn is_ws(b: u8) -> bool {
    matches!(b, b' ' | b'\n' | b'\r' | b'\t')
}

fn bytes_to_f32_exact_into(b: &[u8], little: bool, want: usize) -> Vec<f32> {
    let len = want.min(b.len() / 4);
    let mut out = Vec::with_capacity(len);
    let words = &b[..len * 4];
    for c in words.chunks_exact(4) {
        let bits = if little {
            u32::from_le_bytes([c[0], c[1], c[2], c[3]])
        } else {
            u32::from_be_bytes([c[0], c[1], c[2], c[3]])
        };
        out.push(f32::from_bits(bits));
    }
    out
}

fn bytes_to_f64_as_f32_exact_into(b: &[u8], little: bool, want: usize) -> Vec<f32> {
    let len = want.min(b.len() / 8);
    let mut out = Vec::with_capacity(len);
    let bytes = &b[..len * 8];
    for c in bytes.chunks_exact(8) {
        let bits = if little {
            u64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]])
        } else {
            u64::from_be_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]])
        };
        out.push(f64::from_bits(bits) as f32);
    }
    out
}

fn bytes_to_f64_exact_into(b: &[u8], little: bool, want: usize) -> Vec<f64> {
    let len = want.min(b.len() / 8);
    let mut out = Vec::with_capacity(len);
    let bytes = &b[..len * 8];
    if little {
        for c in bytes.chunks_exact(8) {
            let bits = u64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]);
            out.push(f64::from_bits(bits));
        }
    } else {
        for c in bytes.chunks_exact(8) {
            let bits = u64::from_be_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]);
            out.push(f64::from_bits(bits));
        }
    }
    out
}

fn bytes_to_f32_as_f64_exact_into(b: &[u8], little: bool, want: usize) -> Vec<f64> {
    let len = want.min(b.len() / 4);
    let mut out = Vec::with_capacity(len);
    let words = &b[..len * 4];
    for c in words.chunks_exact(4) {
        let bits = if little {
            u32::from_le_bytes([c[0], c[1], c[2], c[3]])
        } else {
            u32::from_be_bytes([c[0], c[1], c[2], c[3]])
        };
        out.push(f32::from_bits(bits) as f64);
    }
    out
}

fn strip_ws(s: &[u8]) -> &[u8] {
    let mut a = 0;
    let mut b = s.len();
    while a < b && is_ws(s[a]) {
        a += 1;
    }
    while b > a && is_ws(s[b - 1]) {
        b -= 1;
    }
    &s[a..b]
}

fn parse_u64_ascii(s: &[u8]) -> Option<u64> {
    let t = strip_ws(s);
    if t.is_empty() {
        return None;
    }
    let mut v: u64 = 0;
    for &c in t {
        if c < b'0' || c > b'9' {
            return None;
        }
        v = v.checked_mul(10)?.checked_add((c - b'0') as u64)?;
    }
    Some(v)
}

fn find_section<'a>(xml: &'a [u8], name: &[u8]) -> Option<&'a [u8]> {
    let mut open = Vec::with_capacity(name.len() + 1);
    open.push(b'<');
    open.extend_from_slice(name);
    let s = memmem::find(xml, &open)?;
    let mut close = Vec::with_capacity(name.len() + 3);
    close.extend_from_slice(b"</");
    close.extend_from_slice(name);
    close.push(b'>');
    let e_rel = memmem::find(&xml[s..], &close)?;
    Some(&xml[s..s + e_rel + close.len()])
}

fn collect_cv_params(block: &[u8]) -> Vec<CvPair> {
    let mut out = Vec::new();
    let mut cur = 0usize;
    let f = memmem::Finder::new(b"<cvParam ");
    while let Some(p) = f.find(&block[cur..]) {
        let from = cur + p;
        let gt = match mc_memchr(b'>', &block[from..]) {
            Some(x) => from + x,
            None => break,
        };
        let head = &block[from..gt];
        let name = b2s(find_attr_value_in_tag(head, b"name")).unwrap_or_default();
        let value = b2s(find_attr_value_in_tag(head, b"value"));
        let accession = b2s(find_attr_value_in_tag(head, b"accession"));
        let cv_ref = b2s(find_attr_value_in_tag(head, b"cvRef"));
        let unit_name = b2s(find_attr_value_in_tag(head, b"unitName"));
        out.push(CvPair {
            name,
            value,
            accession,
            cv_ref,
            unit_name,
        });
        cur = gt + 1;
    }
    out
}

fn b2s(opt: Option<&[u8]>) -> Option<String> {
    opt.and_then(|b| std::str::from_utf8(b).ok())
        .map(|s| s.to_owned())
}
